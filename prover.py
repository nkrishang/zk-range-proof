from py_ecc.bn128 import curve_order as p
from utils import *

class Prover:
    """Generates a zero knowledge proof (bulletproof) that v < 2^n"""

    def __init__(self, v, n, G_vec, H_vec, Q, B):
        # Public factors
        vec_1n = get_vec_1n(n)
        vec_2n = get_vec_2n(n)

        self.vec_1n = vec_1n
        self.vec_2n = vec_2n

        # Inputs
        self.v = v
        self.n = n

        # Constant term scalars
        aL = get_binary_as_array(v)
        assert len(aL) <= n, "Binary(v) has more bits than n. So, (v < 2^n) cannot be proved since it's incorrect."
        aL = aL + [0]*(n - len(aL))

        aR = mod_vec_add(aL, mod_scalar_mul(vec_1n, -1, p), p)

        self.aL = aL
        self.aR = aR

        # Linear term scalars
        sL = random_scalar_vector(n)
        sR = random_scalar_vector(n)

        self.sL = sL
        self.sR = sR

        # Random EC points
        self.G_vec = G_vec
        self.H_vec = H_vec
        self.Q = Q
        self.B = B

        # Prover blinding terms
        alpha = random_field_element()
        beta = random_field_element()
        gamma = random_field_element()
        tau_1 = random_field_element()
        tau_2 = random_field_element()

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.tau_1 = tau_1
        self.tau_2 = tau_2

        # Commitments A, S, V
        A = add_points(vector_commit(G_vec, aL), vector_commit(H_vec, aR), multiply(B, alpha))
        S = add_points(vector_commit(G_vec, sL), vector_commit(H_vec, sR), multiply(B, beta))
        V = add_points(multiply(Q, v), multiply(B, gamma))

        self.A = A
        self.S = S
        self.V = V
    
    # Generates commitments (A, S, V)
    def generate_commitments_A_S_V(self):        
        return (self.A, self.S, self.V)

    # Computes lu, ru, tu polynomial coefficients
    def compute_polynomial_coefficients(self, y, z):
        aL = self.aL
        aR = self.aR
        sL = self.sL
        sR = self.sR
        n = self.n

        vec_1n = self.vec_1n
        vec_2n = self.vec_2n
        vec_yn = get_vec_yn(y, n)
        vec_yn_inv = get_vec_yn_inv(y, n)

        l_u_term1 = mod_vec_add(aL, mod_scalar_mul(vec_1n, -z, p), p)
        l_u_term2 = sL

        r_u_term1 = mod_vec_add(mod_vec_add(mod_vec_mul(vec_yn, aR, p), mod_scalar_mul(vec_yn, z, p), p), mod_scalar_mul(vec_2n, pow(z, 2, p), p), p)
        r_u_term2 = mod_vec_mul(vec_yn, sR, p)

        tu_term1 = mod_inner(l_u_term1, r_u_term1, p)
        tu_term2 = (mod_inner(l_u_term1, r_u_term2, p) + mod_inner(l_u_term2, r_u_term1, p)) % p
        tu_term3 = mod_inner(l_u_term2, r_u_term2, p)

        self.l_u_term1 = l_u_term1
        self.l_u_term2 = l_u_term2
        self.r_u_term1 = r_u_term1
        self.r_u_term2 = r_u_term2
        self.tu_term1 = tu_term1
        self.tu_term2 = tu_term2
        self.tu_term3 = tu_term3

        self.vec_yn = vec_yn
        self.vec_yn_inv = vec_yn_inv
        
    
    # Generates commitments (T1, T2)
    def generate_commitments_T1_T2(self):

        Q = self.Q
        B = self.B
        
        tau_1 = self.tau_1
        tau_2 = self.tau_2

        tu_term2 = self.tu_term2
        tu_term3 = self.tu_term3

        T1 = add_points(multiply(Q, tu_term2), multiply(B, tau_1))
        T2 = add_points(multiply(Q, tu_term3), multiply(B, tau_2))

        return (T1, T2)

    # Evaluate polynomials at verifier generated random point u
    def evaluate_polynomials(self, u):

        l_u_term1 = self.l_u_term1
        l_u_term2 = self.l_u_term2

        r_u_term1 = self.r_u_term1
        r_u_term2 = self.r_u_term2

        tu_term1 = self.tu_term1
        tu_term2 = self.tu_term2
        tu_term3 = self.tu_term3

        lu = mod_vec_add(l_u_term1, mod_scalar_mul(l_u_term2, u, p), p)
        ru = mod_vec_add(r_u_term1, mod_scalar_mul(r_u_term2, u, p), p)
        tu = (tu_term1 + tu_term2 * u + tu_term3 * pow(u, 2, p)) % p

        self.lu = lu
        self.ru = ru
        self.tu = tu

        self.a = lu
        self.b = ru
        self.b_inv = ru


    # Generates commitments (C, tu, pi_lr, pi_t)
    def generate_commitments_C_tu_pilr_pit(self, z, u):
        Q = self.Q
        B = self.B
        G_vec = self.G_vec
        H_vec = self.H_vec

        alpha = self.alpha
        beta = self.beta
        gamma = self.gamma
        tau_1 = self.tau_1
        tau_2 = self.tau_2

        lu = self.lu
        ru = self.ru
        tu = self.tu

        vec_yn_inv = self.vec_yn_inv
        
        # Compute new basis vector
        H_y_inv = points_vec_mul(H_vec, vec_yn_inv)

        C = add_points(vector_commit(G_vec, lu), vector_commit(H_y_inv, ru))
        pi_lr = alpha + (beta * u % p) % p
        pi_t = ((gamma * pow(z, 2, p) % p) + (tau_1 * u % p) + (tau_2 * pow(u, 2, p) % p)) % p

        # Use new basis vector for the rest of the prover computations
        self.H_vec = H_y_inv

        return (C, multiply(Q, tu), multiply(B, pi_lr), multiply(B, pi_t))
    
    # Generates L and R off-diagonal commitments
    def generate_off_diagonals_L_R(self):
        Q = self.Q
        G_vec = self.G_vec
        H_vec = self.H_vec

        a = self.a
        b = self.b
        b_inv = self.b_inv

        # Compute L and R    
        L_a, R_a = compute_secondary_diagonal(G_vec, a)
        L_b, R_b = compute_secondary_diagonal(H_vec, b)
        L_q, R_q = compute_secondary_diagonal_scalar(b_inv, a)

        L = add_points(L_a, L_b, multiply(Q, L_q))
        R = add_points(R_a, R_b, multiply(Q, R_q))

        self.L = L
        self.R = R
        
        return L, R
    
    # Gets witness to send to verifier for proof verification
    def get_witness(self):
        return self.a, self.b, self.b_inv
    
    # Update (fold) witness for log-size proof process
    def update_witness(self, u):
        # Fold
        aprime = fold(self.a, u)
        bprime = fold(self.b, u)
        bprime_inv = fold(self.b_inv, pow(u, -1, p))

        Gprime = fold_points(self.G_vec, pow(u, -1, p))
        Hprime = fold_points(self.H_vec, pow(u, -1, p))

        self.a = aprime
        self.b = bprime
        self.b_inv = bprime_inv

        self.G_vec = Gprime
        self.H_vec = Hprime
