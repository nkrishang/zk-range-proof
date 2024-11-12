from py_ecc.bn128 import multiply, eq, neg, curve_order as p
from utils import *

class Verifier:
    """Verifies a zero knowledge proof (bulletproof) that v < 2^n"""

    def __init__(self, A, S, V, G_vec, H_vec, Q, B):
        # Verifier generated randomness
        y = random_field_element()
        z = random_field_element()
        
        self.y = y
        self.z = z
        
        # Public factors
        n = len(G_vec)

        vec_1n = get_vec_1n(n)
        vec_2n = get_vec_2n(n)
        vec_yn = get_vec_yn(y, n)
        vec_yn_inv = get_vec_yn_inv(y, n)

        self.vec_1n = vec_1n
        self.vec_2n = vec_2n
        self.vec_yn = vec_yn

        # Prover commitments (A, S, V)
        self.A = A
        self.S = S
        self.V = V

        # Random EC points
        H_y_inv = points_vec_mul(H_vec, vec_yn_inv)

        self.G_vec = G_vec
        self.Gprime_vec = G_vec
        self.H_vec = H_y_inv
        self.Hprime_vec = H_y_inv
        
        self.Q = Q
        self.B = B

    # Store commitments (T1, T2) and generate randomness u
    def commit_T1_T2(self, T1, T2):
        self.T1 = T1
        self.T2 = T2

        self.u = random_field_element()
    
    # Store commitments (C, Q.tu, B.pi_lr, B.pi_t) and calculate commitment P
    def commit_C_tu_pilr_pit(self, C, Qtu, Bpi_lr, Bpi_t):
        P = add_points(C, Qtu)

        self.C = C
        self.Qtu = Qtu
        self.P = P

        self.Bpi_lr = Bpi_lr
        self.Bpi_t = Bpi_t
    
    # Compute P' using off-diagonal commitments and fold vars
    def commit_L_R(self, L, R):
        
        G_vec = self.Gprime_vec
        H_vec = self.Hprime_vec
        P = self.P

        u = random_field_element()

        Pprime = add_points(multiply(L, pow(u, 2, p)), P, multiply(R, pow(u, -2, p)))
        
        self.interactive_proof_randomness = u
        
        self.Gprime_vec = fold_points(G_vec, pow(u, -1, p))
        self.Hprime_vec = fold_points(H_vec, pow(u, -1, p))

        self.P = Pprime

    def verify(self, a, b, b_inv):

        z = self.z
        u = self.u

        Q = self.Q
        Gprime_vec = self.Gprime_vec
        Hprime_vec = self.Hprime_vec
        G_vec = self.G_vec
        H_vec = self.H_vec

        A = self.A
        S = self.S
        V = self.V
        T1 = self.T1
        T2 = self.T2

        C = self.C
        Qtu = self.Qtu
        Bpi_lr = self.Bpi_lr
        Bpi_t = self.Bpi_t

        P = self.P

        # Public factors
        vec_1n = self.vec_1n
        vec_2n = self.vec_2n
        vec_yn = self.vec_yn

        # -z.1n
        j = vector_commit(G_vec, mod_scalar_mul(vec_1n, -z, p))
        # z.yn + z^2.2n
        k = vector_commit(H_vec, mod_vec_add(mod_scalar_mul(vec_yn, z, p), mod_scalar_mul(vec_2n, pow(z, 2, p), p), p))
        
        delta = (
            # (z - z^2)⟨1^n, y^n⟩
            mod_inner([((z - pow(z, 2, p)) % p)], [mod_inner(vec_1n, vec_yn, p)], p) - 
            # z^3⟨1^n, 2^n⟩
            mod_inner([pow(z, 3, p)], [mod_inner(vec_1n, vec_2n, p)], p)
        ) % p

        # Verification check 1: P == aG + bH + abQ
        check = eq(P, add_points(vector_commit(Gprime_vec, a), vector_commit(Hprime_vec, b), multiply(Q, mod_inner(a, b_inv, p))))
        verification = check

        print("1: ", check)

        # Verification check 2: C == A + Su + j + k - B*pi_lr
        check = eq(C, add_points(A, multiply(S, u), j, k, neg(Bpi_lr)))
        verification = verification and check

        print("2: ", check)

        # Verification check 3: Q.tu + B.pi_t == V.z^2 + Q.delta + T1.u + T2.u^2
        check = eq(add_points(Qtu, Bpi_t), add_points(multiply(V, pow(z, 2, p)), multiply(Q, delta), multiply(T1, u), multiply(T2, pow(u, 2, p))))
        verification = verification and check

        print("3: ", check)
        
        return verification

    def get_y(self):
        return self.y
    
    def get_z(self):
        return self.z
    
    def get_u(self):
        return self.u

    def get_interactive_proof_randomness(self):
        return self.interactive_proof_randomness