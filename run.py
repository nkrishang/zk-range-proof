from prover import Prover
from verifier import Verifier
from utils import generateRandomECPointVec, next_power_of_2, log2

v = 50000
n = 23

assert n > 0

# Generate random EC point basis vectors G,H and a random EC point Q
G_vec = generateRandomECPointVec(n)
H_vec = generateRandomECPointVec(n)
Q = generateRandomECPointVec(1)[0]
B = generateRandomECPointVec(1)[0]

# Setup prover
prover = Prover(v, n, G_vec, H_vec, Q, B)

# Generate commitments A, S, V
A, S, V = prover.generate_commitments_A_S_V()

# Instaniate verifier with commitments A, S, V.
verifier = Verifier(A, S, V, G_vec, H_vec, Q, B)

# Compute lu, ru, tu polynomial coefficients
prover.compute_polynomial_coefficients(verifier.get_y(), verifier.get_z())

# Send verifier-generated randomness y,z to prover and generate commitments T1, T2
T1, T2 = prover.generate_commitments_T1_T2()

# Send commitments T1, T2 to verifier
verifier.commit_T1_T2(T1, T2)

# Evaluate polynomials at verifier generated randomness u
prover.evaluate_polynomials(verifier.get_u())

# Send verifier-generated randomness u to prover and generate commitments (C, Qtu, Bpi_lr, Bpi_t)
C, Qtu, Bpi_lr, Bpi_t = prover.generate_commitments_C_tu_pilr_pit(verifier.get_z(), verifier.get_u())

# Commit (C, Qtu, Bpi_lr, Bpi_t) to verifier
verifier.commit_C_tu_pilr_pit(C, Qtu, Bpi_lr, Bpi_t)

# Get witness a,b, b_inv
a, b, b_inv = prover.get_witness()
assert len(a) == len(b) == len(b_inv) == n

# Verification
verification = False

if n == 1:
    verification = verifier.verify(a, b, b_inv)
else:
    # Interactive proof
    u = 1
    ctr = log2(next_power_of_2(n))

    while ctr > 0:
        # Prover commits off-diagonals L, R to verifier
        L, R = prover.generate_off_diagonals_L_R()

        # Commit L,R to verifier
        verifier.commit_L_R(L, R)

        # Verifier returns randomness u
        u = verifier.get_interactive_proof_randomness()

        # Update witness
        prover.update_witness(u)

        ctr -= 1
    
    # Get witness a,b, b_inv
    a, b, b_inv = prover.get_witness()
    assert len(a) == len(b) == len(b_inv) == 1

    verification = verifier.verify(a, b, b_inv)

assert verification, "range proof failed"

print("accepted")

