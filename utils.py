import random
import math
from libnum import has_sqrtmod_prime_power, sqrtmod_prime_power
from py_ecc.bn128 import is_on_curve, FQ, multiply, add, Z1, curve_order as p
from py_ecc.fields import field_properties
from functools import reduce


field_mod = field_properties["bn128"]["field_modulus"]

####################### EC & VECTOR OPERATIONS ###########################

# Adds up multiple EC points
def add_points(*points):
    return reduce(add, points, Z1)

# Inner product of an EC point vector and a scalar vector
def vector_commit(points, scalars):
    return reduce(add, [multiply(P, i) for P, i in zip(points, scalars)], Z1)

# Element-wise scalar and point vector multiplication mod curve order p
def points_vec_mul(points, scalars):
    return [multiply(x, y) for x, y in zip(points, scalars)]

# Inner product of two scalar vectors mod curve order p
def mod_inner(a, b, p):
    return sum((x * y) % p for x, y in zip(a, b)) % p

# Scalar multiplication of a scalar vector by some factor mod curve order p
def mod_scalar_mul(arr, scalar, p):
    return [(x * scalar) % p for x in arr]

# Scalar vector addition mod curve order p
def mod_vec_add(a, b, p):
    return [(x + y) % p for x, y in zip(a, b)]

# Element-wise scalar vector multiplication mod curve order p
def mod_vec_mul(a, b, p):
    return [(x * y) % p for x, y in zip(a, b)]

# Returns a random element from the scalar field of the bn128 elliptic curve.
def random_field_element():
    return random.randint(0, p)

# Returns a random scalar vector of length n
def random_scalar_vector(n):
    return [random_field_element() for _ in range(n)]

# Returns an array of the binary representation (little endian) of n
def get_binary_as_array(n):
    if n < 0:
        raise ValueError("Input must be a positive integer")
    if n == 0:
        return [0]
        
    binary = []
    while n > 0:
        binary.append(n & 1)  # Get least significant bit
        n >>= 1              # Right shift by 1 (divide by 2)
    return binary

# Generates a random EC point vector of length n
def generateRandomECPointVec(n):
    b = 3 # for bn128, y^2 = x^3 + 3

    x = random_field_element()

    entropy = 0
    vector_basis = []

    while len(vector_basis) < n:
        while not has_sqrtmod_prime_power((x**3 + b) % field_mod, field_mod, 1):
            # increment x, so hopefully we are on the curve
            x = (x + 1) % field_mod
            entropy = entropy + 1

        # pick the upper or lower point depending on if entropy is even or odd
        y = list(sqrtmod_prime_power((x**3 + b) % field_mod, field_mod, 1))[entropy & 1 == 0]
        point = (FQ(x), FQ(y))
        assert is_on_curve(point, b), "sanity check"
        vector_basis.append(point)

        # new x value
        x = random_field_element()
    
    return vector_basis

# Fold a scalar vector
def fold(scalar_vec, u):
    
    length = len(scalar_vec)
    n = next_power_of_2(length)

    if not is_power_of_2(length):
        scalar_vec = scalar_vec + [0]*(n-length)

    i = 0
    vec = []
    while i < len(scalar_vec):
        vec.append((mod_inner([scalar_vec[i]], [u], p) + mod_inner([scalar_vec[i+1]], [pow(u, -1, p)], p)) % p)
        i += 2
    return vec

# Fold an EC points vector
def fold_points(point_vec, u):
    length = len(point_vec)
    n = next_power_of_2(length)

    if not is_power_of_2(length):
        point_vec = point_vec + [Z1]*(n-length)

    i = 0
    vec = []
    while i < len(point_vec):
        vec.append(add_points(multiply(point_vec[i], u), multiply(point_vec[i+1], pow(u, -1, p))))
        i += 2
    return vec

# Compute the secondary diagonal L,R for a scalar vector and EC point vector
def compute_secondary_diagonal(G_vec, a):
    R = Z1
    L = Z1

    length = len(a)
    n = next_power_of_2(length)

    if not is_power_of_2(length):
        G_vec = G_vec + [Z1]*(n-length)
        a = a + [0]*(n-length)

    for i in range(n):
        if i % 2 == 0:
            R = add_points(R, multiply(G_vec[i], a[i+1]))
        else:
            L = add_points(L, multiply(G_vec[i], a[i-1]))

    return L, R

# Compute the secondary diagonal L,R for a scalar vectors
def compute_secondary_diagonal_scalar(b, a):
    R = 0
    L = 0

    length = len(a)
    n = next_power_of_2(length)
    
    if not is_power_of_2(length):
        b = b + [0]*(n-length)
        a = a + [0]*(n-length)

    for i in range(len(a)):
        if i % 2 == 0:
            R = (R + (b[i] * a[i+1] % p)) % p
        else:
            L = (L + (b[i] * a[i-1] % p)) % p

    return L, R

def get_vec_1n(n):
    return [1] * n

def get_vec_2n(n):
    vec_2n = []
    for i in range(n):
        vec_2n.append(pow(2, i, p))
    
    return vec_2n

def get_vec_yn(y, n):
    vec_yn = []
    for i in range(n):
        vec_yn.append(pow(y, i, p))
    return vec_yn

def get_vec_yn_inv(y, n):
    vec_yn_inv = []
    for i in range(n):
        vec_yn_inv.append(pow(y, -i, p))
    return vec_yn_inv

# If n is a power of 2, returns n. Else, returns the next power of 2 greater than n.
def next_power_of_2(n):
    if n > 0 and (n & (n - 1)) == 0:
        return n
        
    power = 1
    while power <= n:
        power *= 2
    return power

def is_power_of_2(n):
    return n > 0 and (n & (n - 1)) == 0
def log2(n):
    return int(math.log2(n))