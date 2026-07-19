"""
Certificates for the Divisor Extraction identities:
  Part 1:  n = 240k + 121,  r = 19  (k = 19j + 1)
  Part 2:  n = 240k + 169,  r = 7   (k = 7j + 3)  and  r = 11 (k = 11j + 2)

Every identity is checked (a) symbolically in the free variable j (sympy),
and (b) numerically for j = 0..25 with exact Fraction arithmetic.
"""
import sympy as sp
from fractions import Fraction

j = sp.symbols('j')


def s_of(r):
    """s = (r+1)/4, integral iff d = r == 3 (mod 4)."""
    assert r % 4 == 3
    return (r + 1) // 4


def solve_k0(r0, r):
    """Unique k0 mod r with 240*k0 + r0 = 0 (mod r). Requires gcd(240, r) = 1."""
    a, b = 240 % r, r0 % r
    k0 = (pow(a, -1, r) * (-b)) % r
    assert (240 * k0 + r0) % r == 0
    return k0


def two_term_coeffs(r):
    """4/n = 1/(s m) + 1/(r s m), n = r m."""
    s = s_of(r)
    return (s, r * s)


def three_term_coeffs(r):
    """Split 1/(s m) = 1/((s+1) m) + 1/(s(s+1) m)."""
    s = s_of(r)
    return (s + 1, s * (s + 1), r * s)


def verify_family(r0, r, expanded_two, expanded_three):
    k0 = solve_k0(r0, r)
    n0 = 240 * k0 + r0
    M = 240 * r                    # sub-modulus
    m0 = n0 // r                   # n0 = r*m0
    s = s_of(r)
    m = 240 * j + m0
    n = M * j + n0

    # (a) symbolic checks in j, both coefficient forms
    t2 = two_term_coeffs(r)
    t3 = three_term_coeffs(r)
    assert sp.simplify(sum(1 / (c * m) for c in t2) - 4 / n) == 0, "two-term FAIL"
    assert sp.simplify(sum(1 / (c * m) for c in t3) - 4 / n) == 0, "three-term FAIL"

    # (b) the explicit expanded polynomials stated in the paper
    for poly in expanded_two + expanded_three:
        poly_j = poly(j)
        assert all(poly_j.subs(j, v) > 0 for v in range(6)), "positivity FAIL"
    assert sp.simplify(sum(1 / p(j) for p in expanded_two) - 4 / n) == 0
    assert sp.simplify(sum(1 / p(j) for p in expanded_three) - 4 / n) == 0

    # (c) exact numeric checks, j = 0..25
    for v in range(26):
        Nv, mv = M * v + n0, 240 * v + m0
        assert Nv == r * mv
        assert sum(Fraction(1, c * mv) for c in t2) == Fraction(4, Nv)
        assert sum(Fraction(1, c * mv) for c in t3) == Fraction(4, Nv)
    return dict(k0=k0, n0=n0, M=M, m0=m0, s=s)


#  Part 1: r0 = 121, r = 19 
res = verify_family(
    121, 19,
    expanded_two=(lambda j: 1200*j + 95, lambda j: 22800*j + 1805),
    expanded_three=(lambda j: 1440*j + 114, lambda j: 7200*j + 570,
                    lambda j: 22800*j + 1805),
)
print("Part 1 (121 mod 240, r = 19):", res)

# Consistency with previous work: r = 7 extraction inside n = 121 (mod 240)
res_old = verify_family(
    121, 7,
    expanded_two=(lambda j: 480*j + 446, lambda j: 3360*j + 3122),
    expanded_three=(lambda j: 720*j + 669, lambda j: 1440*j + 1338,
                    lambda j: 3360*j + 3122),
)
print("Check  (121 mod 240, r = 7 ):", res_old, " <- recovery of identity (D)")

#  Part 2: r0 = 169, r = 7 and r = 11 
res = verify_family(
    169, 7,
    expanded_two=(lambda j: 480*j + 254, lambda j: 3360*j + 1778),
    expanded_three=(lambda j: 720*j + 381, lambda j: 1440*j + 762,
                    lambda j: 3360*j + 1778),
)
print("Part 2 (169 mod 240, r = 7 ):", res)

res = verify_family(
    169, 11,
    expanded_two=(lambda j: 720*j + 177, lambda j: 7920*j + 1947),
    expanded_three=(lambda j: 960*j + 236, lambda j: 2880*j + 708,
                    lambda j: 7920*j + 1947),
)
print("Part 2 (169 mod 240, r = 11):", res)

#  Covering table: k0 mod r for small primes r = 3 (mod 4) 
primes34 = [p for p in range(7, 72) if sp.isprime(p) and p % 4 == 3]
assert primes34 == [7, 11, 19, 23, 31, 43, 47, 59, 67, 71]
print("\nCovering table (k0 mod r such that r | 240k + r0):")
for r0, name in ((121, "121"), (169, "169")):
    row = []
    for r in primes34[:8]:
        k0 = solve_k0(r0, r)
        n0 = 240 * k0 + r0
        row.append((r, k0, n0, sp.factorint(n0)))
    print(f"  r0 = {name}:")
    for r, k0, n0, fac in row:
        print(f"    r = {r:2d}:  k = {k0:2d} (mod {r:2d}),  n0 = {n0:5d} = {fac}")

print("\nALL DIVISOR-EXTRACTION CERTIFICATES PASSED.")
