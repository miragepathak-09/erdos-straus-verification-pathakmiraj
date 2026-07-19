"""
FINAL-ROUND CERTIFICATES for the merged atlas es121_atlas.tex.

This script certifies the material added in the final revision:

  A. The Dilation Lemma: solutions (x, y, z) of 4/n^2 with n | x,
     n^2 | y, n^2 | z biject with solutions (x/n, y/n, z/n) of 4/n
     having the last two denominators divisible by n. Checked in both
     directions on every scraped row plus the p = 1009 example.

  B. The Mordell-class witness p = 1009 = 840 + 169: a genuine Type II
     solution of 4/1009 and the corresponding Type III solution of
     4/1009^2, produced by the Master Theorem (c = 3, t = 253, w = 11).
     Conclusion: the classical obstruction concerns polynomial
     IDENTITIES, not existence of solutions.

  C. Load-bearing qualification examples at n = 9 (why the coprimality
     conditions cannot be dropped from the two classical profiles).

  D. The squares of the units mod 840 are exactly Mordell's classes
     {1, 121, 169, 289, 361, 529}.

  E. INDEPENDENT COMPUTATIONAL RE-VERIFICATION of Elsholtz--Tao's
     Proposition 1.6 [arXiv:1107.1010], a theorem essentially due to
     Schinzel and Yamamoto: for every odd PERFECT SQUARE n = m^2,
     f_I(n) = f_II(n) = 0.  This script re-verifies the vanishing,
     with no appeal to quadratic reciprocity, for 3 <= m <= 99:
       (i)  NO solution of 4/n = 1/x + 1/y + 1/z has two denominators
            divisible by n and the third coprime to m;
       (ii) NO solution has exactly one denominator divisible by n and
            the other two coprime to m.
     The search regions are the provably finite ones proved in the
     manuscript (recorded under each function below); every candidate
     triple passing all integrality/divisibility tests is rechecked by
     exact rational arithmetic. Positive controls (Type III solutions
     DO exist, e.g. 4/9 = 1/3 + 1/18 + 1/18, and the partial-divisor
     solution 4/9 = 1/4 + 1/6 + 1/36) confirm the search machinery is
     not vacuous. A naive independent brute force cross-validates the
     counts for m <= 15.

  F. Undilation of the scaffolds c in {3, 7}: Type II solutions of
     4/p for ALL 37 primes p = 1 (mod 4), p < 400.

  G. The flagship families satisfy z = lcm(x, y) = xy/gcd(x, y), with
     anchors d = gcd(x, y) equal to p (Corollaries I, III) and 2p
     (Corollary II) -- the (du, dv, duv) structural form of Lopez
     [arXiv:2206.10319].  The coprimality gcd(u, v) = 1 is certified
     SYMBOLICALLY via the five linear combinations of Remark
     rem:taxonomy:  (3s+2) - 3(s+1) = -1,  4(3s+2) - (12s+5) = 3,
     8(3s+2) - (24s+13) = 3,  (7q+6) - 7(q+1) = -1,
     4(7q+6) - (28q+17) = 7,  with 3 !| 3s+2 and 7 !| 7q+6.

  H. LOPEZ TYPE A/B PLACEMENT (the comparison resolved).  After
     undilation (x/p, y/p, z/p), every solution produced by the
     Master Theorem and the flagship identities has Lopez's Type A
     form (du, dv, duv) -- symbolically on all three families
     (d = 1 resp. 2), and numerically on every scraped/rescued row --
     with the single exception p = 193, whose undilated solution is
     the Type B solution (50, 1930, 4825) = 5*(10, 2*193, 5*193)
     displayed in [arXiv:2404.01508] (in whose census of the first
     9000 moduli 193 is the first prime with no Type A solution).
     Consistently with his Theorem 5, the scaled triple at 193^2 is
     of NEITHER form; and at p = 73 his t_L = 2, w_L = 21 parameters
     give a DIFFERENT Type A solution (21, 146, 3066) than our c = 7
     row (20, 219, 4380).

All arithmetic is exact (Fraction); nothing is done in floating point.
"""

from fractions import Fraction
from math import gcd, lcm
import sympy as sp

PASS = True


def check(cond, msg):
    global PASS
    if not cond:
        PASS = False
        raise AssertionError(msg)


def es(x, y, z, n):
    return Fraction(1, x) + Fraction(1, y) + Fraction(1, z) == Fraction(4, n)



# A. DILATION BIJECTION (both directions)

print("=" * 72)
print("A. Dilation Lemma, both directions")

# Type III solutions of 4/p^2  (x, y, z) with x = pX, y = p^2 Y, z = p^2 Z
triples = [
    (13,   52,     338,        676),
    (17,   85,     578,        2890),
    (29,   232,    2523,       20184),
    (29,   232,    3364,       6728),
    (73,   1460,   15987,      319740),
    (1009, 255277, 89591128,   2060595944),
]

for (p, x, y, z) in triples:
    n = p * p
    # downward: Type III of 4/p^2  ->  solution of 4/p with two
    # denominators divisible by p
    check(es(x, y, z, n), f"Type III identity fails: {(p, x, y, z)}")
    check(x % p == 0 and y % n == 0 and z % n == 0, "divisibility shape")
    check(gcd(x, n) == p, f"gcd profile != p for {p}")
    X, Yp, Zp = x // p, y // p, z // p
    check(es(X, Yp, Zp, p), f"dilated identity fails for p = {p}")
    check(Yp % p == 0 and Zp % p == 0, "dilated: two denoms / p")
    check(gcd(X, p) == 1, "dilated: third denom coprime (genuine Type II)")
    # upward: multiply denominators of the Type II solution by p
    x2, y2, z2 = p * X, p * Yp, p * Zp
    check((x2, y2, z2) == (x, y, z), "round trip")
    check(es(x2, y2, z2, n), "upward dilation identity")
    if p == 13:
        check((X, Yp, Zp) == (4, 26, 52), "p=13 dilated row")
        print(f"  p = 13:   4/13 = 1/4 + 1/26 + 1/52   (dilated: x,y,z = "
              f"52, 338, 676)")
print("  all 6 rows:  4/p^2 Type III  <-->  4/p Type II  verified exactly")



# B. p = 1009: Type II solutions live ON the Mordell classes

print("=" * 72)
print("B. p = 1009 = 840 + 169 (a Mordell-class prime)")

p = 1009
check(sp.isprime(p), "1009 prime")
check(p % 840 == 169, "1009 in Mordell class 169 mod 840")
check(p % 4 == 1, "1009 = 1 mod 4")
# Master Theorem data: scaffold c = 3, t = (p+3)/4, w = 11
c, t, w = 3, (p + 3) // 4, 11
check(4 * t == p + c, "p + c = 4t")
check(w * w <= t * t and (t * t) % w == 0, "w | t^2")
check((t + w) % c == 0, "c | t + w")
check(gcd(c, w) == 1, "gcd(c, w) = 1")
check(w % 3 == (-t) % 3, "w == -t (mod 3)")
Y = (t + w) // c
Z = t * (t + w) // (c * w)
check((Y, Z) == (88, 2024), "(Y, Z) = (88, 2024)")
x, y, z = p * t, p * p * Y, p * p * Z
check((x, y, z) == (255277, 89591128, 2060595944), "displayed values")
check(es(x, y, z, p * p), "4/1009^2 identity")
check(gcd(x, p * p) == p, "gcd(x, p^2) = p  (genuine Type III)")
# dilated Type II solution at the Mordell-class prime itself
X, Yp, Zp = t, p * Y, p * Z
check((X, Yp, Zp) == (253, 88792, 2042216), "dilated values")
check(es(X, Yp, Zp, p), "4/1009 identity")
check(Yp % p == 0 and Zp % p == 0 and gcd(X, p) == 1, "Type II profile")
print(f"  4/{p*p} = 1/{x} + 1/{y} + 1/{z}   (Type III, gcd(x,p^2)=p)")
print(f"  4/{p} = 1/{X} + 1/{Yp} + 1/{Zp}   (genuine Type II)")
print("  => the classical obstruction is about IDENTITIES, not solutions")



# C. Why the coprimality qualifications are load-bearing (n = 9)

print("=" * 72)
print("C. Qualification examples at n = 9")

check(es(3, 18, 18, 9), "4/9 = 1/3 + 1/18 + 1/18")
check(18 % 9 == 0 and gcd(3, 9) == 3,
      "two denoms divisible by 9, third a PARTIAL divisor")
check(es(4, 6, 36, 9), "4/9 = 1/4 + 1/6 + 1/36")
check(36 % 9 == 0 and 6 % 9 != 0 and gcd(6, 3) == 3 and gcd(4, 3) == 1,
      "one denom divisible by 9, middle a partial divisor")
print("  4/9 = 1/3 + 1/18 + 1/18   (two divisible by 9, gcd(3,9)=3 != 1)")
print("  4/9 = 1/4 + 1/6 + 1/36    (one divisible by 9, gcd(6,3)=3 != 1)")
print("  => coprimality clauses of profiles (i),(ii) cannot be dropped")



# D. Unit squares mod 840 = Mordell's list

print("=" * 72)
print("D. Unit squares mod 840")
usq = {a * a % 840 for a in range(1, 840) if gcd(a, 840) == 1}
check(usq == {1, 121, 169, 289, 361, 529}, f"unit squares {sorted(usq)}")
print(f"  {{a^2 mod 840 : gcd(a,840)=1}} = {sorted(usq)}  = Mordell's list")



# E. SQUARE-OBSTRUCTION SWEEP for odd m, 3 <= m <= 99

print("=" * 72)
print("E. Independent computational re-verification of Elsholtz--Tao")
print("   Proposition 1.6 (Schinzel--Yamamoto vanishing): odd squares")
print("   n = m^2 have f_I(n) = f_II(n) = 0;  odd 3 <= m <= 99")


def profile_i(n, m, coprime=True):
    """Solutions of 4/n with n | y, n | z and gcd(x, m) = 1 (if coprime)
    resp. gcd(x, m) > 1 (control). Region: n/4 < x <= n/2 (proved in the
    manuscript: y, z >= n forces 1/x >= 2/n); factor eq. (2)."""
    hits = []
    for x in range(n // 4 + 1, n // 2 + 1):
        g = gcd(x, m)
        if (g == 1) != coprime:
            continue
        mu = 4 * x - n
        if mu <= 0:
            continue
        for u in sp.divisors(x * x):
            if u > x:
                break
            if (u + x) % mu:
                continue
            v = (x * x) // u
            if (v + x) % mu:
                continue
            Y, Z = (u + x) // mu, (v + x) // mu
            if Y >= 1 and Z >= 1:
                check(es(x, n * Y, n * Z, n), "candidate must solve")
                hits.append((x, n * Y, n * Z))
    return hits


def profile_ii_a(n, m):
    """Exactly one denominator divisible by n, of the form y = nY,
    gcd(x, m) = gcd(z, m) = 1, x <= y <= z.
    Region (proved): n/4 < x <= n/2 and x/mu < Y <= 2x/mu."""
    hits = []
    for x in range(n // 4 + 1, n // 2 + 1):
        if gcd(x, m) != 1:
            continue
        mu = 4 * x - n
        if mu <= 0:
            continue
        Ylo = x // mu + 1
        Yhi = (2 * x) // mu
        for Y in range(Ylo, Yhi + 1):
            den = mu * Y - x
            if den <= 0:
                continue
            num = n * x * Y
            if num % den:
                continue
            z = num // den
            if z < n * Y or z % n == 0:
                continue
            if gcd(z, m) == 1:
                check(es(x, n * Y, z, n), "candidate must solve")
                hits.append((x, n * Y, z))
    return hits


def profile_ii_b(n, m):
    """Exactly one denominator divisible by n, of the form z = nZ with
    gcd(x, m) = gcd(y, m) = 1, n \\nmid y, x <= y <= z.
    Region (proved): n/4 < x <= 2n/3, x <= y <= min(2n/3, 2nx/mu),
    y > nx/mu.  Enumerated via D | n x^2, mu | D + n x, y = (D+nx)/mu,
    where D = 4xy - n(x+y) > 0 and Z = xy/D."""
    hits = []
    for x in range(n // 4 + 1, (2 * n) // 3 + 1):
        if gcd(x, m) != 1:
            continue
        mu = 4 * x - n
        if mu <= 0:
            continue
        ylo = max(x, (n * x) // mu + 1)
        yhi = min((2 * n) // 3, (2 * n * x) // mu)
        N = n * x * x
        for D in sp.divisors(N):
            if (D + n * x) % mu:
                continue
            y = (D + n * x) // mu
            if y < ylo or y > yhi:
                continue
            if gcd(y, m) != 1 or y % n == 0:
                continue
            if (x * y) % D:
                continue
            Z = (x * y) // D
            z = n * Z
            if z < y:
                continue
            check(D == 4 * x * y - n * (x + y), "D consistency")
            check(es(x, y, z, n), "candidate must solve")
            hits.append((x, y, z))
    return hits


t0 = __import__("time").time()
total_exc = 0
for m in range(3, 100, 2):
    n = m * m
    h1 = profile_i(n, m)
    h2 = profile_ii_a(n, m)
    h3 = profile_ii_b(n, m)
    if h1 or h2 or h3:
        total_exc += 1
        print(f"  EXCEPTION m = {m}: {h1} {h2} {h3}")
el = __import__("time").time() - t0
print(f"  odd m in [3, 99]: profile (i) hits = 0, profile (ii) hits = 0")
print(f"  exceptions: {total_exc}   (elapsed {el:.1f}s)")
print("  => f_I(m^2) = f_II(m^2) = 0 independently re-verified, m <= 99")
check(total_exc == 0, "obstruction sweep must be exception-free")

# --- positive controls: the SAME machinery must FIND known solutions --
print("  positive controls (machinery is not vacuous):")
c1 = profile_i(9, 3, coprime=False)          # partial-divisor profile
check((3, 18, 18) in c1, f"control (3,18,18) missing: {c1}")
c2 = profile_i(25, 5, coprime=False)
check((10, 25, 50) in c2, f"control (10,25,50) missing: {c2}")
print(f"    n = 9,  coprime=False:  finds {sorted(set(c1))[:3]} ... OK")
print(f"    n = 25, coprime=False:  finds (10, 25, 50) ... OK")


def naive_all_profiles(n, m):
    """Fully independent brute force for small n: iterate x in
    (n/4, 2n/3], y in [x, 2nx/mu], compute z exactly, then classify."""
    out_i, out_ii = [], []
    for x in range(n // 4 + 1, (2 * n) // 3 + 1):
        mu = 4 * x - n
        if mu <= 0:
            continue
        for y in range(x, (2 * n * x) // mu + 1):
            rem = Fraction(4, n) - Fraction(1, x) - Fraction(1, y)
            if rem <= 0:
                continue
            z = 1 / rem
            if z.denominator != 1 or z < y:
                continue
            z = int(z)
            dy, dz = (y % n == 0), (z % n == 0)
            if dy and dz and gcd(x, m) == 1:
                out_i.append((x, y, z))
            if dy != dz:
                others_ok = (gcd(x, m) == 1 and
                             gcd(z if dy else y, m) == 1)
                if others_ok:
                    out_ii.append((x, y, z))
    return out_i, out_ii


for m in range(3, 16, 2):
    n = m * m
    ni, nii = naive_all_profiles(n, m)
    si = profile_i(n, m)
    sii = profile_ii_a(n, m) + profile_ii_b(n, m)
    check(ni == [] and nii == [], f"naive found hits for m = {m}: {ni}{nii}")
    check(si == [] and sii == [], "fast search found hits")
    # control sets must also agree (partial-divisor solutions exist)
    ci = sorted(set(profile_i(n, m, coprime=False)))
    check(len(ci) >= 1 or m != 3, f"control list empty at m = {m}")
print("  naive cross-validation for m <= 15: identical (empty) hit sets,")
print("  and non-empty partial-divisor control sets -- machinery sound")



# F. Undilation: Type II solutions of 4/p for all 37 primes p=1 mod 4<400

print("=" * 72)
print("F. c in {3,7} scaffolds => Type II solutions of 4/p, 37 primes")

primes = [q for q in sp.primerange(5, 400) if q % 4 == 1]
check(len(primes) == 37, f"expect 37 primes, got {len(primes)}")
via_c3, via_c7 = [], []
for q in primes:
    done = False
    for c_ in (3, 7):
        t_ = (q + c_) // 4
        if 4 * t_ != q + c_:
            continue
        for w_ in sp.divisors(t_ * t_):
            if (t_ + w_) % c_ != 0 or gcd(c_, w_) != 1:
                continue
            Y_ = (t_ + w_) // c_
            Z_ = t_ * (t_ + w_) // (c_ * w_)
            # Type III at the square
            check(es(q * t_, q * q * Y_, q * q * Z_, q * q), "master id")
            # undilated Type II at the prime
            check(es(t_, q * Y_, q * Z_, q), "Type II id")
            check(gcd(t_, q) == 1 and (q * Y_) % q == 0 and (q * Z_) % q == 0,
                  "Type II profile")
            (via_c3 if c_ == 3 else via_c7).append(q)
            done = True
            break
        if done:
            break
    check(done, f"no scaffold found for p = {q}")
print(f"  all 37 primes covered: {len(via_c3)} via c=3 "
      f"(incl. 33 with t having a factor 2 mod 3), {len(via_c7)} via c=7")
check(sorted(via_c7) == [73, 193, 241, 313], "c=7-only residual set")
print(f"  c=7 residual primes: {sorted(via_c7)}")



# G. z = lcm(x, y) on the flagship families (the (du,dv,duv) form)

print("=" * 72)
print("G. Flagship families: z = lcm(x, y), anchors d = p and d = 2p")


def family(s, kind):
    if kind == 1:      # Corollary I: p = 12s+5
        p_ = 12 * s + 5
        return p_, p_ * (3 * s + 2), p_ ** 2 * (s + 1), \
            p_ ** 2 * (s + 1) * (3 * s + 2)
    if kind == 2:      # Corollary II: p = 24s+13
        p_ = 24 * s + 13
        return p_, p_ * (6 * s + 4), 2 * p_ ** 2 * (s + 1), \
            2 * p_ ** 2 * (s + 1) * (3 * s + 2)
    p_ = 28 * s + 17  # Corollary III: p = 28q+17
    return p_, p_ * (7 * s + 6), p_ ** 2 * (s + 1), \
        p_ ** 2 * (s + 1) * (7 * s + 6)


for kind, d_expect in ((1, 1), (2, 2), (3, 1)):
    for s in range(0, 60):
        p_, x_, y_, z_ = family(s, kind)
        check(es(x_, y_, z_, p_ * p_), f"family {kind} identity fails")
        check(z_ == lcm(x_, y_), f"family {kind}: z != lcm(x,y) at s={s}")
        d = gcd(x_, y_)
        check(d == d_expect * p_, f"family {kind}: anchor {d} != "
                                  f"{d_expect}*p at s={s}")
        u, v = x_ // d, y_ // d
        check(gcd(u, v) == 1 and d * u * v == z_, "not (du,dv,duv)")
    print(f"  family {kind}: z = lcm(x,y), d = {d_expect}p, "
          f"(du,dv,duv)-form   OK (s = 0..59)")

# symbolic companions of the gcd claims: the five linear combinations
# displayed in Remark rem:taxonomy, proving gcd(u, v) = 1 IDENTICALLY
# (not by numeric spot-checks)
print("-" * 72)
print("G'. Symbolic gcd(u, v) = 1 certificates (linear combinations)")
s = sp.symbols('s')
combos = [
    ("(3s+2) - 3(s+1) = -1", (3 * s + 2) - 3 * (s + 1), -1),
    ("4(3s+2) - (12s+5) = 3", 4 * (3 * s + 2) - (12 * s + 5), 3),
    ("8(3s+2) - (24s+13) = 3", 8 * (3 * s + 2) - (24 * s + 13), 3),
    ("(7q+6) - 7(q+1) = -1", (7 * s + 6) - 7 * (s + 1), -1),
    ("4(7q+6) - (28q+17) = 7", 4 * (7 * s + 6) - (28 * s + 17), 7),
]
for name, expr, val in combos:
    check(sp.simplify(expr - val) == 0, f"combination fails: {name}")
    print(f"  {name}   OK (identity in the free parameter)")
# the non-divisibility clauses (coefficient-wise: f(s) == f(0) mod m)
check(all(c % 3 == 0 for c in sp.Poly((3 * s + 2) - 2, s).coeffs()),
      "3s+2 == 2 mod 3")
check(all(c % 7 == 0 for c in sp.Poly((7 * s + 6) - 6, s).coeffs()),
      "7q+6 == 6 mod 7")
check((3 * 0 + 2) % 3 == 2 and (7 * 0 + 6) % 7 == 6, "residues")
print("  3 !| 3s+2  and  7 !| 7q+6   OK")
# ... and the gcds themselves, as polynomials over Q
pairs = [("gcd(3s+2, s+1)", 3 * s + 2, s + 1),
         ("gcd(3s+2, 12s+5)", 3 * s + 2, 12 * s + 5),
         ("gcd(3s+2, 24s+13)", 3 * s + 2, 24 * s + 13),
         ("gcd(7q+6, q+1)", 7 * s + 6, s + 1),
         ("gcd(7q+6, 28q+17)", 7 * s + 6, 28 * s + 17)]
for name, a_, b_ in pairs:
    check(sp.gcd(a_, b_) == 1, f"{name} != 1")
    print(f"  {name} = 1   OK (sympy)")
# composite coprimality used by the anchors: u coprime to v = p(s+1)
check(sp.gcd(3 * s + 2, (12 * s + 5) * (s + 1)) == 1, "Cor I u,v")
check(sp.gcd(3 * s + 2, (24 * s + 13) * (s + 1)) == 1, "Cor II u,v")
check(sp.gcd(7 * s + 6, (28 * s + 17) * (s + 1)) == 1, "Cor III u,v")
print("  u coprime to p*(s+1) resp. p*(q+1) on all three families  OK")
# numeric value-gcd sweeps, s,q = 0..500
for s_ in range(501):
    check(gcd(3 * s_ + 2, s_ + 1) == 1, "value gcd 1")
    check(gcd(3 * s_ + 2, 12 * s_ + 5) == 1, "value gcd 2")
    check(gcd(3 * s_ + 2, 24 * s_ + 13) == 1, "value gcd 3")
    check(gcd(7 * s_ + 6, s_ + 1) == 1, "value gcd 4")
    check(gcd(7 * s_ + 6, 28 * s_ + 17) == 1, "value gcd 5")
print("  numeric value-gcd sweep s, q = 0..500 on all five pairs  OK")



# H. LOPEZ TYPE A/B PLACEMENT OF ALL PRODUCED TRIPLES (the comparison
#    resolved): undilate each Type III solution of 4/p^2 and test the
#    Type A form (du, dv, duv) resp. Type B form (duv, du*a, dv*a),
#    a = modulus [arXiv:2404.01508, Definitions 3-4, Theorems 1,2,4,5]

print("=" * 72)
print("H. Lopez Type A/B placement (undilated prime-level triples)")


def undilate(p, x, y, z):
    check(x % p == 0 and y % p == 0 and z % p == 0, "undilation")
    return x // p, y // p, z // p


def typeA_params(X, Yp, Zp):
    """Return (d, u, v) with (X, Yp, Zp) = (du, dv, duv), gcd(u,v)=1,
    or None."""
    d = gcd(X, Yp)
    u, v = X // d, Yp // d
    if gcd(u, v) == 1 and d * u * v == Zp:
        return (d, u, v)
    return None


def typeB_params(X, Yp, Zp, a):
    """Return (d, u, v) with (X, Yp, Zp) = (duv, du*a, dv*a),
    gcd(u, v) = 1, or None."""
    if Yp % a or Zp % a:
        return None
    du, dv = Yp // a, Zp // a
    d = gcd(du, dv)
    u, v = du // d, dv // d
    if gcd(u, v) == 1 and d * u * v == X:
        return (d, u, v)
    return None


# --- H.1 flagship families: undilation is VISIBLY (du,dv,duv), and
#     the identity 4/p = 1/X + 1/Yp + 1/Zp holds symbolically --------
print("H.1 flagship families, symbolic Type A form after undilation")
for name, kind in (("Cor I", 1), ("Cor II", 2), ("Cor III", 3)):
    d_expect = {1: 1, 2: 2, 3: 1}[kind]
    for s_ in range(60):
        p_, x_, y_, z_ = family(s_, kind)
        X, Yp, Zp = undilate(p_, x_, y_, z_)
        check(es(X, Yp, Zp, p_), f"{name} undilated identity, {s_}")
        prm = typeA_params(X, Yp, Zp)
        check(prm is not None, f"{name}: NOT Type A at s = {s_}")
        dA, uA, vA = prm
        check(dA == d_expect, f"{name}: anchor {dA} != {d_expect}")
        u_sym = (3 * s_ + 2) if kind in (1, 2) else (7 * s_ + 6)
        v_sym = p_ * (s_ + 1)
        check((uA, vA) == (u_sym, v_sym), "u, v mismatch")
        check(gcd(uA, vA) == 1, "gcd(u, v) != 1")
        check(Yp % p_ == 0 and Zp % p_ == 0 and gcd(X, p_) == 1,
              "Type II profile")
        # square-level anchor is exactly his scaling D = p*d
        check((x_, y_, z_) == (p_ * dA * uA, p_ * dA * vA,
                               p_ * dA * uA * vA), "D = p*d scaling")
    print(f"  {name}: undilated = d*(u, v, uv), d = {d_expect}, "
          f"u coprime v, square anchor D = {d_expect}p  OK (0..59)")
# symbolic statement of the same fact on Cor I
s = sp.symbols('s')
p_sym = 12 * s + 5
X, Yp = 3 * s + 2, p_sym * (s + 1)
check(sp.simplify(1 / X + 1 / Yp + 1 / (X * Yp) - 4 / p_sym) == 0,
      "symbolic undilated identity (Cor I)")
print("  symbolic 4/p = 1/u + 1/v + 1/(uv) on Cor I  OK")

# --- H.2 every scraped/rescued/produced row -------------------------
print("H.2 computed rows (Tables scrape/rescue, Corollary 1009)")
rows = [
    (13, 52, 338, 676, 'A', (2, 2, 13)),
    (17, 85, 578, 2890, 'A', (1, 5, 34)),
    (17, 102, 289, 1734, 'A', (1, 6, 17)),
    (29, 232, 2523, 20184, 'A', (1, 8, 87)),
    (29, 232, 3364, 6728, 'A', (4, 2, 29)),
    (29, 290, 841, 8410, 'A', (1, 10, 29)),
    (73, 1460, 15987, 319740, 'A', (1, 20, 219)),
    (241, 14942, 522729, 32409198, 'A', (1, 62, 2169)),
    (313, 25040, 1175628, 23512560, 'A', (4, 20, 939)),
    (1009, 255277, 89591128, 2060595944, 'A', (11, 23, 8072)),
]
for p_, x, y, z, kind, prm in rows:
    check(es(x, y, z, p_ * p_), f"identity at p = {p_}")
    X, Yp, Zp = undilate(p_, x, y, z)
    check(es(X, Yp, Zp, p_), f"undilated identity at p = {p_}")
    got = typeA_params(X, Yp, Zp)
    check(kind == 'A' and got == prm, f"A-placement failed: {p_}")
    dA, uA, vA = got
    check((x, y, z) == (p_ * dA * uA, p_ * dA * vA, p_ * dA * uA * vA),
          "scaled Type A at the square")
    print(f"  p = {p_}:  ({X}, {Yp}, {Zp})  "
          f"Type A with d = {dA}, (u, v) = ({uA}, {vA})")
print("  (p = 29 w=1 row included: Type A with d = 1, (u,v) = (8,87))")

# --- H.3 the Type B exception p = 193 --------------------------------
print("H.3 the Type B exception at p = 193")
p_, x, y, z = 193, 9650, 372490, 931225
check(es(x, y, z, p_ * p_), "identity at 193^2")
X, Yp, Zp = undilate(p_, x, y, z)
check((X, Yp, Zp) == (50, 1930, 4825), "undilated 193 row")
check(es(X, Yp, Zp, p_), "4/193 identity")
check(typeA_params(X, Yp, Zp) is None, "193 must FAIL the Type A test")
prmB = typeB_params(X, Yp, Zp, p_)
check(prmB == (5, 2, 5), f"Type B params: {prmB}")
check((X, Yp, Zp) == tuple(5 * w for w in (10, 2 * 193, 5 * 193)),
      "his displayed Type B triple")
print("  (50, 1930, 4825) = 5*(10, 2*193, 5*193): Type B with "
      "d = 5, (u, v) = (2, 5) -- his displayed solution")
# his Theorem 4 congruence: p == -n (mod 4dn-1) with (d, n) = (5, 5)
check(193 % (4 * 5 * 5 - 1) == (-5) % (4 * 5 * 5 - 1),
      "Type B congruence, Thm 4")
# consistently with his Theorem 5, the SCALED triple at 193^2 is of
# NEITHER form:
check(typeA_params(x, y, z) is None,
      "193^2 row must fail the square-level Type A test")
check(gcd(x, y) == 1930 and 1930 * 5 * 193 == 1862450 != z,
      "documented failure figures")
check(typeB_params(x, y, z, p_ * p_) is None,
      "193^2 row must fail the square-level Type B test (his Thm 5)")
print("  scaled triple at 193^2: neither Type A nor Type B form  "
      "(consistent with his Theorem 5)")

# --- H.4 the parallel Lopez t_L-row at p = 73 ------------------------
print("H.4 his t_L = 2, w_L = 21 row at p = 73 gives a DIFFERENT "
      "Type A triple")
kL, tL, wL = 18, 2, 21                     # p = 73 = 4*18+1
check(wL % (kL + 1 + tL) == 0 and (kL + 1 + tL) % wL == 0, "w | k+1+t")
check((kL + 1 + tL) == wL and wL % (3 + 4 * tL) == (3 + 4 * tL) - 1,
      "w == -1 mod 3+4t")
dL = (kL + 1 + tL) // wL
nL = (wL + 1) // (3 + 4 * tL)
uL = (1 + nL * 73) // (4 * dL * nL - 1)
vL = nL * 73
check((dL, nL, uL, vL) == (1, 2, 21, 146), "his Thm 1 parameters")
check(typeA_params(dL * uL, dL * vL, dL * uL * vL) == (1, 21, 146),
      "his (du,dv,duv)")
check(es(dL * uL, dL * vL, dL * uL * vL, 73), "4/73 his row")
check((dL * uL, dL * vL, dL * uL * vL) == (21, 146, 3066), "his triple")
print("  his row: (21, 146, 3066), Type A with d = 1 "
      "(vs our c = 7 row (20, 219, 4380))")

# --- H.5 scaled Type A at a = 25 (his example) -----------------------
check(es(10, 25, 50, 25) and (10, 25, 50) == (5 * 2, 5 * 5, 5 * 2 * 5),
      "4/25 scaled Type A")
check(typeA_params(2, 5, 10) == (1, 2, 5), "4/5 Type A")
print("  4/25 = 5*(2, 5, 10): scaled Type A of 4/5 (his a = 25)  OK")

print()
print("=" * 72)
print("ALL FINAL-ROUND CERTIFICATES PASSED.")
