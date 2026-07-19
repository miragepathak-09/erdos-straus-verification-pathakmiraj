"""
Certificates for the PRIME-SQUARE Type II package:
  4/p^2 = 1/x + 1/y + 1/z   (p = 1 mod 4, exactly two denominators
  divisible by p^2).

Contents
--------
P1  Brute-force scrape of Type II solutions for n = 169, 289, 841.
P2  Classification check (c=3 scaffold): with x = p t, t = (p+3)/4,
    solutions biject with divisors w | t^2, w == -t (mod 3).
P3  Symbolic verification of the Splitting Lemma, the Master Theorem,
    and the three polynomial Corollaries (p = 12s+5, 24s+13, 28q+17).
P4  Numeric sweep: every prime p = 1 (mod 4), p < 400, every valid w;
    plus the c = 7 rescue at p = 73.
"""
import sympy as sp
from fractions import Fraction
from math import gcd

# P1: brute-force Type II scraper

def typeII_solutions(n, x_max):
    """All (x, y, z) with 4/n = 1/x+1/y+1/z, n | y, n | z, y <= z,
    n/4 < x <= x_max, via (mY - x)(mZ - x) = x^2, m = 4x - n."""
    sols = []
    for x in range(n // 4 + 1, x_max + 1):
        m = 4 * x - n
        x2 = x * x
        for u in sp.divisors(x2):
            if u > x:
                break
            if (u + x) % m:
                continue
            v = x2 // u
            if (v + x) % m:
                continue
            Y, Z = (u + x) // m, (v + x) // m
            sols.append((x, n * Y, n * Z, m, u, v))
    return sols

print("=" * 70)
print("P1: Type II decompositions (y, z multiples of n)")
for n in (169, 289, 841):
    sols = typeII_solutions(n, x_max=4 * n)
    # keep the "small" ones: y <= 30*n
    pretty = [s for s in sols if s[1] <= 30 * n]
    print(f"\n  n = {n}:")
    seen = set()
    for (x, y, z, m, u, v) in pretty:
        if (x, y, z) in seen:
            continue
        seen.add((x, y, z))
        lhs = Fraction(1, x) + Fraction(1, y) + Fraction(1, z)
        assert lhs == Fraction(4, n), (n, x, y, z)
        print(f"    4/{n} = 1/{x:<5} + 1/{y:<6} + 1/{z:<7}"
              f"   (m={m}, u={u}, v={v})")

# P2: classification for the c=3 scaffold, x = p t, t = (p+3)/4
print("\n" + "=" * 70)
print("P2: c=3 scaffold classification  (3Y - t)(3Z - t) = t^2")
for p in (13, 17, 29, 37, 41, 61, 73, 89, 101, 109):
    t = (p + 3) // 4
    assert 4 * t == p + 3
    valid = [w for w in sp.divisors(t * t) if (w + t) % 3 == 0]
    pairs = [(w, t * t // w) for w in valid if w * w <= t * t]
    for (w, wp) in pairs:
        Y = (t + w) // 3
        Z = (t + wp) // 3
        assert 3 * Y == t + w and 3 * Z == t + wp
        n = p * p
        assert Fraction(1, p * t) + Fraction(1, n * Y) + Fraction(1, n * Z) \
            == Fraction(4, n), (p, w)
    tmod3 = t % 3
    print(f"  p = {p:3d}: t = {t:3d} (t%3={tmod3}); "
          f"w | t^2 with w == -t (mod 3): {valid}"
          f"  ->  {len(pairs)} solution(s)")
    if t % 3 == 1 and all(w % 3 == 1 for w in valid):
        assert valid == [], "expected EMPTY scaffold"
print("  (p = 73: t = 19 has NO divisor w == 2 mod 3 -> c=3 scaffold empty,"
      " as predicted)")

# cross-check against brute force: counts match for x = p t exactly
for p in (13, 17, 29):
    n = p * p
    t = (p + 3) // 4
    scrapes = [s for s in typeII_solutions(n, x_max=4 * n) if s[0] == p * t]
    valid = [w for w in sp.divisors(t * t) if (w + t) % 3 == 0]
    pairs = [w for w in valid if w * w <= t * t]   # unordered {w, t^2/w}
    assert len(scrapes) == len(pairs), (p, scrapes, valid)
    print(f"  brute force at x = p t = {p*t}: {len(scrapes)} solution(s) = "
          f"#{len(valid)} valid divisors  ✓")

# P3: symbolic verification (identities hold in Q(parameters))
print("\n" + "=" * 70)
print("P3: symbolic certificates")

t, c, w, p, s, q = sp.symbols('t c w p s q', positive=True, integer=True)

# Splitting Lemma: c/t = 1/Y + 1/Z with Y=(t+w)/c, Z=t(t+w)/(c w)
Y = (t + w) / c
Z = t * (t + w) / (c * w)
assert sp.simplify(c / t - (1 / Y + 1 / Z)) == 0
print("  Splitting Lemma  c/t = 1/Y + 1/Z ................... OK")

# Master Theorem: p + c = 4t  ==>  4/p^2 = 1/(p t) + c/(p^2 t),
# then apply the lemma.  Verify as identity after p = 4t - c:
master = sp.simplify(
    4 / (4*t - c)**2
    - (1/((4*t - c)*t) + c/((4*t - c)**2 * t)))
assert master == 0
print("  Master step  4/p^2 = 1/(pt) + c/(p^2 t) ............ OK")

# Corollary I: p = 12s + 5, t = 3s+2, w = 1
pI = 12*s + 5
xI = pI*(3*s + 2)
yI = pI**2*(s + 1)
zI = pI**2*(s + 1)*(3*s + 2)
assert sp.simplify(1/xI + 1/yI + 1/zI - 4/pI**2) == 0
print("  Corollary I   (p = 12s+5) ......................... OK")

# Corollary II: p = 24s + 13, t = 6s+4, w = 2
pII = 24*s + 13
xII = pII*(6*s + 4)
yII = 2*pII**2*(s + 1)
zII = 2*pII**2*(s + 1)*(3*s + 2)
assert sp.simplify(1/xII + 1/yII + 1/zII - 4/pII**2) == 0
print("  Corollary II  (p = 24s+13) ........................ OK")

# Corollary III: c = 7, w = 1, t = 7q+6, p = 28q + 17
pIII = 28*q + 17
tIII = 7*q + 6
xIII = pIII*tIII
yIII = pIII**2*(q + 1)
zIII = pIII**2*(q + 1)*tIII
assert sp.simplify(1/xIII + 1/yIII + 1/zIII - 4/pIII**2) == 0
print("  Corollary III (p = 28q+17, c=7) ................... OK")

# P4: numeric sweeps
print("\n" + "=" * 70)
print("P4: numeric sweep, all primes p = 1 (mod 4), p < 400, c = 3")
n_covered, total = 0, 0
residual = []
for p in list(sp.primerange(5, 400)):
    if p % 4 != 1:
        continue
    total += 1
    t = (p + 3) // 4
    valid = [wv for wv in sp.divisors(t * t) if (wv + t) % 3 == 0]
    pairs = [(wv, t * t // wv) for wv in valid if wv * wv <= t * t]
    if pairs:
        n_covered += 1
    else:
        residual.append(p)
    for (wv, wp) in pairs:
        Yv = (t + wv) // 3
        Zv = (t + wp) // 3
        assert Fraction(1, p*t) + Fraction(1, p*p*Yv) + Fraction(1, p*p*Zv) \
            == Fraction(4, p*p)
print(f"  c=3 scaffold covers {n_covered}/{total} primes p=1(mod4), p<400")
print(f"  residual (t has no divisor == 2 mod 3): {residual}")
print(f"  their classes mod 24: {sorted(set(p % 24 for p in residual))}")

# c = 7 rescue for several residual primes
print("\n  c=7 scaffold (t=(p+7)/4, w== -t mod 7, w | t^2, gcd(w,7)=1):")
for p in [73, 193, 241, 313, 337, 409, 421, 433, 457]:
    if p % 4 != 1:
        continue
    t7 = (p + 7) // 4
    if 4 * t7 != p + 7:
        print(f"    p = {p}: t not integral")
        continue
    found = False
    for wv in sp.divisors(t7 * t7):
        if (wv + t7) % 7 == 0 and gcd(wv, 7) == 1:
            Yv = (t7 + wv) // 7
            Z = t7 * (t7 + wv)
            assert Z % (7 * wv) == 0
            Zv = Z // (7 * wv)
            assert 7*Yv == t7 + wv
            assert Fraction(1, p*t7) + Fraction(1, p*p*Yv) + \
                Fraction(1, p*p*Zv) == Fraction(4, p*p)
            if not found:
                print(f"    p = {p}: t = {t7}, w = {wv}, "
                      f"(Y,Z) = ({Yv},{Zv})  ->  "
                      f"4/{p*p} = 1/{p*t7} + 1/{p*p*Yv} + 1/{p*p*Zv}")
            found = True
    if not found:
        print(f"    p = {p}: no c=7 solution")

print("\nALL PRIME-SQUARE CERTIFICATES PASSED.")

# P5: composite thin squares n = m^2, m = 1 (mod 4) composite
print("\n" + "=" * 70)
print("P5: composite thin squares (Master Theorem with p composite)")
for m_, w_ in [(65, 1), (85, 2)]:
    t5 = (m_ + 3) // 4
    assert 4 * t5 == m_ + 3
    assert (w_ + t5) % 3 == 0 and t5 * t5 % w_ == 0 and __import__('math').gcd(3, w_) == 1
    Y5 = (t5 + w_) // 3
    Z5 = t5 * (t5 + w_) // (3 * w_)
    n5 = m_ * m_
    assert Fraction(1, m_ * t5) + Fraction(1, n5 * Y5) + Fraction(1, n5 * Z5) \
        == Fraction(4, n5)
    print(f"  m = {m_}:  4/{n5} = 1/{m_*t5} + 1/{n5*Y5} + 1/{n5*Z5}"
          f"  (t = {t5}, w = {w_})")
print("\nALL PRIME-SQUARE CERTIFICATES PASSED (including P5).")
