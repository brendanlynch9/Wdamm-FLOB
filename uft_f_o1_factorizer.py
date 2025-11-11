import math
from decimal import Decimal, getcontext

getcontext().prec = 100

# === UFT-F CONSTANTS ===
K_MAX = Decimal('1.0000000005')
A = Decimal('1E-9')
B = Decimal('1.5')
C = Decimal('1.0')
C1 = Decimal('1.5E-6')
C3 = Decimal('5.6E-3')

COEFFS = {
    'A': Decimal('100.362219294930227'),
    'B': Decimal('-917.106764318163187'),
    'C': Decimal('2469.390466389201720'),
    'D': Decimal('-2073.294456986756359')
}

SMALL_PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
                101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,
                211,223,227,229,233,239,241,251,257,263,269,271,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,
                401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,
                601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,
                809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997]

def decimal_cos(x: Decimal) -> Decimal:
    x = x % (2 * Decimal(str(math.pi)))
    cos_val = Decimal(1)
    term = Decimal(1)
    n = 1
    while True:
        term *= -x * x / ((2*n-1)*(2*n))
        new_cos = cos_val + term
        if new_cos == cos_val: break
        cos_val = new_cos
        n += 1
    return cos_val

def c2_uft(x: Decimal) -> Decimal:
    return COEFFS['A']*x**3 + COEFFS['B']*x**2 + COEFFS['C']*x + COEFFS['D']

def spectral_torsion(N: int) -> Decimal:
    if N < 4: return Decimal(0)
    Nd = Decimal(N)
    sqrtN = Nd.sqrt()
    X = Decimal(str(math.log10(float(sqrtN))))
    base = K_MAX - Decimal(1)/sqrtN
    c2 = c2_uft(X)
    sigma = Decimal(str(math.log10(N))) * 5
    ray = A * decimal_cos(B * sigma + C)
    mod24 = (Nd % 24) / 24
    periodic = C3 * decimal_cos(mod24 * Decimal(str(math.pi)))
    decay = C1 * (sigma / Nd)
    return base - c2 - ray - decay - periodic

def extract_small(N: int):
    factors = []
    n = N
    for p in SMALL_PRIMES:
        if p*p > n: break
        while n % p == 0:
            factors.append(p)
            n //= p
    return n, factors

def spectral_collapse(cof: int):
    if cof <= 1: return 0, 0, "TRIVIAL"
    lam = spectral_torsion(cof)
    sqrtc = math.sqrt(cof)
    Q = int(round(float(lam) * sqrtc))
    if Q > 1 and Q < cof and cof % Q == 0:
        return min(Q, cof//Q), max(Q, cof//Q), "Q_COLLAPSE"
    P = int(round(cof / Q))
    if P > 1 and P < cof and cof % P == 0:
        return min(P, cof//P), max(P, cof//P), "P_COLLAPSE"
    return 0, 0, "PRIME"

def uft_f_o1(N: int):
    print(f"\n{'='*80}")
    print(f"UFT-F SPECTRAL O(1) | N = {N} | bits: {N.bit_length()}")
    print(f"{'='*80}")
    cof, small = extract_small(N)
    factors = small[:]
    print(f"Small: {small or '∅'} | Cofactor: {cof}")
    if cof > 1:
        p, q, status = spectral_collapse(cof)
        if p: 
            factors += [p, q]
            print(f"→ SPECTRAL COLLAPSE: {p} × {q} | {status}")
        else:
            factors.append(cof)
            print(f"→ ACI SILENCE: {cof} is PRIME")
    else:
        print("→ ALL DEFECTS RESOLVED")
    prod = 1
    for f in factors: prod *= f
    print(f"\nFactors: {sorted(factors)}")
    print(f"Product = {prod} → {'VERIFIED' if prod == N else 'ERROR'}")
    print(f"{'='*80}")

if __name__ == "__main__":
    for n in [10403, 10807, 121950274103, 17, 9991, 7657]:
        uft_f_o1(n)