import collections

# MOD = 998244353

# Combination under a modulo: nCk % p
# Modulo must be a prime number.
# Time complexity: O(n + log(p))
# Intermediate numerators and denominators are cached and reused for
# the same `n`.
# TODO: Share cached denominators among different `n`s.
def combmod(n, k, p=MOD, cache=collections.defaultdict(lambda: [(1, 1)])):
    if k > n:
        return 0
    if k > n - k:
        k = n - k
    # Reuse precalculated numerators and denominators.
    numden = cache[n]
    if k < len(numden):
        num, den = numden[k]
    else:
        # Resume from len(numden)-1.
        num, den = numden[-1]
        for i in range(len(numden) - 1, k):
            num = num * (n - i) % p
            den = den * (i + 1) % p
            numden.append((num, den))
    # a^(p-1) ≡ 1 mod p (by Fermat's little theorem)
    # i.e. a^(p-2) ≡ a^(-1) mod p
    inv_dev = pow(den, p - 2, p)
    return num * inv_dev % p


def comb_table(n, p=MOD):
    C = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        C[i][0] = C[i][i] = 1
    for i in range(1, n + 1):
        for j in range(1, i):
            C[i][j] = (C[i - 1][j - 1] + C[i - 1][j]) % p
    return C
