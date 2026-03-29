import mpmath as mp

mp.mp.dps = 80  # high precision so differences are trustworthy

def forward_differences(seq):
    """Return list diffs[k] = k-th forward differences of seq."""
    diffs = [list(seq)]
    while len(diffs[-1]) > 1:
        prev = diffs[-1]
        diffs.append([prev[j+1] - prev[j] for j in range(len(prev) - 1)])
    return diffs

def first_violation_total_monotone(seq):
    """
    Totally monotone (aka completely monotone sequence convention):
      (-1)^k * Δ^k a_n >= 0 for all k,n.
    Return (k, n, value) for the first violation, or None if none found.
    """
    diffs = forward_differences(seq)
    for k, dk in enumerate(diffs):
        sgn = (-1)**k
        for n, val in enumerate(dk):
            if sgn * val < 0:
                return (k, n, val)
    return None

# --- Define a completely monotone function f(x) as a positive mixture of exponentials ---
# f(x) = sum_j w_j * s_j^x, with 0<s_j<1, w_j>0.
ws = [mp.mpf("0.1"), mp.mpf("1.0")]
ss = [mp.mpf("0.9"), mp.mpf("0.1")]

logs = [-mp.log(s) for s in ss]  # positive numbers

def f(x):
    return mp.fsum([w * (s**x) for w, s in zip(ws, ss)])

def minus_fprime(x):
    # -f'(x) = sum w * (-log s) * s^x
    return mp.fsum([w * (s**x) * lg for w, s, lg in zip(ws, ss, logs)])

def rhs(i):
    return f(i) - f(i+1)

def solve_xi(i, max_hi=mp.mpf("1e6")):
    """
    Solve -f'(x_i) = f(i)-f(i+1) using bisection.
    For these mixtures, -f'(x) is positive and strictly decreasing, so the solution is unique.
    """
    r = rhs(i)
    lo = mp.mpf("0")
    if minus_fprime(lo) < r:
        return None  # should not happen in normal cases

    hi = mp.mpf("1")
    while minus_fprime(hi) > r:
        hi *= 2
        if hi > max_hi:
            return None

    for _ in range(250):  # plenty with mp precision
        mid = (lo + hi) / 2
        if minus_fprime(mid) > r:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

# --- Generate x_i and a_i = f(x_i) ---
N = 20
xs = [solve_xi(i) for i in range(1, N+1)]
a  = [f(x) for x in xs]

print("First few x_i:")
for i in range(5):
    print(i+1, xs[i])

print("\nFirst few a_i = f(x_i):")
for i in range(5):
    print(i+1, a[i])

# --- Check total monotonicity ---
viol = first_violation_total_monotone(a)
print("\nFirst violation (k, n, Δ^k a_n):", viol)

# Optionally print Δ^k a_1 for k=0..10
diffs = forward_differences(a)
print("\nΔ^k a_1 for k=0..10:")
for k in range(0, min(11, len(diffs))):
    print(k, diffs[k][0])
