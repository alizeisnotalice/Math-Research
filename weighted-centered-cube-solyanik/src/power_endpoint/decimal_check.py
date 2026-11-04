"""Finite Decimal cross-checks of complete 1D envelopes; not a set-theoretic proof."""
from decimal import Decimal as D, getcontext

getcontext().prec = 180


def power(x, a):
    return (a * x.ln()).exp()


threshold = (D(5) - D(17).sqrt()) / 2
print("tau threshold:", threshold)
for tau, feasible in [(D('.2'), True), (D('.4'), True), (D('.44'), False)]:
    left, right = tau / 2, (1 - 2 * tau) / (1 - tau)
    assert (left < right) == feasible
    print("window", tau, "nonempty:", feasible)
    if not feasible:
        continue
    exponent = (left + right) / 2
    for k in [20, 64, 128, 256, 512]:
        t = D(2) ** -k
        z = power(t, tau)
        cutoff = power(t, -exponent)
        ac = power(1 + 1 / cutoff, tau) + power(1 - 1 / cutoff, tau)
        assert cutoff > 2 and ac > 2 * z
        lower = z / (power(1 + t, tau) - z)
        long_bound = 2 * z / (ac - 2 * z)
        short_bound = (2 * t * power(2 * (1 + cutoff), 1 - tau)
                       / ((1 - t) * (1 - power(D(2), tau - 1))))
        upper = long_bound + short_bound
        assert 0 < lower <= upper
        normalized = ((lower - z) / z**2, (upper - z) / z**2)
        print("tau,k,second-scale brackets:", tau, k,
              *(f"{x:.18E}" for x in normalized))
    assert all(abs(x - 1) < D('.01') for x in normalized)
print("PASS: 3 window decisions, 10 complete-envelope comparisons, 2 tail checks.")
