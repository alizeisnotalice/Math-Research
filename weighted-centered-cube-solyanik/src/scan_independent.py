#!/usr/bin/env python3
"""Independent double-precision scan for endpoint and critical-side arithmetic."""
from math import exp, log

MS = (2, 3, 4, 8)
NS = (1, 2, 3)
EPS = (1e-9, 1e-6, 1e-3)
ALPHAS = (0.5, 0.75, 0.9, 0.99, 0.999)


def sm(m: int, lam: float) -> float:
    return m * (1.0 - exp(log(lam) / m))


def main() -> None:
    print("FLOAT SCAN: critical threshold")
    print("m,side,epsilon,lambda,s,Phi,expected_sign,observed")
    for m in MS:
        star = ((m - 1.0) / m) ** m
        for eps in EPS:
            for side, lam in (("below", star * (1-eps)), ("above", star + (1-star)*eps)):
                s = sm(m, lam)
                phi = 1 - s
                expected = "Phi<0" if side == "below" else "Phi>0"
                observed = "PASS" if ((phi < 0) == (side == "below")) else "FAIL"
                print(f"{m},{side},{eps:.0e},{lam:.16g},{s:.16g},{phi:.16g},{expected},{observed}")

    print("\nFLOAT SCAN: star construction approach to target")
    print("n,alpha,eta,ratio,target,error,strict_average")
    for n in NS:
        for alpha in ALPHAS:
            target = 1 + 2**n * (1-alpha)/alpha
            for eta in (1e-2, 1e-4, 1e-6):
                q = (1-eta)*(1-alpha)/alpha
                ratio = 1 + 2**n*q
                strict_average = 1/(1+q)
                print(f"{n},{alpha:.12g},{eta:.0e},{ratio:.16g},{target:.16g},{target-ratio:.3e},{strict_average:.16g}")


if __name__ == "__main__":
    main()
