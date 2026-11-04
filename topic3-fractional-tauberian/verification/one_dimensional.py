#!/usr/bin/env python3
"""Reproducible checks for the one-dimensional fractional halo formulas.

This script is diagnostic evidence only.  The statements in the manuscript are
proved analytically; the numerical calculations below test the formulae, the
transition root, and the asymptotic coefficients at high precision.
"""

from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp
from scipy.optimize import brentq


mp.mp.dps = 80


def inner(t: mp.mpf, beta: mp.mpf) -> mp.mpf:
    return 1 - (1 - t) ** (1 / beta)


def outer(t: mp.mpf, beta: mp.mpf) -> mp.mpf:
    return (1 - t) ** (-1 / (1 - beta)) - 1


def transition_mp(t: mp.mpf) -> mp.mpf:
    """The unique beta where inner(t,beta)=outer(t,beta), for 0<t<1/2."""
    fn = lambda beta: inner(t, beta) - outer(t, beta)
    return mp.findroot(fn, (mp.mpf("0.4"), mp.mpf("0.5")))


def transition_scipy(t: float) -> float:
    return brentq(lambda beta: float(inner(mp.mpf(t), mp.mpf(beta)) - outer(mp.mpf(t), mp.mpf(beta))), 1e-8, 1 - 1e-8)


@dataclass(frozen=True)
class Row:
    t: str
    beta_star: str
    beta_series_error: str
    scipy_gap: str


def beta_series(t: mp.mpf) -> mp.mpf:
    # beta_*(t)=1/2-t/2-t^2/4-t^3/3-3t^4/8+O(t^5).
    return mp.mpf("0.5") - t / 2 - t**2 / 4 - t**3 / 3 - 3 * t**4 / 8


def main() -> None:
    rows: list[Row] = []
    for raw_t in ("1e-2", "1e-3", "1e-4", "1e-5"):
        t = mp.mpf(raw_t)
        root = transition_mp(t)
        scipy_root = mp.mpf(transition_scipy(float(t)))
        rows.append(
            Row(
                raw_t,
                mp.nstr(root, 24),
                mp.nstr(abs(root - beta_series(t)), 8),
                mp.nstr(abs(root - scipy_root), 8),
            )
        )
    print("transition-root checks")
    print("t,beta_star,abs(error to fourth-order t series),abs(mp-scipy)")
    for row in rows:
        print(f"{row.t},{row.beta_star},{row.beta_series_error},{row.scipy_gap}")

    # At beta=1/2 the outer branch is exactly dominant for every t in (0,1),
    # and has expansion 2t+3t^2+4t^3+5t^4+... .
    print("critical beta=1/2 checks")
    for raw_t in ("1e-2", "1e-3", "1e-4"):
        t = mp.mpf(raw_t)
        exact = outer(t, mp.mpf("0.5"))
        fourth = 2 * t + 3 * t**2 + 4 * t**3 + 5 * t**4
        print(raw_t, mp.nstr(exact, 20), mp.nstr(abs(exact - fourth), 8))


if __name__ == "__main__":
    main()
