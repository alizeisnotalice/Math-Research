"""Thin-hole band model for w(x,y)=|x|^{-gamma} in R^2.

For disjoint tangential bands B_i and radii r_i, H_i is the scaled normal
half-width of the hole |x|<t H_i.  The t->0 outer lower coefficient is J.
The module separates geometry, floating optimization, and rational checks.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Sequence

import numpy as np
from scipy.optimize import minimize


@dataclass(frozen=True)
class BandGeometry:
    bands: tuple[tuple[float, float], ...]
    radii: tuple[float, ...]
    gamma: float

    def __post_init__(self):
        if len(self.bands) != len(self.radii) or not self.bands:
            raise ValueError("bands/radii mismatch")
        if not (0.0 < self.gamma < 1.0) or any(r <= 0 for r in self.radii):
            raise ValueError("invalid gamma or radius")
        ordered = sorted(self.bands)
        if tuple(ordered) != self.bands or any(a >= b for a, b in self.bands):
            raise ValueError("bands must be ordered nonempty intervals")
        if any(self.bands[k][1] > self.bands[k + 1][0] for k in range(len(self.bands)-1)):
            raise ValueError("bands must be disjoint")

    @property
    def beta(self) -> float:
        return 1.0 - self.gamma

    def constraint_matrix(self) -> np.ndarray:
        """All exact extrema of y -> int_[y-r_i,y+r_i] H(s) ds on B_i."""
        rows = []
        ends = [z for band in self.bands for z in band]
        for i, ((left, right), radius) in enumerate(zip(self.bands, self.radii)):
            candidates = {left, right}
            for endpoint in ends:
                for y in (endpoint - radius, endpoint + radius):
                    if left <= y <= right:
                        candidates.add(y)
            for y in candidates:
                rows.append([max(0.0, min(y + radius, b) - max(y - radius, a))
                             for a, b in self.bands])
        return np.asarray(rows, dtype=float)

    def bounds(self) -> np.ndarray:
        # A row arises from a witness of some i; its allowed integral is 2 r_i^2.
        out = []
        ends = [z for band in self.bands for z in band]
        for (left, right), radius in zip(self.bands, self.radii):
            candidates = {left, right}
            for endpoint in ends:
                for y in (endpoint - radius, endpoint + radius):
                    if left <= y <= right:
                        candidates.add(y)
            out.extend([2.0 * radius * radius] * len(candidates))
        return np.asarray(out, dtype=float)

    def denominator(self) -> float:
        points = sorted({z + shift for (a, b), r in zip(self.bands, self.radii)
                         for z in (a, b) for shift in (-r, r)})
        ans = 0.0
        for left, right in zip(points[:-1], points[1:]):
            mid = (left + right) / 2.0
            height = max((r for (a, b), r in zip(self.bands, self.radii)
                          if a-r <= mid <= b+r), default=0.0)
            ans += (right-left) * height ** self.beta
        return ans

    def objective(self, h: np.ndarray) -> float:
        lengths = np.asarray([b-a for a, b in self.bands])
        return float(np.dot(lengths, np.maximum(h, 0.0) ** self.beta) / self.denominator())

    def optimize_heights(self, margin: float = 0.0) -> dict:
        matrix, bound = self.constraint_matrix(), self.bounds() * (1.0-margin)
        lengths = np.asarray([b-a for a, b in self.bands])
        # maximize concave objective subject to linear constraints: equivalent convex minimization.
        def fun(h): return -float(np.dot(lengths, np.maximum(h, 0.0) ** self.beta))
        def jac(h): return -lengths * self.beta * np.maximum(h, 1e-15) ** (self.beta-1.0)
        start = np.full(len(self.bands), min(self.radii) / 4.0)
        cons = {'type': 'ineq', 'fun': lambda h: bound - matrix @ h,
                'jac': lambda h: -matrix}
        sol = minimize(fun, start, jac=jac, bounds=[(1e-12, None)]*len(start),
                       constraints=cons, method='SLSQP', options={'ftol':1e-12,'maxiter':1000})
        return {'success': bool(sol.success), 'message': str(sol.message), 'heights': sol.x.tolist(),
                'J': self.objective(sol.x), 'max_constraint_ratio': float(np.max((matrix@sol.x)/bound))}


def fraction_constraint_check(bands: Sequence[tuple[Fraction, Fraction]], radii: Sequence[Fraction],
                              heights: Sequence[Fraction], margin: Fraction) -> dict:
    """Exact verifier for all finite extrema of the window constraints."""
    endpoints = [z for interval in bands for z in interval]
    worst = Fraction(0)
    for (left, right), radius in zip(bands, radii):
        ys = {left, right}
        for endpoint in endpoints:
            for y in (endpoint-radius, endpoint+radius):
                if left <= y <= right:
                    ys.add(y)
        for y in ys:
            integral = sum(max(Fraction(0), min(y+radius,b)-max(y-radius,a))*h
                           for (a,b),h in zip(bands,heights))
            if integral > 2*(1-margin)*radius*radius:
                raise AssertionError(f'constraint fails at y={y}: {integral}')
            worst = max(worst, integral/(2*radius*radius))
    return {'status':'PASS','worst_missing_coefficient':str(worst),
            'strict_upper':str(1-margin)}
