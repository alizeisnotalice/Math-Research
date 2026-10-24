"""Accurate A_alpha of Psi = (r^-alpha - R^-alpha)_+ via A_alpha r^-alpha = C_alpha r^-2alpha (Fourier identity)
and Psi = r^-alpha - min(r^-alpha, R^-alpha):  A Psi = C_alpha r^-2alpha - A g,  g = min(r^-alpha, R^-alpha)."""
import numpy as np
from math import gamma, cos, pi
from scipy.integrate import quad


def C_alpha(alpha):
    if abs(alpha - 1) < 1e-12:
        return 0.0
    return -8 * gamma(1 - alpha) * gamma(2 * alpha) * cos(pi * alpha / 2) ** 2 / gamma(alpha + 1)


def H_g(alpha, beta, R, u, v, Tcut=2000.0, tol=1e-12):
    Rb = R ** (-beta)
    g = lambda rho: Rb if rho <= R else rho ** (-beta)
    r = u + v
    g0 = g(r)
    f = lambda t: (g(abs(u + t) + v) + g(abs(u - t) + v) - 2 * g0) * t ** (-1 - alpha)
    pts = {u}
    for t in (R - v - u, u + R - v, u - (R - v)):
        if t > 0:
            pts.add(t)
    pts = sorted(p for p in pts if 1e-15 < p < Tcut)
    edges = [0.0] + pts + [Tcut]
    tot = 0.0
    kinkset = set(pts)
    # first interval [0, edges[1]]: no kink inside, second difference = g(r+t)+g(r-t)-2g(r) (t < u).
    # Near t = 0 use the Taylor expansion (avoids catastrophic cancellation).
    t0 = edges[1]
    if r < R:
        d2 = d4 = 0.0
    else:
        d2 = beta * (beta + 1) * r ** (-beta - 2)
        d4 = beta * (beta + 1) * (beta + 2) * (beta + 3) * r ** (-beta - 4)
    eps = min(t0, 0.02 * r)
    tot += d2 * eps ** (2 - alpha) / (2 - alpha) + d4 / 12 * eps ** (4 - alpha) / (4 - alpha)
    edges = [eps] + edges[1:] if eps < t0 else edges[1:]
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < 1e-15:
            continue
        sub = [a, b]
        # a kink p close to t = 0 makes the integrand peak like p^-alpha just right of p: refine geometrically
        if a in kinkset and a < 0.05:
            g_ = [a * (1 + 2.0 ** k) for k in range(-6, 80) if a * (1 + 2.0 ** k) < b]
            sub = [a] + g_ + [b]
        for a2, b2 in zip(sub[:-1], sub[1:]):
            val, err = quad(f, a2, b2, epsabs=tol, epsrel=1e-12, limit=500)
            tot += val
    # tail: g(|u+-t|+v) ~ (t + u + v)^-beta, (t - u + v)^-beta for t > Tcut
    tail = quad(lambda t: ((t + u + v) ** (-beta) + (t - u + v) ** (-beta)) * t ** (-1 - alpha), Tcut, np.inf,
                epsabs=1e-16, limit=200)[0]
    tot += tail - 2 * g0 * Tcut ** (-alpha) / alpha
    return 2 * tot


def A_psi(alpha, R, u, v):
    u, v = abs(u), abs(v)
    r = u + v
    if abs(r - R) < 1e-12:
        return np.inf          # convex kink of Psi at r = R: A Psi = +inf there
    return C_alpha(alpha) * r ** (-2 * alpha) - (H_g(alpha, alpha, R, u, v) + H_g(alpha, alpha, R, v, u))
