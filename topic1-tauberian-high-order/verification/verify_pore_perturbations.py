#!/usr/bin/env python3
"""Exact finite checks for the mass bounds in lower_round003.md.

This validates transcription on explicit measurable periodic sets, not a
global optimization and not a replacement for the analytic proof.
"""

from fractions import Fraction as F
from itertools import product


def one_dimensional_mass(x, radius, left, right):
    lo, hi = x - radius, x + radius
    k_min = (lo - right).__floor__() - 1
    k_max = (hi - left).__ceil__() + 1
    return sum(
        (max(F(0), min(hi, k + right) - max(lo, k + left))
         for k in range(k_min, k_max + 1)),
        F(0),
    )


def periodic_box_mass(x, radius, box):
    result = F(1)
    for coordinate, (left, right) in zip(x, box):
        result *= one_dimensional_mass(coordinate, radius, left, right)
    return result


def main():
    intermediate_checks = 0
    large_checks = 0
    small_density_checks = 0
    for n in (2, 3, 4, 5):
        for u in (F(1, 100), F(1, 500)):
            t = u**n
            delta = u / (1 + 2 * u)
            m0 = (2 * delta)**n
            width = delta * u / 2
            core = [(-delta, delta)] * n
            removed = [(delta - width, delta)] + [(-delta, delta)] * (n - 1)
            satellite = [(F(1, 4), F(1, 4) + width)] + [(-delta, delta)] * (n - 1)
            a = width * (2 * delta)**(n - 1)
            p = a
            q = m0 - a + p
            contact_slope = m0 / (2 * delta)
            za = 2 * delta * (3**n - 1) * a / m0
            assert delta <= F(1, 6 * n)
            assert q <= 2 * m0
            assert m0 - a > t
            assert F(1, 4) > delta

            # Reflection symmetry is not assumed after removal/addition.
            for x in product((-delta, F(0), delta), repeat=n):
                r0 = 1 - delta - max(abs(v) for v in x)
                for radius in (delta / 2, delta, 2 * delta):
                    baseline = periodic_box_mass(x, radius, core)
                    assert baseline >= radius**n
                    small_density_checks += 1

                radii = (
                    F(1, 2),
                    F(3, 4),
                    r0 - delta / 2,
                    r0,
                    r0 + za / 2,
                    r0 + za,
                    r0 + delta / 2,
                    F(1),
                )
                for radius in radii:
                    if not F(1, 2) <= radius <= 1:
                        continue
                    base_mass = periodic_box_mass(x, radius, core)
                    removed_mass = periodic_box_mass(x, radius, removed)
                    satellite_mass = periodic_box_mass(x, radius, satellite)
                    actual = base_mass - removed_mass + satellite_mass
                    assert satellite_mass >= p
                    z = max(F(0), radius - r0)
                    expected_mass = q + contact_slope * max(F(0), z - za)
                    assert actual >= expected_mass, (n, u, x, radius, actual, expected_mass)
                    assert actual / (2 * radius)**n >= q / (2 * (r0 + za))**n
                    intermediate_checks += 1

                for radius in (F(1), F(3, 2), F(9, 4), F(21, 4)):
                    base_mass = periodic_box_mass(x, radius, core)
                    removed_mass = periodic_box_mass(x, radius, removed)
                    satellite_mass = periodic_box_mass(x, radius, satellite)
                    actual_fraction = (base_mass - removed_mass + satellite_mass) / (2 * radius)**n
                    expected_fraction = (3 * delta / 2)**n - 2**n * a
                    assert actual_fraction >= expected_fraction
                    assert expected_fraction > t
                    large_checks += 1

            print(f"explicit n={n}, u={u}, a=p=m0*u/4: PASS")

    print(f"small-radius baseline density: {small_density_checks} exact samples PASS")
    print(f"intermediate-radius perturbed mass bounds: {intermediate_checks} exact samples PASS")
    print(f"all-large-radius lower-bound samples: {large_checks} exact samples PASS")
    print("The general continuum and all-radius proofs are in 历史提交 2316b612 的 paper/sections/perturbations.tex.")


if __name__ == "__main__":
    main()
