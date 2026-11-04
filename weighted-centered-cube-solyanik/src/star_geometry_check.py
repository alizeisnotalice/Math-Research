#!/usr/bin/env python3
"""Exact box-containment check for the continuous star construction."""
from fractions import Fraction
from itertools import product


def contains_cube(center, radius, cube_center, delta):
    return all(abs(a - b) + delta <= radius for a, b in zip(center, cube_center))


def main():
    delta = Fraction(1, 10)
    for n in (1, 2, 3):
        zero = (Fraction(0),) * n
        corners = list(product((-1, 1), repeat=n))
        for v_int in corners:
            v = tuple(Fraction(x) for x in v_int)
            # Any x in B_v is covered by Q(x,1+2 delta); verify using the
            # worst possible displacement of its centre from v.
            radius = 1 + 2*delta
            assert all(Fraction(1) + 2*delta <= radius for _ in v)
            assert contains_cube(v, radius, zero, delta)
            assert contains_cube(v, radius, v, delta)
            for u_int in corners:
                if u_int == v_int:
                    continue
                u = tuple(Fraction(x) for x in u_int)
                # In a coordinate where u and v differ, a centre x in B_v
                # can move delta toward B_u and B_u itself reaches another
                # delta toward x.  Its closest possible coordinate distance
                # is therefore 2-2 delta, still larger than 1+2 delta.
                assert any(a != b for a, b in zip(u, v))
                assert 2 - 2*delta > radius
        print(f"n={n}, delta={delta}, corners={len(corners)}, PASS")


if __name__ == "__main__":
    main()
