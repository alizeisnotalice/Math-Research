"""Exact non-affine 3D grid-volume and crossed-edge check for R018."""

from fractions import Fraction as F
from functools import lru_cache
from itertools import product


period = (2, 3, 2)
indices = tuple(product(*(range(p) for p in period)))
vertices = tuple(product((0, 1), repeat=3))
quadrature = tuple(product(((F(0), F(1, 6)), (F(1, 2), F(2, 3)), (F(1), F(1, 6))), repeat=3))
shear = ((0, F(2, 3), F(-1, 5)), (F(3, 7), 0, F(1, 2)), (F(-2, 7), F(1, 4), 0))
diagonal = (F(1, 7), F(-2, 11), F(1, 8))


def add_axis(k, i):
    return tuple(value + (j == i) for j, value in enumerate(k))


def det(m):
    a, b, c = m
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def run(u):
    @lru_cache(None)
    def centre(k):
        h = tuple(x % p for x, p in zip(k, period))
        z = []
        for i in range(3):
            j, l = (i+1) % 3, (i+2) % 3
            cross = ((2+i)*h[j] + (3-i)*h[l]+(5+2*i)*h[j]*h[l]) * F(1, 7)
            wobble = (h[i] + 2*h[j]*h[l] + (3+i)*h[0]*h[1]*h[2])*F(1, 5+i)
            z.append(k[i] + u**3*(sum(shear[i][j]*k[j] for j in range(3))+cross)
                     + u**6*(diagonal[i]*k[i]+wobble))
        return tuple(z)

    @lru_cache(None)
    def edge(k, j):
        return tuple(y-x for x, y in zip(centre(k), centre(add_axis(k, j))))

    volume = F(0)
    pair_local_difference = F(0)
    pair_global_difference = F(0)
    for k in indices:
        cell_volume = F(0)
        mixed_pair = F(0)
        for quad in quadrature:
            theta = tuple(x for x, weight in quad)
            w = quad[0][1]*quad[1][1]*quad[2][1]
            cols = []
            for j in range(3):
                col = [F(0)]*3
                for v in vertices:
                    if v[j]:
                        continue
                    axial_weight = F(1)
                    for l in range(3):
                        if l != j:
                            axial_weight *= theta[l] if v[l] else 1-theta[l]
                    for i, delta in enumerate(edge(tuple(x+b for x,b in zip(k,v)), j)):
                        col[i] += axial_weight*delta
                cols.append(col)
            cell_volume += w*det(cols)
            mixed_pair += w*cols[1][0]*cols[0][1]/u**6
        base_pair = edge(k, 1)[0]*edge(k, 0)[1]/u**6
        pair_local_difference = max(pair_local_difference, abs(mixed_pair-base_pair))
        pair_global_difference += mixed_pair-base_pair
        volume += cell_volume

    mat = [[(1+diagonal[i]*u**6 if i == j else shear[i][j]*u**3) for j in range(3)] for i in range(3)]
    assert volume / len(indices) == det(mat)
    assert pair_local_difference > F(1, 10)
    assert abs(pair_global_difference / len(indices)) <= 10*u**3
    print(f"u={u} cells={len(indices)} Simpson-average-Jacobian minus period determinant=0")
    print("  largest single-cell mixed (01) product discrepancy =", pair_local_difference)
    print("  averaged signed discrepancy =", pair_global_difference/len(indices))


if __name__ == "__main__":
    for epsilon in (F(1, 16), F(1, 32), F(1, 64)):
        run(epsilon)
