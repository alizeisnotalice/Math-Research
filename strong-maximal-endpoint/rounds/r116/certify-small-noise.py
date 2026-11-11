"""Directed integer certificate for two arbitrary actions, b in [0,1], s <= 1/3."""
from concurrent.futures import ProcessPoolExecutor

SCALE = 10**9
M = 256
J = 1024


def logistic_tables():
    # Exact floor of F(k log(101/100)) in fixed-point units.
    grid = []
    a = b = 1
    for _ in range(J + 1):
        grid.append(SCALE * (a - b) // (a + b))
        a *= 101
        b *= 100
    # The previous CDF is at most (floor+1)/SCALE.
    mass = [max(0, grid[j] - grid[j - 1] - 1) for j in range(1, J + 1)]
    return grid, mass


def certify(block):
    grid, mass = logistic_tables()
    values = []
    for i in range(*block):
        numer = 0
        for j, weight in enumerate(mass, 1):
            k = (128 * (M + i) - (i + 1) * j) // M
            if k <= 0:
                break
            numer += weight * grid[k]
        values.append((numer, i))
    return min(values)


if __name__ == '__main__':
    blocks = [(i, min(i + 16, M)) for i in range(0, M, 16)]
    with ProcessPoolExecutor(max_workers=12) as pool:
        worst, index = min(pool.map(certify, blocks))
    assert 90 * worst > 46 * SCALE**2
    print('minimum_lower_numerator', worst, 'denominator', SCALE**2, 'bin', index, 'of', M)
    print('after_smoothing', 90 * worst - SCALE**2, '/', 90 * SCALE**2)
