"""Compare the two directed scans against all start-endpoint constraints."""
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
from itertools import product

def scan(widths, gaps, a):
    out = [a * widths[0]]
    for gap, width in zip(gaps, widths[1:]):
        out.append(a * width + max(0, out[-1] - gap))
    return out

def test(args):
    widths, a = args
    count = 0
    for gaps in product([F(1,4), F(1), F(3,2), F(3)], repeat=len(widths)-1):
        for ws, gs in ((widths, gaps), (widths[::-1], gaps[::-1])):
            direct = [max(a*sum(ws[h:i+1])-sum(gs[h:i]) for h in range(i+1)) for i in range(len(ws))]
            assert scan(ws, gs, a) == direct
        count += 1
    return count

if __name__ == '__main__':
    jobs = [(w,a) for w in product([1,2],repeat=4) for a in [F(1,3),F(1),F(3,2)]]
    with ProcessPoolExecutor(max_workers=12) as pool:
        print('exact catalogues',sum(pool.map(test,jobs)))
