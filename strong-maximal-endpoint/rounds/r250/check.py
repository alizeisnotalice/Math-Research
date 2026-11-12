from fractions import Fraction as F

def envelope(y, a, v, w):
    left, right = max(a-1, y-v-w), min(a, y-v)
    if left > right:
        return F(0)
    return 1-max(left, -right, F(0))

cases = 0
for k in (1, 2, 3, 7):
    for j in range(k):
        v, w = F(j,k), F(1,k)
        for h in range(21):
            a = F(h,20)
            edges = sorted({a-1+v, v, v+w, a+v+w})
            area = sum((r-l)*envelope((r+l)/2, a, v, w)
                       for l,r in zip(edges, edges[1:]))
            assert area == w + F(1,2) + a*(1-a)
            cases += 1
print(f'{cases} exact envelope integrals passed')
