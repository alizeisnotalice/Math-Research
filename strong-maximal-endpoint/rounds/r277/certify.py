from fractions import Fraction as Q
for den in range(2,15):
    vals=[Q(i,den) for i in range(den+1)]
    weights=[Q(i+1) for i in range(den+1)]
    for s in range(1,9):
        by_levels=sum((hi**s-lo**s)*sum(w for v,w in zip(vals,weights) if v>=hi)
                      for lo,hi in zip(vals,vals[1:]))
        assert by_levels==sum(w*v**s for v,w in zip(vals,weights))
    for i in range(1,den):
        p=vals[i]
        for a in vals[:i]:
            window=sum(w*min(1,max(0,(v-a)/(p-a))) for v,w in zip(vals,weights))
            assert sum(w for v,w in zip(vals,weights) if v>p)<=window
            assert window<=sum(w for v,w in zip(vals,weights) if v>a)
print('Exact layer-cake and strict-window comparisons passed for 13 discrete laws')
