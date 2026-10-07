"""Three exact guards for the positive cell oracle; soft test kernel is not phi."""
from fractions import Fraction as F
from pathlib import Path
import json

def product(seq):
    out=F(1)
    for v in seq: out*=v
    return out

def response(s,L,source,weights):
    return sum((w*product((1-s)*int(abs(v)<=L/2)+s*F(1,2)*int(abs(v)<=L)
                          for v in z) for z,w in zip(source,weights)),F(0))/L**len(source[0])

rounds=[]
for n in (8,32,128):
    # Correlated atoms, all coordinates retained; no coordinate-product input.
    source=[[F((i*3+j*5)%11-5,8) for i in range(n)] for j in range(7)]
    source[0]=[F(1,2)]*n  # Exact closed face at L=1.
    weights=[F(v,28) for v in range(1,8)]
    sgrid=list(map(F,[0]))+[F(1,8),F(1,2),F(7,8),F(1)]
    lgrid=[F(1),F(9,8),F(5,4),F(3,2),F(2)]
    records=[]
    for p,q in zip(sgrid,sgrid[1:]):
        star=q/(1+q-p)
        assert p<=star<=q
        for l,u in zip(lgrid,lgrid[1:]):
            point=response(star,u,source,weights)
            factor=((1+q-p)*u/l)**n
            upper=factor*point
            worst=F(0)
            for s in [p,(p+q)/2,q]:
                for L in [l,(l+u)/2,u]:
                    value=response(s,L,source,weights)
                    assert value<=upper
                    worst=max(worst,value)
            records.append(dict(p=str(p),q=str(q),l=str(l),u=str(u),star=str(star),
                                upper=str(upper),sample_max=str(worst),factor=str(factor)))
    # Degenerate softness intervals and boundary dimensions.
    for p in (F(0),F(1,2),F(1)):
        assert p/(1+p-p)==p
    closed=response(F(0),F(1),[source[0]],[F(1)])
    assert closed==1
    assert response(F(0),F(999,1000),[source[0]],[F(1)])==0
    rounds.append(dict(n=n,parameter_cells=len(records),sample_inequalities=9*len(records),
                       closed_face_jump='PASS',records=records))

# Endpoint-only softness scan misses an interior maximum, even for a positive kernel.
src=[[F(3,4),F(0),F(0),F(0)]]
end=max(response(F(0),F(1),src,[F(1)]),response(F(1),F(1),src,[F(1)]))
mid=response(F(1,2),F(1),src,[F(1)])
assert mid>end
out=dict(status='PASS',seed=None,rounds=rounds,
         endpoint_scan_guard=dict(endpoint=str(end),interior=str(mid)),
         scope='Exact positive even decreasing surrogate kernel; not the original phi, not FIRST, no weak-type budget.')
dest=Path(__file__).with_suffix('.json');dest.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(status='PASS',cells=48,inequalities=432,endpoint=str(end),interior=str(mid))))
