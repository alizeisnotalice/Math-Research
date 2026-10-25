"""Read-only exact geometry/budget check; not the missing H or I certificate."""
from pathlib import Path
from fractions import Fraction as F
import json
from flint import arb,ctx
OUT=Path(__file__).resolve().parent

def af(v):v=F(v);return arb(v.numerator)/v.denominator

def run():
    r=json.loads((OUT/'cb11_integral_results.json').read_text());assert r['status'].startswith('finite diagnostic incomplete')
    boxes=[tuple(map(F,b['box'])) for b in r['runs'][0]['remaining_boxes']];removed=list(map(lambda b:tuple(map(F,b)),r['discarded_corner_boxes']));allb=boxes+removed
    assert len(boxes)==247 and len(removed)==3
    assert sum(((b[1]-b[0])*(b[3]-b[2]) for b in allb),F(0))==1
    for i,a in enumerate(allb):
        assert 0<=a[0]<a[1]<=1 and 0<=a[2]<a[3]<=1
        for b in allb[i+1:]:assert min(a[1],b[1])<=max(a[0],b[0]) or min(a[3],b[3])<=max(a[2],b[2])
    assert 2*F(1,2**28)**2<F(1,10**8)**2
    for row in r['runs']:
        assert row['space_time_items']==3952<=4096 and len(row['remaining_boxes'])==247 and not row['target_certified']
    ctx.prec=128;pi=arb.pi();A=af(55200)/pi**3;J=3*pi.sqrt()*af(F(3,2)).zeta()/8;K=1600*pi*A*A.sqrt()*J
    pw=lambda x:(x.log()*af(F(2,3))).exp()
    value=pw(af(F(9,16000)))+pw(K*af(F(1,2**26))**2);bound=value*value.sqrt();assert bound<af(F(57,100000))
    print(json.dumps({'status':'PASS exact finite cover and conditional budget only','I_certified':False,'space_boxes':247,'space_time_items_per_representation':3952,'conditional_I_ball':str(bound),'scope':'no time/space grid re-evaluation; no global H/I or g/new bound claim'},indent=2))
if __name__=='__main__':run()
