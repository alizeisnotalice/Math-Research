"""Exact full96 formal asymptotic Hermite calculation, no author executable."""
from collections import Counter
from fractions import Fraction as Q
from math import prod
import hashlib,json
from recheck import geometry,ROOT,OUT,read,dot

def add(x,y,a=Q(1)):
    out=dict(x)
    for k,v in y.items():out[k]=out.get(k,Q(0))+a*v
    return {k:v for k,v in out.items() if v}
def scale(x,a):return {k:v*a for k,v in x.items() if v*a}
def expansion(bs):
    t=tuple(sorted(bs))
    if t[0]==t[-1]:
        b=t[0];n=len(t)-1
        if n==0:return {'log_eps':Q(1)} if not b else {str(abs(b)):Q(1)}
        if n==1:return {} if not b else {'constant':1/b}
        assert b
        return {'constant':-1/(2*b*b)}
    return scale(add(expansion(t[1:]),expansion(t[:-1]),-1),1/(t[-1]-t[0]))

cones,_=geometry();out=[];types=set()
for case in read(ROOT/'audit/032_cases.json'):
    x=tuple(map(Q,case['x']));tau=Q(case['tau'])
    if not any(x):
        out.append({**case,'classification':'infinity_origin_power'});continue
    total={};count=Counter()
    for rays,_,_ in cones:
        rates=[]
        for r in rays:
            a,b,c=sorted(map(abs,r),reverse=True)
            rates.append(max(4*a,2*(a+b+c))+tau*(4*a+2*b))
        bs=tuple(dot(r,x)/k for r,k in zip(rays,rates));zeros=bs.count(0)
        label=(zeros,len(set(bs)));types.add(label);count[str(label)]+=1
        total=add(total,scale(expansion(bs),Q(1)/prod(rates)))
    gamma=total.get('log_eps',Q(0));assert gamma<=0
    assert gamma+sum((v for k,v in total.items() if k not in ('constant','log_eps')),Q(0))==0
    if x==(Q(1),Q(0),Q(0)):assert gamma==0
    if x==(Q(1),Q(1),Q(0)):assert gamma==(Q(0) if tau==0 else -Q(79,612))
    if x==(Q(1),Q(1),Q(1)):assert gamma==-Q(6,17)
    old=next(r for r in read(ROOT/'audit/032_results.json')['fixed_cases'] if r['x']==case['x'] and r['tau']==case['tau'])
    assert gamma==Q(old['global_Gamma_without_8pi3'])
    assert total.get('constant',Q(0))==Q(old['constant_without_8pi3'])
    assert {k:str(v) for k,v in total.items() if k not in ('constant','log_eps')}==old['log_coefficients_without_8pi3']
    out.append({**case,'node_types_full96':dict(count),'gamma':str(gamma),'formal_terms':{k:str(v) for k,v in total.items()},'classification':'infinity_log' if gamma<0 else 'finite'})
assert types=={(0,1),(0,2),(0,3),(1,2),(1,3),(2,2)}
result={'status':'PASS all6 zero/repeated types and full96 signed asymptotic sums','cases':out,'new_Green_nodes':0,'checker_sha256':hashlib.sha256(__import__('pathlib').Path(__file__).read_bytes()).hexdigest()}
(OUT/'degenerate_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['status'])
