"""Independent integer GF products and exact rational Walsh audits; no MC."""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_right
import math,json,hashlib

HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'sidon'
r=json.loads((SRC/'results.json').read_text());reg=json.loads((SRC/'registration.json').read_text())
assert hashlib.sha256((SRC/'probe.py').read_bytes()).hexdigest()==reg['script_sha256']
assert hashlib.sha256((SRC/'registration.json').read_bytes()).hexdigest()==r['registration_sha256']
fields=[]
for row in r['fields']:
    B=row['B'];q=row['q'];K=q-1;poly=int(row['polynomial_hex'],16)
    def mul(x,y):
        z=0
        while y:
            if y&1:z^=x
            y>>=1;x<<=1
            if x&q:x^=poly
        return z
    def power(x,n):
        y=1
        while n:
            if n&1:y=mul(y,x)
            x=mul(x,x);n>>=1
        return y
    W=row['W'];w=[x for x in W if x]
    assert len(W)==B and {a^b for a in W for b in W}==set(W)
    ps=[mul(a,b)for i,a in enumerate(w)for b in w[i:]]
    assert len(set(ps))==len(ps)==row['Sidon_product_count']
    inv={x:power(x,q-2)for x in w}
    D={mul(a,inv[b])for a in w for b in w if a!=b}
    assert len(D)==(B-1)*(B-2)==row['distinct_nontrivial_ratio_count']
    s=row['minimizing_primitive_exponent'];alpha=power(2,s)
    assert math.gcd(s,K)==1 and alpha not in D and mul(alpha,alpha)not in D
    overlap=len(D & {mul(alpha,x)for x in D})
    assert overlap==row['minimum_adjacent_defect_overlap_among_short_good_primitive']
    fields.append({'B':B,'independent_product_count':len(ps),'ratio_count':len(D),
                   'short_lag_check':True,'selected_adjacent_overlap':overlap})

walsh=[];high_checks=0
for B in [4,8]:
    cuts=[F(2*y+1,2*B)for y in range(B)]
    prefixes={}
    for a in range(1,B):
        seq=[F(0)]
        for y in range(B):seq.append(seq[-1]+F(1-2*((a&y).bit_count()%2),B))
        prefixes[a]=seq
    def f(a,x):return prefixes[a][bisect_right(cuts,x)]
    def integral(a,c,lo,hi,shift=F(0)):
        if lo==hi:return F(0)
        grid=sorted({lo,hi,*[x for x in cuts if lo<x<hi],*[x-shift for x in cuts if lo<x-shift<hi]})
        return sum(((v-u)*f(a,(u+v)/2)*f(c,(u+v)/2+shift)for u,v in zip(grid,grid[1:])),F(0))
    h=B//2;EH=F(1,12)+F(1,6*B*B);EL=F(1,48)+F(1,6*B*B)
    assert integral(h,h,F(0),F(1))==EH
    for a in range(1,B):
        qlow=a&-a
        assert integral(a,a,F(0),F(1))==F(qlow*qlow,3*B*B)+F(1,6*B*B)
    worst=F(0)
    for j in range(18):
        t=F(j,17)
        e0=integral(h,h,F(0),1-t);e1=integral(h,h,1-t,F(1))
        et0=integral(h,h,F(0),t);et1=integral(h,h,t,F(1))
        assert e0+e1==et0+et1==e0+et0==et1+e1==EH
        high_checks+=4
        first=max(abs(integral(a,c,F(0),1-t,t))for a in range(1,B)for c in range(1,B)if a!=c)
        second=max(abs(integral(a,c,1-t,F(1),t-1))for a in range(1,B)for c in range(1,B)if a!=c)
        total=first+second
        assert total<=EL or (total-EL)**2<=EH*EL
        worst=max(worst,total)
    walsh.append({'B':B,'theta_grid_denominator':17,'theta_count':18,
                  'all_two_pair_bounds_passed':True,'maximum_exact_sum':str(worst)})

bud=json.loads((SRC/'budget_results.json').read_text())
for row in bud['rows']:
    B=row['B'];K=B**3-1;D=(B-1)*(B-2)
    assert row['primitive_count_sufficient_exact']==(K**3>(6*D)**4)
up=(F(1,4)+F(34,261))/F(100,261)
lo=-(F(1,4)+F(1,116))/F(100,261)
assert up==F(397,400)and lo==F(-27,40)
out={'scope':'No random/source simulation; independent binary-polynomial GF multiplication and Fraction Walsh integration.',
     'source_hash_checks_passed':True,'fields':fields,'Walsh':walsh,
     'high_character_complement_checks':high_checks,
     'far_upper_exact':str(up),'far_lower_exact':str(lo),
     'note':'f_h reflection is a.e.; right-continuous jump values need not match pointwise.'}
(HERE/'sidon_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
