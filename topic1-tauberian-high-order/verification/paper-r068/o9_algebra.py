"""Finite Fraction checks of the upgraded O(av^9) ledger and normal peak.
This is arithmetic corroboration, not a proof of geometric source estimates.
"""
from fractions import Fraction as F
import datetime, json, os, sys

start=datetime.datetime.now().astimezone().isoformat()
def add(a,b): return [x+y for x,y in zip(a,b)]
def scale(a,c): return [x*c for x in a]
def mul(a,b):
    return [sum((a[j]*b[i-j] for j in range(i+1)),F(0)) for i in range(3)]
one=[F(1),F(0),F(0)]
v=[F(0),F(1),F(0)]
inv=[F(1),F(-1),F(1)]
r=[F(1,3),F(20,81),F(0)]
r2=mul(r,r);r3=mul(r2,r);r4=mul(r2,r2);r5=mul(r4,r)
rows=[]
for n in (2,3,4,8,32):
    base=add(r2,scale(r3,-2))
    q4=scale(r4,F(5,3));q5=scale(r5,F(13,4))
    phi=add(base,scale(mul(mul(v,inv),q4),2))
    phi=add(phi,scale(mul(mul(v,v),q5),-2))
    phi=add(phi,scale(mul(mul(v,v),r4),-n*(n-2)))
    expected=[F(1,27),F(10,243),F(-91+324*n-162*n*n,13122)]
    assert phi==expected
    peak=scale(mul(mul(inv,inv),phi),n*n)
    k=F(n*n,27);gamma=F(10,9)*k
    c8=-F(13,18)-F(n*(n-2),3)+F(400,243)
    qn=[k,-2*k+gamma,3*k-3*gamma+k*c8]
    assert peak==qn
    # Whole slot and the actual-M denominator through relative order two.
    Kn=F(n*(n-1),2)*k
    local=[Kn,-F(26,9)*Kn,(F(40,9)+c8)*Kn]
    if n==2: local[2]+=F(8,27) # independent inside fee
    denominator=[F(1),F(-n),F(n*(n+1),2)+(2 if n==2 else 0)]
    final=mul(local,denominator)
    Jhat=F((n-1)*n**3*(2609+1971*n+81*n*n),26244)+(F(16,27) if n==2 else 0)
    assert final==[Kn,-(n+F(26,9))*Kn,Jhat]
    assert 1+F(9,n)<6
    rows.append({'n':n,'Phi_peak_to_v2':list(map(str,phi)),
                 'qN_over_v6_to_v2':list(map(str,qn)),
                 'actual_ratio_over_av6_to_v2':list(map(str,final)),
                 'fifth_degree_sigma_gap':str(6-(1+F(9,n)))})
rates=[{'layer':h,'weighted_order':(9 if h==9 else (8 if h in (7,8) else 6))+9-h,
        'kind':'little-o' if h==6 else 'big-O'} for h in range(1,10)]
assert rates[5]=={'layer':6,'weighted_order':9,'kind':'little-o'}
assert all(row['weighted_order']>=9 for row in rates)
print(json.dumps({'started_at':start,'finished_at':datetime.datetime.now().astimezone().isoformat(),
                  'pid':os.getpid(),'python':sys.version,'rows':rows,
                  'N9_centered_ledger':rates,'initial_endpoint_order':10,
                  'scope':'Only the specified finite exact arithmetic and rate ledger. No geometric inequality or all-input theorem is established by this run.'},indent=2))
