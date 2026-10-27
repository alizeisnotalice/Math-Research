#!/usr/bin/env python3
"""Exactly one preregistered feasible-family probe; no search or optimizer."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib,json
p=Path(__file__).parent
plan=json.loads((p/'probe_plan.json').read_text())
assert plan['candidate_count']==1
n=8;t=Q(2);a=Q(8);r=Q(1);q=1+1/t;b=a*q;tau=Q(1,2**30)
anchor=(Q(7),)*n
centers=[tuple(map(Q,z)) for z in product((-7,7),repeat=n) if z!=(7,)*n]+[anchor]
assert len(centers)==256 and a>2*t and r<t
M=a**(-n)
assert all(max(map(abs,y))<a*r for y in centers)
assert all(abs(y_i)+a*r<16 for y in centers for y_i in y)
checks=0
energy=Q(0)
for j,yj in enumerate(centers):
    incoming=Q(0)
    for i in range(j):
        yi=centers[i]
        distance=max(abs(x-y) for x,y in zip(yi,yj))
        assert distance>2*r  # disjoint small cubes => complete source increments
        captured=distance<b*r
        assert not captured
        incoming+=M*captured
        if j<255:energy+=M*M*captured
        checks+=1
    assert incoming<=tau
assert all(max(abs(x-y) for x,y in zip(anchor,yj))>=b*r for yj in centers[:-1])
d=2**n-1;S=d*M;diag=d*M*M
calibrated=S*S/(2*d)-diag/2
stronger=S*S/(2*(d-1))-diag/2
assert energy==calibrated==0 and stronger>0
result={
'probe_count':1,'stopped_as_registered':True,'n':n,'a':str(a),'t':str(t),'b':str(b),'tau':str(tau),
'family_size':len(centers),'unpaid_size':d,'pair_checks':checks,
'weight_M':str(M),'unpaid_mass_S':str(S),'sum_squared_weights':str(diag),
'exact_energy':str(energy),'calibrated_rhs':str(calibrated),
'false_stronger_rhs_D254':str(stronger),
'input_witness':'f=2*1_{(-16,16)^8}, lambda=1; each selected a-star lies inside (-15,15)^8, hence has average 2.',
'all_checks_passed':True,
'candidate_conclusion':'Any coefficient strictly larger than 1/(2*(2^n-1)) fails with either exact diagonal correction or a^(-n)*S/2; this is a feasible-family energy counterexample, not a lower bound for C_n.',
'probe_plan_sha256':hashlib.sha256((p/'probe_plan.json').read_bytes()).hexdigest()
}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
