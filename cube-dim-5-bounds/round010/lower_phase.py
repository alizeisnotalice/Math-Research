"""Exact phase-block area polynomials, independent of continuous-window cells."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import hashlib,json,os,time
P=Path(__file__).resolve().parents[1];start=time.monotonic();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
reps=[(u,v) for u,v in product(range(25),repeat=2) if (4*u+3*v)%25==0 and (-3*u+4*v)%25==0];assert len(reps)==25
strata={};checks=0
for m in [14,15,16]:
    hist={}
    for e1,e2 in product((0,1),repeat=2):
        h=Counter()
        for a,b in product(range(25),repeat=2):
            count=sum((u-a)%25<m+e1 and (v-b)%25<m+e2 for u,v in reps);h[count]+=1;checks+=25
        assert sum(h.values())==625;hist[(e1,e2)]=h
    counts=set().union(*(set(h) for h in hist.values()));polys={}
    for k in sorted(counts):
        h00,h01,h10,h11=[hist[e].get(k,0) for e in [(0,0),(0,1),(1,0),(1,1)]]
        c0,c1,c2=h00,-2*h00+h01+h10,h00-h01-h10+h11
        # Descending rho powers; normalization from scaled period area625 to25.
        polys[k]=[F(100*c2,25),F(10*c1-20*m*c2,25),F(c0-m*c1+m*m*c2,25)]
    assert [sum(p[j] for p in polys.values()) for j in range(3)]==[0,0,25]
    thresholds={n:[sum((p[j] for k,p in polys.items() if k>=n),F(0)) for j in range(3)] for n in range(1,17)}
    author=json.loads((P/f'lower_round020/stratum_m{m}.json').read_text())
    assert {str(k):list(map(str,p)) for k,p in polys.items()}==author['count_area_polynomials']
    assert {str(k):list(map(str,p)) for k,p in thresholds.items()}==author['threshold_area_polynomials']
    strata[m]={'integer_phase_histograms':{','.join(map(str,e)):dict(h) for e,h in hist.items()},
               'count_area_polynomials':{str(k):list(map(str,p)) for k,p in polys.items()},
               'threshold_area_polynomials':{str(k):list(map(str,p)) for k,p in thresholds.items()}}
cases=json.loads((P/'lower_round020/results.json').read_text())['all_case_exact_formulas']
mass=json.loads((P/'round010/lower_integer.json').read_text())['actual_mass_to_stratum'];verified=[]
for c in cases:
    assert c['K']*c['n']==c['N'] and mass[str(c['N'])]==c['m']
    poly=strata[c['m']]['threshold_area_polynomials'][str(c['n'])]
    assert poly==c['normalized_exact_area_polynomial'];verified.append({'case':[c['K'],c['n']],'N':c['N'],'polynomial_descending':poly})
assert len(verified)==16 and checks==187500 and time.monotonic()-start<1 and os.getpriority(os.PRIO_PROCESS,0)>=10
files=[Path(__file__),P/'round010/lower_phase_registration.md',P/'round010/lower_integer.json',P/'lower_round020/results.json']+[P/f'lower_round020/stratum_m{m}.json' for m in [14,15,16]]
out={'status':'PASS_independent_exact_continuous_radius_area_polynomials','strata':strata,'all_16_cases':verified,
     'physical_area_rule':'For normalized a*rho^2+b*rho+c and rho=w/(1+u), physical area=(a*w^2+b*w*(1+u)+c*(1+u)^2)/4;u^5=2,w^5=K*n.',
     'closed_boundary_note':'Phase endpoints have planar measure zero; independent integer extrema checker separately covers closed-window maxima/minima.',
     'scope':'Full planar area formulas only; no 5D union/parameter maximum/C5 result.',
     'independence':'Uses25 integer residues and four width-block frequency histograms, not producer endpoint arrangement, bitmasks, cells, or polynomial coefficients as input.',
     'source_sha256':{str(p.relative_to(P)):sha(p) for p in files},
     'runtime':{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'membership_checks':checks,'seconds':time.monotonic()-start}}
(P/'round010/lower_phase.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','runtime']},indent=2))
