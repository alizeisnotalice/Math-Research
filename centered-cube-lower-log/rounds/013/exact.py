# Import regenerates deterministic labels and the floating experiment first.
from fractions import Fraction as F
from bisect import bisect_left
from math import lcm
import json
import probe
checks=[];K=probe.K
for row in probe.rows:
 tau=F(str(row['tau']));delta=tau/K;t=delta/2
 L=F(29,25) if row['kind']=='generic' else (K+164+F(1,4)-tau/2)/K
 labels=(probe.configs[row['name']]*4-.5).astype(int)
 atoms=[F(8*j+2*int(y)+1,8*K) for j,y in enumerate(labels)]
 c1=[(a-L/2-t/2)%1 for a in atoms];c2=sorted((a+L/2+t/2)%1 for a in atoms)
 den=lcm(t.denominator,*[a.denominator for a in c1],*[a.denominator for a in c2]);tt=int(t*den)
 one=[int(a*den) for a in c1];two=[int(a*den) for a in c2];ext=[a-den for a in two]+two+[a+den for a in two]
 overlap=0
 for a in one:
  idx=bisect_left(ext,a)
  overlap+=max(0,tt-abs(a-ext[idx-1]))+max(0,tt-abs(a-ext[idx]))
 O=F(overlap,den);inc=K*delta+2*O-(K*delta)**2;ratio=inc/(K*delta)
 assert abs(float(ratio)-row['ratio'])<1e-9
 assert 1-tau<=ratio<=2-tau
 checks.append(dict(name=row['name'],tau=str(tau),kind=row['kind'],ratio=str(ratio),overlap=str(O)))
(probe.P/'exact.json').write_text(json.dumps(checks,indent=2));print('24 rational interval checks passed')
