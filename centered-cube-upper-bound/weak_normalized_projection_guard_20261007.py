from pathlib import Path
from fractions import Fraction as F
import json,hashlib
p=Path(__file__).parent
reg=p/'weak_normalized_projection_registration_20261007.json';registration=json.loads(reg.read_text())
fixture=p/registration['fixture'];assert hashlib.sha256(fixture.read_bytes()).hexdigest()==registration['fixture_sha256']
checks=0;rows=[]
def check(v):
 global checks
 assert v;checks+=1
def mean(a):return sum(a,F(0))/len(a)
def cond(a,j):
 step=2**j;out=[]
 for k in range(0,len(a),step):out.extend([mean(a[k:k+step])]*step)
 return out
for old in json.loads(fixture.read_text())['rounds']:
 J=old['filtration_J'];f=list(map(F,old['complete_terminal_f']));N=len(f);W=mean(f)
 M={j:list(map(F,old['M'][str(j)])) for j in range(1,J+1)}
 peak=[max(M[j][x] for j in M) for x in range(N)]
 winner=[min(j for j in M if M[j][x]==peak[x]) for x in range(N)]
 start=checks
 for alpha in map(F,registration['thresholds']):
  event=[f[x]==0 and peak[x]>alpha for x in range(N)]
  g={j:[alpha/peak[x] if event[x] and winner[x]==j else F(0) for x in range(N)] for j in M}
  raw={j:[F(event[x] and winner[x]==j) for x in range(N)] for j in M}
  pi={j:cond(g[j],j) for j in M};piraw={j:cond(raw[j],j) for j in M}
  check(all(sum(g[j][x] for j in M)<=1 for x in range(N)))
  check(all(0<=pi[j][x]<=piraw[j][x]<=1 for j in M for x in range(N)))
  traffic=sum(mean([g[j][x]*M[j][x] for x in range(N)]) for j in M)
  check(traffic==alpha*F(sum(event),N))
  spi=[sum(pi[j][x] for j in M) for x in range(N)]
  check(traffic==mean([f[x]*spi[x] for x in range(N)]))
  surv=[F(1)]*N;stop=[F(0)]*N
  for j in range(J,0,-1):
   b=[surv[x]*pi[j][x] for x in range(N)]
   check(cond(b,j)==b)
   check(mean([b[x]*M[j][x] for x in range(N)])==mean([b[x]*f[x] for x in range(N)]))
   for x in range(N):stop[x]+=b[x];surv[x]*=1-pi[j][x]
  Ts=mean([f[x]*stop[x] for x in range(N)])
  D=mean([f[x]*(spi[x]-stop[x]) for x in range(N)])
  pair=mean([f[x]*(spi[x]**2-sum(pi[j][x]**2 for j in M))/2 for x in range(N)])
  rawsum=[sum(piraw[j][x] for j in M) for x in range(N)]
  rawpair=mean([f[x]*(rawsum[x]**2-sum(piraw[j][x]**2 for j in M))/2 for x in range(N)])
  check(all(stop[x]==1-surv[x] and 0<=stop[x]<=1 for x in range(N)))
  check(0<=Ts<=W and 0<=D<=pair<=rawpair)
  check(traffic==Ts+D)
  check(sum(mean([g[j][x]*f[x] for x in range(N)]) for j in M)==0)
  rows.append({'J':J,'threshold':str(alpha),'W':str(W),'weak_ratio':str(traffic/W),'stopped_ratio':str(Ts/W),'deficit_ratio':str(D/W),'pair_ratio':str(pair/W),'unweighted_pair_ratio':str(rawpair/W),'winner':winner,'event':event,'g':{str(j):list(map(str,g[j])) for j in g}})
 rows.append({'J':J,'checks':checks-start})
out={'status':'PASS','checks':checks,'rows':rows,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),'scope':registration['scope']}
(p/'weak_normalized_projection_results_20261007.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':checks,'rounds':[x for x in rows if 'checks'in x]}))
