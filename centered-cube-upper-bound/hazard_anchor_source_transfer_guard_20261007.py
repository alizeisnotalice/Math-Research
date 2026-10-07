from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json,hashlib
P=Path(__file__).resolve().parent
REG=P/'hazard_anchor_source_transfer_registration_20261007.json'
reg=json.loads(REG.read_text()); checks=0

def check(x,msg):
 global checks
 checks+=1
 if not x: raise AssertionError(msg)

def d(x):
 x=F(x); return Decimal(x.numerator)/Decimal(x.denominator)

def weak(v):
 m=len(v)
 return max((F(i,m)*x for i,x in enumerate(sorted(v,reverse=True),1)),default=F(0))

rows=[]
with localcontext() as ctx:
 ctx.prec=reg['decimal_precision']
 tol=Decimal('1e-70')
 symbols={s:(Decimal(1)+d(F(s))).ln()/d(F(s)) for s in reg['symbol_lambda']}
 for n in reg['rounds']:
  start=checks; N=2*n
  maxres=Decimal(0); minhold=F(1); triples=0
  anchors=[F(j,4*n) for j in range(N+1)]
  check(anchors[0]==0 and anchors[-1]==F(1,2),'exact early coverage')
  for j in range(N):
   a=anchors[j]; b=anchors[j+1]
   for k in range(5):
    t=a+(b-a)*F(k,4); ca=1-a; ct=1-t
    hold=(ct/ca)**n
    minhold=min(minhold,hold)
    check(hold>=F(1,2),'original no-jump weight')
    check((t-a)/ca<=F(1,2*n),'relative bin width')
    check((1-F(1,2*n))**n>=F(1,2),'Bernoulli lower certificate')
    check(F(1,2)<=ct<=ca<=1,'early softness bounds')
    for s in (t,(1+t)/2,1-F(1,100*n*n),F(1)):
     check(t<=s<=1,'whole future, not just early')
     for ls,m in symbols.items():
      def kernel(end,begin):
       cb=d(1-begin); ce=d(1-end)
       return ((ce+(1-ce)*m)/(cb+(1-cb)*m))**n
      direct=kernel(s,a); split=kernel(s,t)*kernel(t,a)
      res=abs(direct-split); maxres=max(maxres,res)
      check(res<=tol,'original logarithmic cocycle')
      # Formula for the complete n-fold Bernoulli transfer, not just Fourier order.
      r=(s-t)/(1-t); g=m/(d(1-t)+d(t)*m)
      explicit=(1-d(r)+d(r)*g)**n
      check(abs(kernel(s,t)-explicit)<=tol,'Bernoulli symbol identity')
      triples+=1
  # Exact no-repeat entry probability plus positive remaining probability.
  check(F(0)<=1-minhold<=F(1,2),'positive remaining entry-label mass')
  # Full weak norms of step functions on Lebesgue [0,1), each cell length1/N.
  harm=sum((F(1,k) for k in range(1,N+1)),F(0))
  uniform=[F(1,N)]*N
  geom=[F(1,2**(j+1)) for j in range(N)]; z=sum(geom); geom=[x/z for x in geom]
  family=[]
  for name,weights,cyclic in [('aligned harmonic',uniform,False),('cyclic harmonic',uniform,True),('nonuniform cyclic',geom,True)]:
   sums=[F(0)]*N; norms=[]
   for j,w in enumerate(weights):
    column=[N*w/F(1+((k+j)%N) if cyclic else 1+k) for k in range(N)]
    aj=weak(column); check(aj==w,'exact constituent weak norm')
    norms.append(aj); sums=[a+b for a,b in zip(sums,column)]
   A=sum(norms); actual=weak(sums)
   entropy=-sum((d(w)*d(w).ln() for w in weights),Decimal(0))
   bound=d(A)*(Decimal(3)+2*Decimal(2).ln()+2*entropy)
   check(d(actual)<=bound,'entropy weak summation bound')
   check(entropy<=Decimal(N).ln()+tol,'entropy not input amplitude')
   # A rational, looser ceiling: ln(2N)<=ceil(log2(2N)).
   log2ceil=(2*N-1).bit_length()
   check(actual<=(3+2*log2ceil)*A,'rational logarithmic bound')
   if name=='cyclic harmonic':
    check(all(x==harm for x in sums),'cyclic harmonic constant sum')
    check(actual==harm and A==1,'logarithmic-loss pressure witness')
   if name=='aligned harmonic': check(actual==1,'aligned weak cost')
   family.append(dict(name=name,total_constituent_weak=str(A),sum_weak=str(actual),sum_weak_approx=float(actual),entropy=str(entropy)))
  spikes=[N*w for w in geom]
  check(weak(spikes)<=sum(geom),'disjoint spike total budget')
  check(sum(geom)==sum(uniform)==1,'one source mass partition')
  # A source at t=1/2 belongs to the last bin; no extra source receipt.
  bins=[[] for _ in range(N)]
  source_times=[F(0),F(1,2)]+[(anchors[j]+anchors[j+1])/2 for j in range(N)]
  for t in source_times:
   idx=min(N-1,int(t*4*n)); bins[idx].append(t)
  check(sum(len(x) for x in bins)==len(source_times),'single counting including endpoint')
  check(F(1,2) in bins[-1] and F(0) in bins[0],'endpoint assignments')
  rows.append(dict(n=n,bins=N,checks=checks-start,symbol_cases=triples,min_holding=str(minhold),max_cocycle_decimal_residual=str(maxres),families=family,scope='Original symbolic transition and exact measure-theoretic sum components, not spatial FIRST/CPGP'))
result=dict(status='PASS_EXACT_AND_DECIMAL',checks=checks,rounds=rows,registration_sha256=hashlib.sha256(REG.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='No proof of full ordered initial-data A(n)=polylog. No dimension fit. Decimal is not interval.')
(P/'hazard_anchor_source_transfer_guard_results_20261007.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(status=result['status'],checks=checks,rounds=[dict(n=x['n'],checks=x['checks'],maxres=x['max_cocycle_decimal_residual']) for x in rows])))
