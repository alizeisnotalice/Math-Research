#!/usr/bin/env python3
"""Independent exact barycentric fmpq certification of the negative triple.

The C++ diagnostic uses the DD triangular recurrence and binary64. Here the
same established generic Green formula is evaluated with exact rational
barycentric weights and the complete group. No diagnostic float is accepted.
"""
from pathlib import Path
from time import perf_counter
import hashlib,itertools,json,os
from flint import arb,ctx,fmpq

HERE=Path(__file__).resolve().parent
FAN=HERE.parents[2]/'upper/round003/direct_kernel/fan_integer_input.txt'
ctx.prec=384
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse():
 it=iter(FAN.read_text().split());out=[]
 for _ in range(int(next(it))):
  det=int(next(it));rs=[];alpha=fmpq(det)
  for _ in range(6):
   q=int(next(it));ray=tuple(int(next(it))for _ in range(6));rs.append((q,ray));alpha/=q
  out.append((alpha,rs))
 return out
def exact_sign_dd(nodes):
 z=sorted(nodes);assert all(v!=0 for v in z)
 assert all(z[j]!=z[j-1]for j in range(1,6))
 positive=[j for j,v in enumerate(z)if v>0]
 if len(positive)in[0,6]:return fmpq(0)
 if len(positive)<=3:chosen=positive;factor=2
 else:chosen=[j for j,v in enumerate(z)if v<0];factor=-2
 value=fmpq(0)
 for j in chosen:
  product=fmpq(1)
  for k in range(6):
   if k!=j:product*=z[j]-z[k]
  value+=fmpq(factor)/product
 return value
def evaluate(x,den,cones,started):
 tic=perf_counter();total=fmpq(0);terms=0;nonzero=0
 for perm in itertools.permutations(range(6)):
  for mask in range(64):
   y=[(-1 if mask>>j&1 else 1)*x[perm[j]]for j in range(6)]
   for alpha,rs in cones:
    nodes=[fmpq(sum(y[j]*ray[j]for j in range(6)),den*q)for q,ray in rs]
    value=exact_sign_dd(nodes);total+=alpha*value;terms+=1;nonzero+=int(value!=0)
  if perf_counter()-started>1200:raise TimeoutError('registered 20 minute cap')
 assert terms==967680
 q=total/128;g=arb(q.numerator)/arb(q.denominator)/arb.pi()**5
 return q,{'point_numerators':x,'denominator':den,'complete_group_terms':terms,'nonzero_terms':nonzero,
  'all_nodes_nonzero_and_distinct':True,'G_times_pi5_numerator':str(q.numerator),
  'G_times_pi5_denominator':str(q.denominator),'G_arb384':g.str(100,more=True),
  'elapsed_seconds':perf_counter()-tic}
def main():
 started=perf_counter();data=json.loads((HERE/'negative_triple.json').read_text());cones=parse()
 out={'status':'RUNNING_EXACT_THREE_POINT_CERTIFICATION','pid':os.getpid(),'arithmetic':'exact fmpq barycentric DD5 + Arb384 pi',
  'fan_sha256':sha(FAN),'triple_sha256':sha(HERE/'negative_triple.json'),
  'source_sha256':sha(Path(__file__)),'point_order':['minus','center','plus'],'points':[]}
 qs=[]
 for row,label in zip(data['points'],out['point_order']):
  q,result=evaluate(row['numerators'],row['denominator'],cones,started);qs.append(q);out['points'].append(result)
  out['elapsed_seconds']=perf_counter()-started
  (HERE/'exact_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
  print(label,result['G_arb384'][:65],'seconds',result['elapsed_seconds'],flush=True)
 gap=(qs[0]+qs[2])/2-qs[1];assert gap<0
 value=arb(gap.numerator)/arb(gap.denominator)/arb.pi()**5
 assert value<-arb(13)/4000000000
 out.update({'status':'PASS_STRICT_CHAMBER_CONVEXITY_COUNTEREXAMPLE','exact_Jensen_gap_times_pi5_numerator':str(gap.numerator),
  'exact_Jensen_gap_times_pi5_denominator':str(gap.denominator),'Jensen_gap_arb384':value.str(100,more=True),
  'Jensen_gap_strict_upper':'-13/4000000000','chamber_convexity_false':True,'new_C6_bound':False,
  'elapsed_seconds':perf_counter()-started,'stopped_on_first_negative_diagnostic':True})
 (HERE/'exact_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['Jensen_gap_arb384'],'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
