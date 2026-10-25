#!/usr/bin/env python3
"""One bounded low-priority event-interface diagnostic; no new union run."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product,combinations
from pathlib import Path
from math import prod
import hashlib,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;DEN=2**80
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
@lru_cache(None)
def roots(n):
 lo,hi=0,8*DEN
 while lo+1<hi:
  m=(lo+hi)//2
  if m**5<=n*DEN**5:lo=m
  else:hi=m
 if lo**5==n*DEN**5:return Q(lo,DEN),Q(lo,DEN)
 return Q(lo,DEN),Q(lo+1,DEN)
def norm(n):
 c,d=1,1;p=2
 while p*p<=n:
  e=0
  while n%p==0:n//=p;e+=1
  c*=p**(e//5);d*=p**(e%5);p+=1
 if n>1:d*=n
 return c,d
def rb(terms):
 lo=hi=Q(0)
 for d,c in terms:
  a,b=roots(d)
  if c>=0:lo+=c*a;hi+=c*b
  else:lo+=c*b;hi+=c*a
 return lo,hi
def runtime(start):
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def main():
 assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
 start=time.monotonic();nodes=0
 def check():
  nonlocal nodes
  nodes+=1
  if nodes>200000:raise RuntimeError('save checkpoint before continuing event stream')
  if time.monotonic()-start>120:raise TimeoutError('save checkpoint before continuing event stream')
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();prior=json.loads((ROOT/'lower_round016/delivery007/case_evidence.json').read_text());case_map={(r['K'],r['n']):r for r in prior['cases']}
 ul,uh=roots(2);ss=((1+ul)/2,(1+uh)/2);functions={(Q(0),Q(0),1),(Q(1,2),Q(0),1)};raw=[];unknown=set();uniform=0;planarzeros=0
 for lens in product(range(1,5),repeat=3):
  l1,l4,l5=lens;K=prod(lens)
  jmin=-(l1-1)//2 if l1%2 else 1-l1//2;jmax=(l1-1)//2 if l1%2 else l1//2
  for n in range(1,17):
   check();N=K*n;c,d=norm(N);assert c**5*d==N;radius=(roots(N)[0]/2,roots(N)[1]/2)
   lower=(Q(jmax),Q(-c,2),d);upper=(Q(jmin),Q(c,2),d);functions.update((lower,upper))
   reasons=[];sidehi=2*radius[1]
   if sidehi<=(l1-1)*Q(3,2):reasons.append('uniform_first_axis_span')
   if sidehi<=(l4-1)*ss[0] or sidehi<=(l5-1)*ss[0]:reasons.append('fixed_tail_span')
   if reasons:uniform+=1
   old=case_map[K,n];status=old['status'];unevaluated=status.startswith('SKIP_') or status=='CERTIFIED_EMPTY_all_unrotated_shapes'
   planezero=not unevaluated and Q(old['area_interval'][1])==0
   if planezero:planarzeros+=1
   if unevaluated and not reasons:unknown.add((K,n))
   raw.append({'shape':lens,'K':K,'n':n,'N':N,'raw_lower':[str(lower[0]),str(lower[1]),lower[2]],'raw_upper':[str(upper[0]),str(upper[1]),upper[2]],'uniform_geometric_exclusion':reasons,'planar_status':'not_previously_evaluated' if unevaluated else ('certified_zero_volume' if planezero else 'evaluated_nonzero'),'origin_result_sha256':old['original_result_sha256']})
 assert len(raw)==1024
 fs=sorted(functions);catalog={'status':'COMPLETE_all1024_raw_endpoint_cases_before_pruning','window':['3/2','2'],'functions':[{'slope':str(a),'offset_coefficient':str(c),'radical':d} for a,c,d in fs],'raw_cases':raw,'not_previously_evaluated_planar_cases_needed':sorted(unknown),'source_sha256':sha(Path(__file__))};cp=HERE/'raw_endpoint_catalog.json';cp.write_text(json.dumps(catalog,separators=(',',':'))+'\n')
 # Exactly one registered slope-pair stream in this small batch.
 group0=[f for f in fs if f[0]==0];group1=[f for f in fs if f[0]==1];events={};pairs=0
 for f,g in product(group0,group1):
  check();pairs+=1;terms={}
  for d,c in ((g[2],g[1]),(f[2],-f[1])):terms[d]=terms.get(d,Q(0))+c
  terms=tuple(sorted((d,c/(f[0]-g[0])) for d,c in terms.items() if c));lo,hi=rb(terms)
  if hi<Q(3,2) or lo>2:continue
  events[terms]=(lo,hi)
 ordered=sorted(events.items(),key=lambda row:(row[1][0],row[1][1]));clusters=[]
 for expression,bounds in ordered:
  if clusters and bounds[0]<=clusters[-1]['hi']:
   clusters[-1]['expressions'].append(expression);clusters[-1]['hi']=max(clusters[-1]['hi'],bounds[1])
  else:clusters.append({'lo':bounds[0],'hi':bounds[1],'expressions':[expression]})
 cp0=ROOT/'lower_round003/a2_equals_a3_results.json';C=json.loads(cp0.read_text());cl,ch=map(Q,C['ratio_interval']);s4=(ss[0]**4,ss[1]**4);capacity_window=[(Q(3,2)*s4[0],Q(3,2)*s4[1]),(2*s4[0],2*s4[1])]
 cutoff=(1+Q(2,12))*(1+ss[1]/12)**2*(1+Q(14,5)*ss[1]/12+ss[1]**2/144)/(Q(3,2)*s4[0]);assert cutoff<1
 serialize=lambda terms:[[d,str(c)] for d,c in terms]
 stream={'status':'COMPLETE_one_slope_pair_stream_only','slopes':['0','1'],'pairs':pairs,'potential_events':[{'expression':serialize(t),'interval':[str(x) for x in b],'certainly_in_window':b[0]>=Q(3,2) and b[1]<=2} for t,b in ordered],'overlap_clusters':[{'interval':[str(c['lo']),str(c['hi'])],'expressions':[serialize(t) for t in c['expressions']]} for c in clusters if len(c['expressions'])>1],'raw_catalog_sha256':sha(cp),'source_sha256':sha(Path(__file__))};ep=HERE/'event_stream_slopes0_1.json';ep.write_text(json.dumps(stream,separators=(',',':'))+'\n')
 slopes=sorted({f[0] for f in fs});remaining=[(str(a),str(b)) for a,b in combinations(slopes,2) if (a,b)!=(Q(0),Q(1))]
 out={'status':'PASS_small_a1_event_interface_diagnostic_no_union','raw_case_count':len(raw),'unique_raw_functions':len(fs),'slope_group_counts':{str(a):sum(f[0]==a for f in fs) for a in slopes},'one_stream_pairs':pairs,'potential_events_in_this_stream':len(events),'unresolved_overlap_clusters':len(stream['overlap_clusters']),'remaining_slope_pair_streams':remaining,'uniform_geometrically_excluded_raw_cases':uniform,'raw_cases_with_reusable_planar_zero':planarzeros,'new_planar_K_n_cases_needed_for_window':sorted(unknown),'full_period_capacity_at_window_ends':[[str(t) for t in v] for v in capacity_window],'capacity_threshold_a_interval':[str(cl/s4[1]),str(ch/s4[0])],'entire_window_capacity_exceeds_C':capacity_window[0][0]>ch,'uniform_side12_average_upper':str(cutoff),'all_sides_ge12_excluded_uniformly':True,'C_source_sha256':sha(cp0),'catalog_sha256':sha(cp),'event_stream_sha256':sha(ep),'source_sha256':sha(Path(__file__)),'nodes':nodes,'runtime':runtime(start),'scope':'complete raw catalog,but only one event stream; finite familyratio is piecewiseaffine; no newunion,maximization or uniformallscale completenessclaim'}
 (HERE/'endpoint_diagnostic_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
