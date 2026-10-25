#!/usr/bin/env python3
"""One new exact endpoint stream, keeping every earlier checkpoint intact."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import hashlib,importlib.util,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;PREV=HERE.parent;PROJECT=PREV.parent
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def runtime(start):
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def run():
 start=time.monotonic();assert all(os.environ.get(k)=='1' for k in ENV);assert os.getpriority(os.PRIO_PROCESS,0)>=10
 reg=json.loads((HERE/'registration_half2.json').read_text());rawp=PREV/'raw_endpoint_catalog.json';oldp=PREV/'endpoint_diagnostic_results.json';helperp=PREV/'endpoint_diagnostic.py'
 assert sha(rawp)==reg['raw_catalog_sha256'] and sha(oldp)==reg['previous_result_sha256'] and sha(helperp)==reg['exact_helper_sha256']
 spec=importlib.util.spec_from_file_location('round017_exact_helpers',helperp);helpers=importlib.util.module_from_spec(spec);spec.loader.exec_module(helpers)
 raw=json.loads(rawp.read_text());old=json.loads(oldp.read_text());fs=[(Q(f['slope']),Q(f['offset_coefficient']),f['radical']) for f in raw['functions']];assert len(fs)==len(set(fs))==433
 pending=[tuple(Q(x) for x in pair) for pair in old['remaining_slope_pair_streams']];selected=(Q(1,2),Q(2));assert selected in pending
 group={s:[f for f in fs if f[0]==s] for s in {f[0] for f in fs}}
 counts={pair:len(group[pair[0]])*len(group[pair[1]]) for pair in pending};assert counts[selected]==min(counts.values())==70
 rows=[];nodes=0
 for f,g in product(group[selected[0]],group[selected[1]]):
  nodes+=1
  if nodes>200000 or time.monotonic()-start>120:raise RuntimeError('save stream checkpoint')
  terms={}
  for d,c in ((g[2],g[1]),(f[2],-f[1])):terms[d]=terms.get(d,Q(0))+c
  terms=tuple(sorted((d,c/(f[0]-g[0])) for d,c in terms.items() if c));lo,hi=helpers.rb(terms)
  assert hi<=Q(4,3)<Q(3,2)
  rows.append({'formal_terms':[[d,str(c)] for d,c in terms],'interval':[str(lo),str(hi)],'window_class':'strictly_left_of_window'})
 rows.sort(key=lambda r:tuple(Q(x) for x in r['interval']));clusters=[]
 for row in rows:
  lo,hi=map(Q,row['interval'])
  if clusters and lo<=Q(clusters[-1]['interval'][1]):
   clusters[-1]['events'].append(row);clusters[-1]['interval'][1]=str(max(Q(clusters[-1]['interval'][1]),hi))
  else:clusters.append({'interval':[str(lo),str(hi)],'events':[row]})
 assert len(rows)==70 and all(len(c['events'])==1 for c in clusters)
 length4=[c for c in raw['raw_cases'] if c['shape'][0]==4];assert len(length4)==256
 for c in length4:assert helpers.roots(c['N'])[1]<=4<Q(3,2)*3
 cp=PROJECT/'lower_round003/a2_equals_a3_results.json';assert sha(cp)==old['C_source_sha256']
 remaining=[pair for pair in pending if pair!=selected]
 out={'status':'PASS_second_endpoint_stream_exactly_excluded_from_window','selected_stream':['1/2','2'],'raw_function_count':len(fs),'raw_case_count':len(raw['raw_cases']),'examined_pairs':nodes,'examined_events':len(rows),'potential_window_events':0,'isolation_clusters':clusters,'unresolved_clusters':0,
      'analytic_event_bound':'a_event=2*r/3<=4/3<3/2','uniform_first_length4_cases_excluded':256,
      'completed_stream_count':2,'total_stream_count':10,'previous_stream_not_rerun':['0','1'],'remaining_streams':[[str(a),str(b)] for a,b in remaining],
      'remaining_stream_pair_counts':[{'slopes':[str(a),str(b)],'pair_count':counts[a,b]} for a,b in remaining],
      'prior_pending_planar_cases':old['new_planar_K_n_cases_needed_for_window'],'prior_pending_planar_count':18,
      'capacity_diagnosis_unchanged':{'capacity_threshold_a_interval':old['capacity_threshold_a_interval'],'entire_window_capacity_exceeds_C':old['entire_window_capacity_exceeds_C'],'uniform_side12_average_upper':old['uniform_side12_average_upper'],'all_sides_ge12_excluded_uniformly':old['all_sides_ge12_excluded_uniformly']},
      'C_source_sha256':sha(cp),'raw_catalog_sha256':sha(rawp),'previous_result_sha256':sha(oldp),'previous_event_stream_sha256':sha(PREV/'event_stream_slopes0_1.json'),'source_sha256':sha(Path(__file__)),'nodes':nodes,'runtime':runtime(start),'limits':{'seconds':120,'nodes':200000},'scope':'only one newstream; no full eventexhaustion, union, parametergrid or newlowerbound'}
 (HERE/'stream_half2_results.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 d=run();print(json.dumps({k:v for k,v in d.items() if k not in ('isolation_clusters','capacity_diagnosis_unchanged')},indent=2))
