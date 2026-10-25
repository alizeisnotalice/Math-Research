#!/usr/bin/env python3
"""One registered serial endpoint stream, with unresolved clusters retained."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import hashlib,importlib.util,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;PREV=HERE.parent;PROJECT=PREV.parent
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(terms):
 out={}
 for d,c in terms:out[int(d)]=out.get(int(d),Q(0))+Q(c)
 return tuple(sorted((d,c) for d,c in out.items() if c))
def serial(terms):return [[d,str(c)] for d,c in terms]
def runtime(start):
 peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':peak if sys.platform=='darwin' else peak*1024,'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def run(label):
 start=time.monotonic();assert all(os.environ.get(k)=='1' for k in ENV);assert os.getpriority(os.PRIO_PROCESS,0)>=10
 reg=json.loads((HERE/f'registration_{label}.json').read_text());rawp=PREV/'raw_endpoint_catalog.json';hp=PREV/'endpoint_diagnostic.py';assert sha(rawp)==reg['raw_catalog_sha256'] and sha(hp)==reg['exact_helper_sha256']
 progressp=PREV/reg['prior_progress_path'];assert sha(progressp)==reg['prior_progress_sha256'];progress=json.loads(progressp.read_text())
 spec=importlib.util.spec_from_file_location('round017_one_stream_exact_helper',hp);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
 raw=json.loads(rawp.read_text());fs=[(Q(f['slope']),Q(f['offset_coefficient']),f['radical']) for f in raw['functions']];assert len(fs)==len(set(fs))==433
 selected=tuple(Q(t) for t in reg['selected_stream']);pending=[tuple(Q(t) for t in pair) for pair in progress.get('remaining_streams',progress.get('remaining_slope_pair_streams',[]))];assert selected in pending
 completed=[tuple(Q(t) for t in pair) for pair in reg['starting_completed_streams']];assert len(completed)==progress['completed_stream_count'] and selected not in completed
 groups={s:[(i,f) for i,f in enumerate(fs) if f[0]==s] for s in {f[0] for f in fs}}
 costs={pair:len(groups[pair[0]])*len(groups[pair[1]]) for pair in pending};assert costs[selected]==min(costs.values())==reg['expected_pairs']
 events={};pairs=0
 for (fi,f),(gi,g) in product(groups[selected[0]],groups[selected[1]]):
  pairs+=1
  if pairs>200000 or time.monotonic()-start>120:raise RuntimeError('save single-stream checkpoint')
  terms=key(((g[2],g[1]/(f[0]-g[0])),(f[2],-f[1]/(f[0]-g[0]))));bounds=h.rb(terms)
  if terms not in events:events[terms]={'bounds':bounds,'pairs':[]}
  events[terms]['pairs'].append([fi,gi])
 def classify(bounds):
  lo,hi=bounds
  return 'left' if hi<Q(3,2) else ('right' if lo>2 else ('in_window' if lo>=Q(3,2) and hi<=2 else 'boundary_uncertain'))
 def clusters(items):
  out=[]
  for terms,bounds in sorted(items,key=lambda t:(t[1][0],t[1][1])):
   if out and bounds[0]<=out[-1]['hi']:
    out[-1]['expressions'].append(terms);out[-1]['hi']=max(out[-1]['hi'],bounds[1])
   else:out.append({'lo':bounds[0],'hi':bounds[1],'expressions':[terms]})
  return out
 ownclusters=clusters((terms,row['bounds']) for terms,row in events.items());window={t:r for t,r in events.items() if classify(r['bounds']) not in ('left','right')}
 previous={}
 for row in reg['completed_event_inputs']:
  pp=PREV/row['path'];assert sha(pp)==row['sha256'];d=json.loads(pp.read_text())
  if 'window_events' in d:records=d['window_events']
  elif 'potential_events' in d:records=d['potential_events']
  else:
   records=[r for cluster in d.get('isolation_clusters',[]) for r in cluster['events'] if classify(tuple(Q(x) for x in r['interval'])) not in ('left','right')]
  for r in records:
   terms=key(r.get('formal_terms',r.get('expression')));bounds=tuple(Q(x) for x in r['interval'])
   if terms in previous:previous[terms]=(min(previous[terms][0],bounds[0]),max(previous[terms][1],bounds[1]))
   else:previous[terms]=bounds
 duplicates=sum(t in previous for t in window);merged=dict(previous)
 for t,r in window.items():
  b=r['bounds'];merged[t]=(min(merged[t][0],b[0]),max(merged[t][1],b[1])) if t in merged else b
 cross=clusters(merged.items());unresolved=[c for c in cross if len(c['expressions'])>1]
 winrows=[]
 for t,r in sorted(window.items(),key=lambda row:row[1]['bounds'][0]):
  power5=None
  if len(t)==1:power5=str(t[0][1]**5*t[0][0])
  winrows.append({'formal_terms':serial(t),'interval':[str(x) for x in r['bounds']],'classification':classify(r['bounds']),'fifth_power_if_single_radical':power5,'raw_function_pairs':r['pairs'],'formal_duplicate_of_completed':t in previous})
 old=json.loads((PREV/'endpoint_diagnostic_results.json').read_text());qp=PREV/'append_half2/qualification_K8_n1_results.json';assert sha(qp)==reg['plane_qualification_sha256'];qual=json.loads(qp.read_text());cp=PROJECT/'lower_round003/a2_equals_a3_results.json';assert sha(cp)==old['C_source_sha256']
 remaining=[p for p in pending if p!=selected];completed.append(selected)
 out={'status':'PASS_one_additional_stream_saved_with_cluster_boundaries','label':label,'starting_completed_stream_count':len(completed)-1,'raw_function_count':433,'raw_case_count':1024,'selected_stream':[str(t) for t in selected],'examined_pairs':pairs,'distinct_formal_events':len(events),'event_classification_counts':{c:sum(classify(r['bounds'])==c for r in events.values()) for c in ('left','in_window','right','boundary_uncertain')},
      'event_isolation_clusters':[{'interval':[str(c['lo']),str(c['hi'])],'formal_expressions':[serial(t) for t in c['expressions']],'resolved':len(c['expressions'])==1} for c in ownclusters],
      'within_stream_unresolved_clusters':sum(len(c['expressions'])>1 for c in ownclusters),'window_events':winrows,'formal_duplicates_against_completed':duplicates,'new_unique_formal_window_events':len(window)-duplicates,'completed_unique_formal_window_events_total':len(merged),
      'cross_completed_unresolved_clusters':[{'interval':[str(c['lo']),str(c['hi'])],'formal_expressions':[serial(t) for t in c['expressions']]} for c in unresolved],
      'ordering_scope':'only disjoint isolating clusters ordered; no order or inequality claimed inside unresolvedclusters',
      'completed_stream_count':len(completed),'total_stream_count':10,'completed_streams':[[str(a),str(b)] for a,b in completed],'remaining_streams':[[str(a),str(b)] for a,b in remaining],'remaining_stream_pair_counts':[{'slopes':[str(a),str(b)],'pair_count':costs[a,b]} for a,b in remaining],
      'original_pending_materialized_plane_packet_count':18,'analytically_qualified_capture_cases':qual['analytically_resolved_capture_cases'],'unresolved_capture_count':qual['remaining_unresolved_capture_count'],'unresolved_capture_cases':qual['remaining_unresolved_capture_cases'],'no_new_plane_packet_materialized':True,
      'capacity_diagnosis_unchanged':{'capacity_threshold_a_interval':old['capacity_threshold_a_interval'],'entire_window_capacity_exceeds_C':old['entire_window_capacity_exceeds_C'],'uniform_side12_average_upper':old['uniform_side12_average_upper'],'all_sides_ge12_excluded_uniformly':old['all_sides_ge12_excluded_uniformly']},
      'raw_catalog_sha256':sha(rawp),'prior_progress_sha256':sha(progressp),'C_source_sha256':sha(cp),'source_sha256':sha(Path(__file__)),'nodes':pairs,'runtime':runtime(start),'limits':{'seconds':120,'nodes':200000},'scope':'one registered stream only; no full eventexhaustion, parameterunion/grid, fixed-IDpartition or anglesearch'}
 (HERE/f'result_{label}.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 d=run(sys.argv[1]);print(json.dumps({k:v for k,v in d.items() if k not in ('event_isolation_clusters','window_events','capacity_diagnosis_unchanged')},indent=2))
