#!/usr/bin/env python3
"""One serial stream; retain window clusters and certify omitted outside pairs."""
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
def serial(t):return [[d,str(c)] for d,c in t]
def classify(b):
 lo,hi=b
 return 'left' if hi<Q(3,2) else ('right' if lo>2 else ('in_window' if lo>=Q(3,2) and hi<=2 else 'boundary_uncertain'))
def clusters(items):
 out=[]
 for terms,b in sorted(items,key=lambda t:(t[1][0],t[1][1])):
  if out and b[0]<=out[-1]['hi']:out[-1]['expressions'].append(terms);out[-1]['hi']=max(out[-1]['hi'],b[1])
  else:out.append({'lo':b[0],'hi':b[1],'expressions':[terms]})
 return out
def runtime(start):
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def run(label):
 start=time.monotonic();assert all(os.environ.get(k)=='1' for k in ENV);assert os.getpriority(os.PRIO_PROCESS,0)>=10
 reg=json.loads((HERE/f'registration_{label}.json').read_text());rawp=PREV/'raw_endpoint_catalog.json';hp=PREV/'endpoint_diagnostic.py';assert sha(rawp)==reg['raw_catalog_sha256'] and sha(hp)==reg['exact_helper_sha256']
 progressp=PREV/reg['prior_progress_path'];assert sha(progressp)==reg['prior_progress_sha256'];progress=json.loads(progressp.read_text());raw=json.loads(rawp.read_text())
 spec=importlib.util.spec_from_file_location('round017_low_resource_exact',hp);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
 fs=[(Q(f['slope']),Q(f['offset_coefficient']),f['radical']) for f in raw['functions']];assert len(fs)==len(set(fs))==433
 selected=tuple(Q(x) for x in reg['selected_stream']);pending=[tuple(Q(x) for x in pair) for pair in progress['remaining_streams']];completed=[tuple(Q(x) for x in pair) for pair in reg['starting_completed_streams']];assert len(completed)==progress['completed_stream_count'] and selected in pending and selected not in completed
 groups={s:[(i,f) for i,f in enumerate(fs) if f[0]==s] for s in {f[0] for f in fs}};costs={p:len(groups[p[0]])*len(groups[p[1]]) for p in pending};assert costs[selected]==min(costs.values())==reg['expected_pairs']
 counts={k:0 for k in ('left','in_window','right','boundary_uncertain')};events={};pairs=0;trace=hashlib.sha256();outside_bounds={}
 for (fi,f),(gi,g) in product(groups[selected[0]],groups[selected[1]]):
  pairs+=1
  if pairs>200000 or time.monotonic()-start>120:raise RuntimeError('save one stream checkpoint')
  terms=key(((g[2],g[1]/(f[0]-g[0])),(f[2],-f[1]/(f[0]-g[0]))));b=h.rb(terms);c=classify(b);counts[c]+=1
  trace.update((json.dumps([fi,gi,serial(terms),[str(x) for x in b],c],separators=(',',':'))+'\n').encode())
  if c in ('left','right'):
   if c in outside_bounds:outside_bounds[c]=(min(outside_bounds[c][0],b[0]),max(outside_bounds[c][1],b[1]))
   else:outside_bounds[c]=b
   continue
  if terms not in events:events[terms]={'bounds':b,'pairs':[]}
  events[terms]['pairs'].append([fi,gi])
 own=clusters((t,r['bounds']) for t,r in events.items());previous={}
 for row in reg['completed_event_inputs']:
  pp=PREV/row['path'];assert sha(pp)==row['sha256'];d=json.loads(pp.read_text())
  if 'window_events' in d:records=d['window_events']
  elif 'potential_events' in d:records=d['potential_events']
  else:records=[r for c in d.get('isolation_clusters',[]) for r in c['events'] if classify(tuple(Q(x) for x in r['interval'])) not in ('left','right')]
  for r in records:
   t=key(r.get('formal_terms',r.get('expression')));b=tuple(Q(x) for x in r['interval']);previous[t]=(min(previous[t][0],b[0]),max(previous[t][1],b[1])) if t in previous else b
 merged=dict(previous);window=[];links={}
 for c in raw['raw_cases']:
  for name in ('raw_lower','raw_upper'):
   a,b,d=c[name];fk=(Q(a),Q(b),int(d));links.setdefault(fk,[]).append(c)
 for t,r in sorted(events.items(),key=lambda x:x[1]['bounds'][0]):
  b=r['bounds'];merged[t]=(min(merged[t][0],b[0]),max(merged[t][1],b[1])) if t in merged else b
  linked={}
  for fi,gi in r['pairs']:
   for i in (fi,gi):
    for c in links.get(fs[i],[]):linked[(tuple(c['shape']),c['n'])]=c
  existing_zero=bool(linked) and all(c['planar_status']=='certified_zero_volume' for c in linked.values())
  window.append({'formal_terms':serial(t),'interval':[str(x) for x in b],'classification':classify(b),'fifth_power_if_single_radical':str(t[0][1]**5*t[0][0]) if len(t)==1 else None,'raw_function_pairs':r['pairs'],'formal_duplicate_of_completed':t in previous,'all_linked_nondomain_cases_have_existing_planar_zero':existing_zero,'linked_cases':[{'shape':c['shape'],'K':c['K'],'n':c['n'],'planar_status':c['planar_status'],'origin_result_sha256':c['origin_result_sha256']} for c in linked.values()] if t not in previous else []})
 cross=clusters(merged.items());unresolved=[c for c in cross if len(c['expressions'])>1];duplicates=sum(t in previous for t in events)
 old=json.loads((PREV/'endpoint_diagnostic_results.json').read_text());qp=PREV/'append_half2/qualification_K8_n1_results.json';assert sha(qp)==reg['plane_qualification_sha256'];qual=json.loads(qp.read_text());cp=PROJECT/'lower_round003/a2_equals_a3_results.json';assert sha(cp)==old['C_source_sha256']
 remaining=[p for p in pending if p!=selected];completed.append(selected)
 out={'status':'PASS_one_window_stream_with_outside_range_certificate','label':label,'starting_completed_stream_count':len(completed)-1,'raw_function_count':433,'raw_case_count':1024,'selected_stream':[str(x) for x in selected],'examined_pairs':pairs,'pair_classification_counts':counts,'outside_pair_interval_bounds':{k:[str(x) for x in b] for k,b in outside_bounds.items()},'all_pair_exact_trace_sha256':trace.hexdigest(),'outside_event_records_not_materialized':True,
      'window_events':window,'window_isolation_clusters':[{'interval':[str(c['lo']),str(c['hi'])],'formal_expressions':[serial(t) for t in c['expressions']],'resolved':len(c['expressions'])==1} for c in own],'within_window_unresolved_clusters':sum(len(c['expressions'])>1 for c in own),'formal_duplicates_against_completed':duplicates,'new_unique_formal_window_events':len(events)-duplicates,'new_unique_events_with_reusable_planar_zero':sum(not r['formal_duplicate_of_completed'] and r['all_linked_nondomain_cases_have_existing_planar_zero'] for r in window),'completed_unique_formal_window_events_total':len(merged),
      'cross_completed_unresolved_clusters':[{'interval':[str(c['lo']),str(c['hi'])],'formal_expressions':[serial(t) for t in c['expressions']]} for c in unresolved],'ordering_scope':'only separated windowclusters ordered; unresolvedclusters retained unchanged',
      'completed_stream_count':len(completed),'total_stream_count':10,'completed_streams':[[str(a),str(b)] for a,b in completed],'remaining_streams':[[str(a),str(b)] for a,b in remaining],'remaining_stream_pair_counts':[{'slopes':[str(a),str(b)],'pair_count':costs[a,b]} for a,b in remaining],
      'original_pending_materialized_plane_packet_count':18,'analytically_qualified_capture_cases':qual['analytically_resolved_capture_cases'],'unresolved_capture_count':qual['remaining_unresolved_capture_count'],'unresolved_capture_cases':qual['remaining_unresolved_capture_cases'],'no_new_plane_qualification_or_partition':True,
      'capacity_diagnosis_unchanged':{'capacity_threshold_a_interval':old['capacity_threshold_a_interval'],'entire_window_capacity_exceeds_C':old['entire_window_capacity_exceeds_C'],'uniform_side12_average_upper':old['uniform_side12_average_upper'],'all_sides_ge12_excluded_uniformly':old['all_sides_ge12_excluded_uniformly']},'raw_catalog_sha256':sha(rawp),'prior_progress_sha256':sha(progressp),'C_source_sha256':sha(cp),'source_sha256':sha(Path(__file__)),'nodes':pairs,'runtime':runtime(start),'limits':{'seconds':120,'nodes':200000},'scope':'one eventstream, no fullunion or grid; outsidepairs certified by exact intervals and tracehash rather than archivedindividually'}
 (HERE/f'result_{label}.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 d=run(sys.argv[1]);print(json.dumps({k:v for k,v in d.items() if k not in ('window_events','window_isolation_clusters','capacity_diagnosis_unchanged')},indent=2))
