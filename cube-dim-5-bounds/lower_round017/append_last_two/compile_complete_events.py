#!/usr/bin/env python3
"""Compile saved stream certificates; do not rerun a stream or a union."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import hashlib,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;PREV=HERE.parent
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(t):
 out={}
 for d,c in t:out[int(d)]=out.get(int(d),Q(0))+Q(c)
 return tuple(sorted((d,c) for d,c in out.items() if c))
def main():
 start=time.monotonic();nodes=0;assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
 reg=json.loads((HERE/'compile_registration.json').read_text());rp=PREV/'raw_endpoint_catalog.json';assert sha(rp)==reg['raw_catalog_sha256'];raw=json.loads(rp.read_text());fs=[(Q(f['slope']),Q(f['offset_coefficient']),f['radical']) for f in raw['functions']];assert len(fs)==len(set(fs))==433 and len(raw['raw_cases'])==1024
 groups={s:[(i,f) for i,f in enumerate(fs) if f[0]==s] for s in {f[0] for f in fs}};expected=set(combinations(sorted(groups),2));coverage={};events={};parallel=[]
 for s,rows in sorted(groups.items()):
  ordered=sorted((f[1]**5*f[2],i) for i,f in rows);assert len({p for p,i in ordered})==len(ordered)
  parallel.append({'slope':str(s),'function_count':len(rows),'function_indices_in_constant_offset_order':[i for power,i in ordered],'ordering_proof':'the real fifthpower is strictly increasing; offset^5=coefficient^5*radicand is rational'})
  nodes+=len(rows)
 def inside(b):return b[1]>=Q(3,2) and b[0]<=2
 for stream_index,row in enumerate(reg['stream_inputs']):
  pp=PREV/row['path'];assert sha(pp)==row['sha256'];d=json.loads(pp.read_text());pair=tuple(Q(t) for t in d.get('selected_stream',d.get('slopes')));assert pair in expected and pair not in coverage
  count=d.get('examined_pairs',d.get('pairs'));assert count==len(groups[pair[0]])*len(groups[pair[1]])
  coverage[pair]={'stream_index':stream_index,'slopes':[str(x) for x in pair],'examined_pairs':count,'input_path':row['path'],'input_sha256':sha(pp),'exact_trace_sha256':d.get('all_pair_exact_trace_sha256'),'runtime':d.get('runtime')}
  if 'execution_path' in row:
   ep=PREV/row['execution_path'];assert sha(ep)==row['execution_sha256'];execution=json.loads(ep.read_text());coverage[pair]['runtime']=execution['runtime'];assert execution['one_stream_pairs']==count
  if 'window_events' in d:records=d['window_events']
  elif 'potential_events' in d:records=d['potential_events']
  else:records=[r for c in d.get('isolation_clusters',[]) for r in c['events'] if inside(tuple(Q(x) for x in r['interval']))]
  for record in records:
   nodes+=1
   if nodes>200000 or time.monotonic()-start>120:raise RuntimeError('saved-certificate compilation checkpoint')
   t=canon(record.get('formal_terms',record.get('expression')));b=tuple(Q(x) for x in record['interval']);assert inside(b)
   if t not in events:events[t]={'bounds':b,'streams':set()}
   else:events[t]['bounds']=(min(events[t]['bounds'][0],b[0]),max(events[t]['bounds'][1],b[1]))
   events[t]['streams'].add(stream_index)
 assert set(coverage)==expected and len(coverage)==10 and sum(r['examined_pairs'] for r in coverage.values())==68351
 ordered=sorted(events.items(),key=lambda row:(row[1]['bounds'][0],row[1]['bounds'][1]));clusters=[];out=[]
 for i,(t,row) in enumerate(ordered):
  b=row['bounds'];out.append({'event_index':i,'formal_terms':[[d,str(c)] for d,c in t],'interval':[str(x) for x in b],'supporting_stream_indices':sorted(row['streams'])})
  if clusters and b[0]<=clusters[-1]['hi']:clusters[-1]['indices'].append(i);clusters[-1]['hi']=max(clusters[-1]['hi'],b[1])
  else:clusters.append({'lo':b[0],'hi':b[1],'indices':[i]})
 unresolved=[c for c in clusters if len(c['indices'])>1];within=all(row['bounds'][0]>=Q(3,2) and row['bounds'][1]<=2 for t,row in ordered)
 last=json.loads((HERE/'checkpoint_10of10.json').read_text());assert last['completed_stream_count']==10 and not last['remaining_streams'] and len(events)==last['completed_unique_formal_window_events_total']
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 result={'status':'COMPLETE_raw_endpoint_event_coverage_only','a1_window':['3/2','2'],'finite_raw_family':'l1,l4,l5 in1..4;n in1..16;fixedtailrotation4/5,3/5','raw_case_count':1024,'raw_function_count':433,'slope_group_count':5,'completed_stream_count':10,'expected_stream_count':10,'all_different_slope_pairs_examined':68351,'parallel_offset_order_verified':parallel,'stream_coverage':sorted(coverage.values(),key=lambda row:row['stream_index']),'unique_formal_window_event_count':len(events),'ordered_events':out,'unresolved_clusters':[{'interval':[str(c['lo']),str(c['hi'])],'event_indices':c['indices']} for c in unresolved],'all_event_enclosures_inside_window':within,'order_of_distinct_event_positions_certified':not unresolved and within,
         'endpoint_permutations_on_each_component_not_evaluated':True,'remaining_event_streams':[],'original_materialized_plane_packet_count':18,'prior_analytic_capture':[[8,1]],'unresolved_capture_count':17,'unresolved_capture_cases':last['unresolved_capture_cases'],'raw_catalog_sha256':sha(rp),'source_sha256':sha(Path(__file__)),'runtime':{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'wall_seconds':time.monotonic()-start},'nodes':nodes,'scope':'complete a1 ordering-event coverage for only the current1024candidates/433endpoints; not a parameteroptimum, all-scaleparametercompleteness, unionvalue orC5 bound'}
 (HERE/'complete_events_10of10.json').write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('ordered_events','parallel_offset_order_verified','stream_coverage')},indent=2))
if __name__=='__main__':main()
