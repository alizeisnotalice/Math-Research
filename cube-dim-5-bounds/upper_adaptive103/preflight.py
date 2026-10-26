"""Finite2048-parent candidate continuation from frozen095; no newH."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
from math import factorial,prod,lcm
import gzip,json,hashlib,time,importlib.util
H=Path(__file__).resolve().parent;O=H.parent/'upper_adaptive095';B=H.parent/'upper_global080';A=H.parent/'upper_adapt082';R=H.parent/'upper_range092';S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive103');G=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_mixed098');GS=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_mixed098')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def gz(p,x):p.write_bytes(gzip.compress(serial(x).encode(),mtime=0))
def key(L,i):return str(L)+':'+','.join(map(str,i))
def weight(i):return factorial(4)//prod(factorial(v) for v in Counter(i).values())
spec=importlib.util.spec_from_file_location('pre095',O/'preflight.py');P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)
def main():
 t=time.monotonic();assert not (H/'children_plan.json.gz').exists();assert sha(O/'manifest.json')=='ad4e78f0e484d217c1d14bd13f585cd8f85240039b13915ffbd16b96b8d2f433';rc=json.loads((O/'run_receipt.json').read_text());sp=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive095/state.json.gz');assert sha(sp)==rc['state_sha256'];state=json.loads(gzip.decompress(sp.read_bytes()));nodes=state['nodes'];selected=sorted(state['frontier'],key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:2048];parents=[{**{x:nodes[k][x] for x in ('key','level','indices','weight','selfL','selfU')},'gap_units':nodes[k]['selfU']-nodes[k]['selfL']} for k in selected];children={};edges={}
 for p in parents:
  assert nodes[p['key']]['children'] is None;L=p['level']+1;cs=sorted({tuple(sorted((2*i+b for i,b in zip(p['indices'],bits)),reverse=True)) for bits in product((0,1),repeat=4)});assert sum(weight(c) for c in cs)==16*p['weight'];edges[p['key']]=[]
  for ii in cs:
   k=key(L,ii);assert k not in children and k not in nodes and tuple(i//2 for i in ii)==tuple(p['indices']);children[k]={'key':k,'level':L,'indices':ii,'weight':weight(ii),'parent':p['key']};edges[p['key']].append(k)
 assert len(children)<=32768;reg={'status':'finite_candidate_continuation_registered_before_H','source095_manifest_sha256':sha(O/'manifest.json'),'source095_state_sha256':sha(sp),'source095_pending_independent102':True,'parents':parents,'children_by_parent':edges,'parent_count':2048,'child_count':len(children),'candidate_rule':'same frozen095 finite046 family exact280 eligibility;true fullbox','source_sha256':sha(__file__),'new_H_calls':0};(H/'preflight_registration.json').write_text(json.dumps(reg,indent=2)+'\n');table=json.loads(gzip.decompress((B/'geometry017.json.gz').read_bytes()));points=set()
 for r in children.values():r['point'],r['candidate']=P.qualify(r['level'],r['indices'],table['vertices']);points.add(','.join(r['point']))
 packet={r['key']:(r['packet_path'],r['packet_sha256']) for r in json.loads(gzip.decompress((R/'nodes.json.gz').read_bytes()))};h73={}
 for r in json.loads(gzip.decompress((B/'roots.json.gz').read_bytes())):h73[','.join(r['point'])]=packet[key(5,r['indices'])]
 for f in sorted(A.glob('batch*_jobs.json.gz')):
  for r in json.loads(gzip.decompress(f.read_bytes())):h73[','.join(r['point'])]=packet[r['key']]
 cert={r['key']:r for r in json.loads(gzip.decompress((O/'children_certificates.json.gz').read_bytes()))}
 for r in json.loads(gzip.decompress((O/'jobs.json.gz').read_bytes())):
  k=r['key'];path=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive095/evidence')/(k.replace(':','_').replace(',','-')+'.json.gz');h73[','.join(r['point'])]=(str(path),cert[k]['evidence_sha256'])
 assert len(h73)==69545;mf=GS/'production_manifest.json';assert sha(mf)=='b318e998a41d3b8d78ab663e8a457244e1b87b03c380b05f7c27664fdd0f557d';m=json.loads(mf.read_text());gf=G/'final_state.json';assert sha(gf)==m['external_files'][str(gf)];cache=json.loads(gf.read_text())['cache_by_epsilon'];assert len(cache['1/64'])==184796 and len(cache['1/128'])==194;reuse={};hits=Counter()
 for r in children.values():
  pk=','.join(r['point']);r['reuse_H73']=pk in h73;r['reuse_H64']=pk in cache['1/64']
  if pk in h73:r['H73_packet_path'],r['H73_packet_sha256']=h73[pk];hits['H73']+=1
  for ts,tab in cache.items():
   if pk in tab:
    entry=tab[pk];assert entry['point']==r['point'];reuse.setdefault(pk,{})[ts]={'point':r['point'],'t':ts,'H':{'lower':entry['values'][0][0],'upper':entry['values'][0][1]},'source':'098.frozen_complete_same_model_cache','source_state_sha256':sha(gf)};hits[ts]+=1
 gz(H/'children_plan.json.gz',[children[k] for k in sorted(children)]);gz(H/'reuse_values.json.gz',reuse);summary={'status':'preflight_complete_no_H','parents':2048,'canonical_children':len(children),'child_levels':dict(Counter(r['level'] for r in children.values())),'candidate_counts':dict(Counter(r['candidate'] for r in children.values())),'prior_H73_point_count':69545,'exact_cache_hits':dict(hits),'cache_H64_count':184796,'cache_H128_count':194,'source098_production_manifest_sha256':sha(mf),'source098_final_state_sha256':sha(gf),'source095_manifest_sha256':sha(O/'manifest.json'),'source095_state_sha256':sha(sp),'children_sha256':sha(H/'children_plan.json.gz'),'reuse_values_sha256':sha(H/'reuse_values.json.gz'),'source_preflight_registration_sha256':sha(H/'preflight_registration.json'),'new_H_calls':0,'seconds':time.monotonic()-t};(H/'preflight_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
