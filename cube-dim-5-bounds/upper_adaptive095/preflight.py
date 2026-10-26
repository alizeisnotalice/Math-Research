"""Finite1024-parent geometry/point-key preparation ONLY; no H or ratio call."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
from math import factorial,prod,lcm
import gzip,json,hashlib,time
H=Path(__file__).resolve().parent;O=H.parent/'upper_range092';A=H.parent/'upper_adapt082';B=H.parent/'upper_global080';S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive095');G=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089');GS=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_g105089')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def gz(p,x):p.write_bytes(gzip.compress(serial(x).encode(),mtime=0))
def key(L,i):return str(L)+':'+','.join(map(str,i))
def weight(i):return factorial(4)//prod(factorial(v) for v in Counter(i).values())
def qualify(L,i,verts):
 mid=[F(2*a+1,2**(L+1)) for a in i];eta=F(1,2**L*4096);cand=[mid]+[[x+eta*b**j for j,x in enumerate(mid)] for b in (F(1,8),F(1,16),F(1,32))];cand += [[x+F(v,3351328) for x,v in zip(mid,jit)] for jit in ((-1,-2,-4,-8),(-1,-3,-7,-11),(-2,-5,-11,-19),(-3,-7,-17,-29),(-2,-11,-23,-43),(-5,-13,-31,-47),(-7,-19,-37,-49))]
 for index,tail in enumerate(cand):
  p=[F(1)]+tail;D=lcm(*(x.denominator for x in p));nodes=[sum(int(x*D)*z for x,z in zip(p,v)) for v in verts]
  if len(set(nodes))==280 and all(F(a,2**L)<x<F(a+1,2**L) for a,x in zip(i,tail)):return list(map(str,p)),index
 raise ValueError(('no_finite_generic_candidate',L,i))
def main():
 began=time.monotonic();assert not (H/'children_plan.json.gz').exists();assert sha(O/'manifest.json')=='05cd878504800f4a66673d9a408666c6480b92a1db5828de01f6ba5ec0cd801c';sp=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_range092/state.json.gz');assert sha(sp)==json.loads((O/'summary.json').read_text())['state_sha256'];state=json.loads(gzip.decompress(sp.read_bytes()));nodes=state['nodes'];selected=sorted(state['frontier'],key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:1024];parents=[{'key':k,'level':nodes[k]['level'],'indices':nodes[k]['indices'],'weight':nodes[k]['weight'],'selfL':nodes[k]['selfL'],'selfU':nodes[k]['selfU'],'gap_units':nodes[k]['selfU']-nodes[k]['selfL']} for k in selected];children={};edges={}
 for p in parents:
  L=p['level']+1;cs=sorted({tuple(sorted((2*i+b for i,b in zip(p['indices'],bits)),reverse=True)) for bits in product((0,1),repeat=4)});assert sum(weight(c) for c in cs)==16*p['weight'];keys=[]
  for c in cs:
   k=key(L,c);assert k not in children and k not in nodes and tuple(a//2 for a in c)==tuple(p['indices']);children[k]={'key':k,'level':L,'indices':c,'weight':weight(c),'parent':p['key']};keys.append(k)
  edges[p['key']]=keys
 assert len(children)<=16384;pending_file=G/'batch0009/point_manifest.json';pending=json.loads(pending_file.read_text());pendingkeys={','.join(r['point']) for r in pending['points']};assert len(pendingkeys)==12
 reg={'status':'preparation_only_not_heat_execution','source092_manifest_sha256':sha(O/'manifest.json'),'source092_state_sha256':sha(sp),'parents':parents,'children_by_parent':edges,'parent_count':1024,'child_count':len(children),'candidate_rule':'midpoint +3 frozen046 eta beta candidates +7 fixed negative fallbacks, first exact280 eligible','pending9_manifest_path':str(pending_file),'pending9_manifest_sha256':sha(pending_file),'pending9_raw_point_keys':sorted(pendingkeys),'new_H_registered_for_execution':False,'main_execution_confirmation':False,'source_sha256':sha(__file__)};(H/'preflight_registration.json').write_text(json.dumps(reg,indent=2)+'\n');table=json.loads(gzip.decompress((B/'geometry017.json.gz').read_bytes()));points=set()
 for r in children.values():r['point'],r['candidate']=qualify(r['level'],r['indices'],table['vertices']);points.add(','.join(r['point']))
 pending_hit=sorted(points&pendingkeys);assert not pending_hit,'do not evaluate duplicated pending H64; wait for point qualification or reuse complete interval'
 # Metadata index of all previously completed H73 values; actual intervals read only if hit.
 irows=json.loads(gzip.decompress((O/'nodes.json.gz').read_bytes()));packetbykey={r['key']:(r['packet_path'],r['packet_sha256']) for r in irows};h73={}
 for r in json.loads(gzip.decompress((B/'roots.json.gz').read_bytes())):h73[','.join(r['point'])]=packetbykey[key(5,r['indices'])]
 for p in sorted(A.glob('batch*_jobs.json.gz')):
  for r in json.loads(gzip.decompress(p.read_bytes())):h73[','.join(r['point'])]=packetbykey[r['key']]
 assert len(h73)==56139;hit73=points&set(h73)
 # Accepted H64 base180040 plus frozen author/independently checked batch8 additions.
 greg=json.loads((G/'registration.json').read_text());base=Path(greg['source']);assert sha(base)==greg['source_sha256'];overlay_path=G/'batch0008/state.json';mf=GS/'batch0008_manifest.json';m=json.loads(mf.read_text());assert sha(overlay_path)==m['files'][str(overlay_path)];overlay=json.loads(overlay_path.read_text());print('overlay keys',list(overlay),flush=True);baseobj=json.loads(base.read_text());assert len(baseobj['cache'])==180040;cache=baseobj['cache'];new=overlay.get('new_cache',{});assert isinstance(new,dict);assert len(new)==4744;reused={};hits64=0
 for r in children.values():
  pk=','.join(r['point']);entry=new.get(pk,cache.get(pk));r['reuse_H73']=pk in h73;r['reuse_H64']=entry is not None
  if pk in h73:r['H73_packet_path'],r['H73_packet_sha256']=h73[pk]
  if entry is not None:
   assert entry['point']==r['point'];reused[pk]={'point':r['point'],'t':'1/64','H':{'lower':entry['values'][0][0],'upper':entry['values'][0][1]},'source':'089.batch8.new_cache' if pk in new else '089.base180040'};hits64+=1
 gz(H/'children_plan.json.gz',[children[k] for k in sorted(children)]);gz(H/'H64_reuse_values.json.gz',reused);summary={'status':'geometry_and_exact_point_registry_prepared_no_H','parents':1024,'canonical_children':len(children),'child_levels':dict(Counter(r['level'] for r in children.values())),'candidate_counts':dict(Counter(r['candidate'] for r in children.values())),'prior_H73_point_count':56139,'H73_exact_point_hits':len(hit73),'accepted_H64_base':180040,'accepted_H64_new_batch8':len(new),'H64_exact_point_hits':hits64,'pending9_point_count':12,'pending9_exact_point_time_intersection':pending_hit,'different_1over128_cache_kept_separate':True,'source089_base_sha256':greg['source_sha256'],'source089_batch8_manifest_sha256':sha(mf),'source089_batch8_overlay_sha256':sha(overlay_path),'source_pending9_manifest_sha256':sha(pending_file),'children_sha256':sha(H/'children_plan.json.gz'),'reuse_values_sha256':sha(H/'H64_reuse_values.json.gz'),'new_H_calls':0,'new_ratio_calls':0,'new_DD_calls':0,'main_execution_confirmation':False,'seconds':time.monotonic()-began};(H/'preflight_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
