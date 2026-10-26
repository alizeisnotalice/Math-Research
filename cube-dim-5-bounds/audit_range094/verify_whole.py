#!/usr/bin/env python3
"""094 all-node exact directed fixed-point envelope / integer-root audit.
No Arb, production import, new H, log, exp, asinh, jet or stencil evaluation.
Old M2 values are consumed from the locked accepted 088 full moment certificate.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,factorial
from collections import Counter
import hashlib,json,gzip,os,time,multiprocessing as mp,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-alt-cont/cube-dim-5-bounds')
R=BASE/'upper_range092';OLD=BASE/'upper_reenvelope088'
SCR=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_range092')
N=1<<160;K=1<<256;SHIFT=1<<96
CAT=None;TS=None;CODE=None;SCOPE=None

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def objsha(x):return hashlib.sha256(serial(x).encode()).hexdigest()
def readgz(p):return json.loads(gzip.decompress(Path(p).read_bytes()))
def save(p,x):
 p=Path(p);t=Path(str(p)+'.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def init():
 global CAT,TS
 CAT={F(k):F(v['upper']) for k,v in json.loads((R/'ratio_catalog.json').read_text()).items()}
 TS=list(map(F,json.loads((R/'registration.json').read_text())['time_endpoints']))
def directed(q):
 a=q.numerator*K;d=q.denominator
 return a//d,-(-a//d)
def ceildiv(a,b):return -(-a//b)
def bins(pair,kind):
 a,b=pair
 if kind=='floor':lo,hi=a//SHIFT,b//SHIFT
 else:lo,hi=ceildiv(a,SHIFT),ceildiv(b,SHIFT)
 return lo,hi

def audit_one(task):
 row,c,old=task;rank=row['rank']
 try:
  # Independent identity comparison: 092 metadata cannot relabel an 088 certificate.
  assert objsha(old)==row['source088_record_sha256']
  assert row['source088_rank']==old['rank']
  for name in ['key','level','indices','weight','M2','mu_infinity','h_infinity']:
   assert row[name]==old[name]==c[name],(name,rank)
  mapping={'old_self_L':'old_self_L','old_self_U':'new_self_U','old_rangeL_i':'rangeL_i','old_rangeU_i':'rangeU_i','old_M_i':'new_mean_i'}
  for k,v in mapping.items():assert row[k]==old[v]==c[k],(k,rank)
  assert row['original_raw_sha256']==old['source_evidence_sha256']==c['original_raw_sha256']
  assert c['rank']==rank and c['status']=='complete' and c['positive_bands']==73
  assert 0<row['old_rangeL_i']<=row['old_M_i']<=row['old_rangeU_i']<=3*N
  assert c['packet_sha256']==row['packet_sha256'] and c['new_H_calls']==c['old_exp_asinh_calls']==0
  # Read bytes of shared original endpoints; this checks copying, not their H production.
  rb=Path(row['original_raw_path']).read_bytes();assert hashlib.sha256(rb).hexdigest()==row['original_raw_sha256']
  raw=json.loads(gzip.decompress(rb));pb=Path(row['packet_path']).read_bytes();assert hashlib.sha256(pb).hexdigest()==row['packet_sha256']
  packet=json.loads(gzip.decompress(pb));assert packet['key']==row['key'] and packet['original_raw_sha256']==row['original_raw_sha256']
  assert packet['point']==raw['point'] and packet['endpoints']==raw['endpoints']
  assert raw['indices']==row['indices'] and raw['weight']==row['weight'] and raw.get('level',5)==row['level']
  point=list(map(F,packet['point']));assert len(point)==5 and point[0]==1
  level=row['level'];idx=row['indices'];offsets=[(F(i,1<<level)-x,F(i+1,1<<level)-x) for i,x in zip(idx,point[1:])]
  h=max(abs(t) for span in offsets for t in span);mu=max(abs((a+b)/2) for a,b in offsets)
  assert h==F(row['h_infinity']) and mu==F(row['mu_infinity'])
  M2=F(row['M2']);assert 0<=M2<=h*h # Full M2 identity consumed from accepted 088, not recomputed.
  ends=packet['endpoints'];assert [F(e['t']) for e in ends]==TS
  masses=[(F(e['H']['lower']),F(e['H']['upper'])) for e in ends];assert all(0<a<=b for a,b in masses)
  pairs=[[0,0],[0,0],[0,0]]
  for i,(a,(lo,up)) in enumerate(zip(TS,masses)):
   nlo,nup=masses[i+1] if i<72 else (F(0),F(0))
   ml,mh=lo-nup,up-nlo;assert 0<ml<=mh
   rat=CAT[h/a];fac=min(rat,1+mu/a+F(25,48)*M2*rat/(a*a))
   for pair,q in zip(pairs,[ml/rat,mh*rat,mh*fac]):
    fl,ce=directed(q);pair[0]+=fl;pair[1]+=ce
  qb=[bins(pairs[0],'floor'),bins(pairs[1],'ceil'),bins(pairs[2],'ceil')]
  # Fixed 256-bit sum intervals identify the exact 160-bit real-arithmetic bins.
  # A rare ambiguous bin falls back to exact rational sums (same finite node).
  fallback=False
  if any(a!=b for a,b in qb):
   fallback=True;sums=[F(0),F(0),F(0)]
   for i,(a,(lo,up)) in enumerate(zip(TS,masses)):
    nlo,nup=masses[i+1] if i<72 else (F(0),F(0));ml,mh=lo-nup,up-nlo;rat=CAT[h/a];fac=min(rat,1+mu/a+F(25,48)*M2*rat/(a*a))
    for j,v in enumerate([ml/rat,mh*rat,mh*fac]):sums[j]+=v
   vals=[(sums[0]*N).__floor__(),(sums[1]*N).__ceil__(),(sums[2]*N).__ceil__()]
  else:vals=[a for a,b in qb]
  el=max(row['old_rangeL_i'],vals[0]);eu=min(row['old_rangeU_i'],3*N,vals[1]);em=min(row['old_M_i'],eu,vals[2])
  li,ui,mi=c['rangeL_i'],c['rangeU_i'],c['mean_i']
  assert 0<li<=mi<=ui
  # Arb may outward-round by one extra quantum; any discrepancy is recorded.
  assert row['old_rangeL_i']<=li<=el and eu<=ui<=min(row['old_rangeU_i'],3*N)
  assert min(row['old_M_i'],ui,vals[2])<=mi<=min(row['old_M_i'],ui)
  assert li<=min(row['old_M_i'],ui,vals[2]),'range/mean incompatible'
  a=isqrt(isqrt(li**5));b=isqrt(isqrt(ui**5))
  assert a**4<=li**5<(a+1)**4 and b**4<=ui**5<(b+1)**4
  au=a if a**4==li**5 else a+1;bu=b if b**4==ui**5 else b+1
  il,iu=c['new_power_lower_i'],c['new_chord_upper_i'];D=1<<(40+4*level)
  assert il*D<=a
  if ui==li:assert iu*D>=au
  else:assert iu*(ui-li)*D>=(ui-mi)*au+(mi-li)*bu
  assert c['new_self_L']==max(row['old_self_L'],row['weight']*il)
  assert c['new_self_U']==min(row['old_self_U'],row['weight']*iu)
  assert c['new_self_L']<=c['new_self_U']
  return {'rank':rank,'key':row['key'],'status':'pass','bands':73,'exact_bins_match':(li,ui,mi)==(el,eu,em),'fallback':fallback,'raw_bytes':len(rb),'packet_bytes':len(pb)}
 except BaseException as e:
  import traceback
  return {'rank':rank,'key':row['key'],'status':'failed','error':repr(e),'traceback':traceback.format_exc()}

def verify_manifest(directory,expected):
 assert sha(directory/'manifest.json')==expected
 data=json.loads((directory/'manifest.json').read_text());items=data.get('files',data.get('artifacts'))
 if isinstance(items,dict):items=[{'path':k,'sha256':v['sha256'] if isinstance(v,dict) else v} for k,v in items.items()]
 assert isinstance(items,list),directory
 for r in items:
  name=r.get('path',r.get('name'));assert sha(directory/name)==r['sha256'],(directory,name)
 return len(items)

def hierarchy(state):
 nodes=state['nodes'];seen=set();visiting=set()
 def dfs(k):
  if k in seen:return nodes[k]['effL'],nodes[k]['effU']
  assert k not in visiting;visiting.add(k);n=nodes[k];L,U=n['selfL'],n['selfU']
  if n['children']:
   child=[dfs(c) for c in n['children']];L=max(L,sum(a for a,b in child));U=min(U,sum(b for a,b in child))
  assert L<=U;n['effL'],n['effU']=L,U;visiting.remove(k);seen.add(k);return L,U
 sums=[dfs(k) for k in state['roots']];assert len(seen)==56139
 return sum(a for a,b in sums),sum(b for a,b in sums)

def verify_budget(budget,ints,reg):
 lo,up=map(lambda x:F(x,N),ints);B,I,C=map(F,[budget['B_upper'],budget['I_upper'],budget['cost_upper']]);E=F(reg['E076_upper']);g=F(reg['valid_g'])
 assert F(budget['LH_truncated_integral'])==lo and F(budget['UH_truncated_integral'])==up
 assert E==F(budget['fixed_E076_upper']) and g==F(budget['valid_g'])==F(1,10000)
 assert I>0 and C>0
 assert B>E and (B-E)**5>=up**4 and I**4>=B**5 and C**4>=(5*B/g)**5
 gate=5**5*B**5<32**4*g**5;improve=I<F(reg['old_complete_I_cap'])
 assert gate==budget['cost_strictly_below32'] and improve==budget['I_strictly_below_old_complete_cap']
 assert budget['no_tube_added'] and budget['physical_frontier_count']==55819
 return {'LH':str(lo),'UH':str(up),'I_upper':str(I),'B_upper':str(B),'cost_upper':str(C),'strictly_improves_old_complete_cap':improve,'strict32_gate':gate}

def main():
 started=time.monotonic();scope=json.loads((P/'whole_scope.json').read_text());assert scope['code_sha256']==sha(__file__)
 assert not (P/'whole_results.json').exists(),'completed audit refuses replay'
 for r in scope['source_files']:assert sha(r['path'])==r['sha256'],r['path']
 manifest_counts={}
 for d,h in scope['manifests'].items():manifest_counts[d]=verify_manifest(Path(d),h)
 assert json.loads((P/'catalog_results.json').read_text())['status']=='pass'
 reg=json.loads((R/'registration.json').read_text())
 Epath=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/upper_shorttime076/results.json')
 Ecase=next(c for c in json.loads(Epath.read_text())['cases'] if c['epsilon']=='1/512')
 assert reg['E076_upper']==Ecase['E_norm']['upper']
 gp=OLD/'globalg_acceptance012.json';assert sha(gp)==reg['g_acceptance_sha256']
 assert reg['valid_g']==json.loads(gp.read_text())['g']=='1/10000'
 assert reg['old_complete_I_cap']==json.loads((OLD/'summary.json').read_text())['I_upper_simple']=='70301737/1000000000000'
 rows=readgz(R/'nodes.json.gz');assert len(rows)==56139 and [r['rank'] for r in rows]==list(range(56139))
 state=readgz(SCR/'initial_state.json.gz');assert len(state['nodes'])==56139 and len(state['frontier'])==55819
 assert sha(SCR/'initial_state.json.gz')==reg['initial_state_sha256']==json.loads((OLD/'summary.json').read_text())['state_sha256']
 finalstate=readgz(SCR/'state.json.gz');assert state['roots']==finalstate['roots'] and state['frontier']==finalstate['frontier']
 immutable=['key','level','indices','weight','children','parent','source','rootid']
 assert set(state['nodes'])==set(finalstate['nodes'])=={r['key'] for r in rows}
 for k,n in state['nodes'].items():
  for field in immutable:assert n.get(field)==finalstate['nodes'][k].get(field)
  assert n['key']==k and n['indices']==sorted(n['indices'],reverse=True)
  assert n['weight']==factorial(4)//__import__('functools').reduce(lambda a,b:a*factorial(b),Counter(n['indices']).values(),1)
  if n['children']:
   assert sum(state['nodes'][c]['weight'] for c in n['children'])==16*n['weight']
   for c in n['children']:assert state['nodes'][c]['parent']==k and state['nodes'][c]['level']==n['level']+1
 assert sum(F(state['nodes'][k]['weight'],1<<(4*state['nodes'][k]['level'])) for k in state['frontier'])==1
 assert len(state['roots'])==52360 and sum(bool(n['children']) for n in state['nodes'].values())==320
 old={};cells=[];blockdata=[]
 for bid in range(14):
  old.update({r['rank']:r for r in readgz(OLD/f'block{bid:04d}_cells.json.gz')})
  cs=readgz(R/f'block{bid:04d}_cells.json.gz');br=json.loads((R/f'block{bid:04d}_registration.json').read_text());rc=json.loads((R/f'block{bid:04d}_receipt.json').read_text())
  start,end=4096*bid,min(56139,4096*(bid+1));assert br['start']==start and br['end']==end and br['count']==end-start
  assert br['selection_sha256']==objsha(rows[start:end]) and br['keys']==[r['key'] for r in rows[start:end]]
  assert [c['rank'] for c in cs]==list(range(start,end)) and rc['count']==end-start and rc['status']=='complete'
  assert rc['registration_sha256']==sha(R/f'block{bid:04d}_registration.json') and rc['certificates_sha256']==sha(R/f'block{bid:04d}_cells.json.gz')
  assert br['source_engine_sha256']==reg['engine_sha256'] and br['source_registration_sha256']==sha(R/'registration.json')
  assert br['parent_receipt_sha256']==(sha(R/f'block{bid-1:04d}_receipt.json') if bid else None)
  cells.extend(cs);blockdata.append((br,rc))
 assert len(old)==len(cells)==56139
 for row in rows:
  n=state['nodes'][row['key']];assert n['selfL']==row['old_self_L'] and n['selfU']==row['old_self_U'] and (row['key'] in state['frontier'])==row['physical_frontier']
  for field in ['indices','level','weight']:assert n[field]==row[field]
 # Each block is a durable finite checkpoint. Resume skips blocks already audited with identical code/scope.
 blockstats=[]
 for bid in range(14):
  cp=P/f'whole_block{bid:04d}.json';start,end=4096*bid,min(56139,4096*(bid+1));signature={'scope_sha256':sha(P/'whole_scope.json'),'code_sha256':sha(__file__),'block_cells_sha256':sha(R/f'block{bid:04d}_cells.json.gz')}
  if cp.exists():
   stats=json.loads(cp.read_text());assert stats['signature']==signature and stats['status']=='pass';blockstats.append(stats);continue
  t=time.monotonic();tasks=[(rows[i],cells[i],old[rows[i]['source088_rank']]) for i in range(start,end)]
  with mp.get_context('spawn').Pool(8,initializer=init) as pool:
   out=list(pool.imap_unordered(audit_one,tasks,chunksize=16))
  failed=[x for x in out if x['status']!='pass']
  if failed:save(P/'failure.json',{'block':bid,'failures':failed});raise AssertionError(f'{len(failed)} failed in {bid}')
  stats={'signature':signature,'block':bid,'status':'pass','nodes':len(out),'bands':sum(x['bands'] for x in out),'exact_bins_match':sum(x['exact_bins_match'] for x in out),'fallbacks':sum(x['fallback'] for x in out),'original_raw_bytes':sum(x['raw_bytes'] for x in out),'packet_bytes':sum(x['packet_bytes'] for x in out),'seconds':time.monotonic()-t}
  save(cp,stats);blockstats.append(stats);print('VERIFIED_BLOCK',bid,stats['nodes'],stats['seconds'],flush=True)
 # A separately ordered recursive proof reconstruction for every complete prefix.
 before=hierarchy(state);assert F(before[0],N)==F(reg['initial_LH']) and F(before[1],N)==F(reg['initial_UH'])
 prefixes=[]
 for bid,(br,rc) in enumerate(blockdata):
  assert verify_budget(br['before'],before,reg)==verify_budget(rc['before'],before,reg)
  start,end=br['start'],br['end']
  for c in cells[start:end]:
   n=state['nodes'][c['key']];assert (n['selfL'],n['selfU'])==(c['old_self_L'],c['old_self_U']);n['selfL'],n['selfU']=c['new_self_L'],c['new_self_U']
  # Explicitly check all yet-unprocessed nodes continue paying their old self interval.
  for row in rows[end:]:
   n=state['nodes'][row['key']];assert (n['selfL'],n['selfU'])==(row['old_self_L'],row['old_self_U'])
  after=hierarchy(state);assert after[0]>=before[0] and after[1]<=before[1]
  b=verify_budget(rc['after'],after,reg);b.update({'block':bid,'processed':end,'unprocessed_paid_old':56139-end});prefixes.append(b);before=after
 for k,n in state['nodes'].items():
  for field in ['selfL','selfU','effL','effU']:assert n[field]==finalstate['nodes'][k][field]
 final=prefixes[-1];I=F(final['I_upper']);cap=F(ceildiv(I.numerator*10**12,I.denominator),10**12)
 assert cap<F(reg['old_complete_I_cap']) and cap**4>=F(final['B_upper'])**5
 out={'status':'pass','pid':os.getpid(),'seconds_including_loading':time.monotonic()-started,'scope_sha256':sha(P/'whole_scope.json'),'code_sha256':sha(__file__),'nodes':56139,'bands':sum(b['bands'] for b in blockstats),'blocks':14,'physical_frontier':55819,'internal_parents':320,'physical_volume':'1','all_new_bins_exact_match':sum(b['exact_bins_match'] for b in blockstats)==56139,'rational_sum_fallbacks':sum(b['fallbacks'] for b in blockstats),'original_raw_bytes_checked':sum(b['original_raw_bytes'] for b in blockstats),'endpoint_packet_bytes_checked':sum(b['packet_bytes'] for b in blockstats),'manifest_file_counts':manifest_counts,'prefixes':prefixes,'I_upper_simple':str(cap),'shared_H_production':True,'shared_accepted_M2':True,'new_H_calls':0,'Arb_log_exp_asinh_calls':0,'new_g_certificate':False,'new_C5_claim':False}
 save(P/'whole_results.json',out);print('ALL_PASS',json.dumps({k:v for k,v in out.items() if k not in ['prefixes','manifest_file_counts']}),flush=True)
if __name__=='__main__':main()
