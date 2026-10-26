"""121 all-node readback, same accepted H inputs, new mean only; no producer imports."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,ceil
from collections import Counter
import json,gzip,hashlib,time,os,sys,multiprocessing as mp,copy
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
T=Path('/Users/zhengzhihao/Desktop/T');H=Path(__file__).resolve().parent
P=T/'Math-Research-cube5-path112-astra/cube-dim-5-bounds/upper_path112'
B=T/'Math-Research-cube5-alt-cont/cube-dim-5-bounds';S=T/'cube-finite-dim-scratch/upper_path112_astra'
N=1<<160;Q=1<<320;SHIFT=1<<160
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ser(x):return json.dumps(x,sort_keys=True,separators=(',',':'),default=str)
def read(p):return json.loads(gzip.decompress(Path(p).read_bytes())) if str(p).endswith('.gz') else json.loads(Path(p).read_text())
def write(p,x):
 p=Path(p);tmp=Path(str(p)+'.tmp');tmp.write_text(ser(x)+'\n');tmp.replace(p)
def key(L,ii):return str(L)+':'+','.join(map(str,ii))
def bind(rows,state):
 # Pure identity linkage to already accepted092/095/103 point/packet/bounds; no old numerical reevaluation.
 base={};pts={key(5,r['indices']):r['point'] for r in read(B/'upper_global080/roots.json.gz')}
 for j in range(2):pts.update({r['key']:r['point'] for r in read(B/'upper_adapt082'/f'batch{j:04d}_jobs.json.gz')})
 for r in read(B/'upper_range092/nodes.json.gz'):
  base[r['key']]={'point':pts[r['key']],'packet_path':r['packet_path'],'packet_sha256':r['packet_sha256']}
 for j in range(14):
  for r in read(B/'upper_range092'/f'block{j:04d}_cells.json.gz'):base[r['key']].update(old_L_i=r['rangeL_i'],old_U_i=r['rangeU_i'],old_M_i=r['mean_i'])
 for name in ('upper_adaptive095','upper_adaptive103'):
  jobs={r['key']:r for r in read(B/name/'jobs.json.gz')}
  for r in read(B/name/'children_certificates.json.gz'):
   k=r['key'];z=r['bounds'];base[k]={'point':jobs[k]['point'],'packet_path':str(T/'cube-finite-dim-scratch'/name/'evidence'/(k.replace(':','_').replace(',','-')+'.json.gz')),'packet_sha256':r['evidence_sha256'],'old_L_i':z['L_i'],'old_U_i':z['U_i'],'old_M_i':z['M_i']}
 geos=read(P/'geometries.json.gz');count=Counter();assert len(rows)==len(base)==len(state['nodes'])==97972
 for r in rows:
  k=r['key'];n=state['nodes'][k]
  for f,v in base[k].items():assert r[f]==v,(k,f)
  for f in ('level','indices','weight'):assert r[f]==n[f]
  assert r['old_selfL']==n['selfL'] and r['old_selfU']==n['selfU']
  assert F(r['point'][0])==1 and 0<r['old_L_i']<=r['old_M_i']<=r['old_U_i']
  off=[(F(i,1<<r['level'])-F(x),F(i+1,1<<r['level'])-F(x)) for i,x in zip(r['indices'],r['point'][1:])]
  assert off==[tuple(map(F,s)) for s in geos[r['geometry']]['offsets']];count[r['geometry']]+=1
 assert count=={g:v['nodes'] for g,v in geos.items()}
 return {'status':'pass','points_packets_bounds_bound':len(rows),'full_offset_identities':len(rows),'geometry_counts':dict(count),'old_geometry_and_range_math_rerun':False}
FAC=None;TIMES=None
def init(factors,times):
 global FAC,TIMES
 FAC={g:{a:F(v['factor_upper']) for a,v in fs.items()} for g,fs in factors.items()};TIMES=times
def one(task):
 row,out=task;k=row['key']
 try:
  assert out['key']==k and out['old_M_i']==row['old_M_i'] and out['old_selfU']==row['old_selfU'] and out['packet_sha256']==row['packet_sha256']
  assert sha(row['packet_path'])==row['packet_sha256'];ev=read(row['packet_path']);assert ev['point']==row['point'];ends=ev['endpoints'];assert [e['t'] for e in ends]==TIMES
  vals=[(F(e['H']['lower']),F(e['H']['upper'])) for e in ends];assert len(vals)==73 and all(0<l<=u for l,u in vals)
  low=high=0;terms=[]
  for j,(lo,up) in enumerate(vals):
   nl,nu=vals[j+1] if j<72 else (F(0),F(0));assert 0<lo-nu<=up-nl
   term=(up-nl)*FAC[row['geometry']][TIMES[j]];num=term.numerator*Q;den=term.denominator
   low+=num//den;high+=-(-num//den);terms.append(term)
  cl,ch=-(-low//SHIFT),-(-high//SHIFT);fallback=cl!=ch
  v=ceil(sum(terms)*N) if fallback else cl
  li,ui=row['old_L_i'],row['old_U_i'];mi=min(row['old_M_i'],ui,v);assert out['new_M_i']==mi and li<=mi<=ui
  # Fourth-root upper bound on 256-bit grid: exact dyadic cancellation before integer roots.
  powers=[li**5<<224,ui**5<<224];ru=[]
  for z in powers:
   a=isqrt(isqrt(z));assert a**4<=z<(a+1)**4;ru.append(a+(a**4<z))
  num=ru[0] if ui==li else (ui-mi)*ru[0]+(mi-li)*ru[1]
  den=1 if ui==li else ui-li
  chord_units=F(num,den*(1<<(96+4*row['level'])))
  independent=min(row['old_selfU'],row['weight']*ceil(chord_units))
  assert row['old_selfL']<=independent<=out['new_selfU']<=row['old_selfU']
  return {'key':k,'status':'pass','mean_exact_bin':not fallback,'root256_same_final_bin':independent==out['new_selfU'],'mean_improved':mi<row['old_M_i'],'self_improved':out['new_selfU']<row['old_selfU']}
 except BaseException as ex:
  import traceback
  return {'key':k,'status':'failed','error':repr(ex),'traceback':traceback.format_exc()}
def budget(b,up,E):
 g=F(21,200000);bv,iv,cv=map(F,[b['B_upper'],b['I_upper'],b['cost_upper']]);ev=F(b['E076']);u=F(up,N)
 assert ev==E and F(b['g105'])==g and F(b['UH'])==u
 assert bv>ev and iv>0 and cv>0 and (bv-ev)**5>=u**4 and iv**4>=bv**5 and cv**4>=(5*bv/g)**5
 margin=32**4*(g/5)**5-iv**4;assert F(b['exact_32_polynomial_margin'])==margin and b['strict32']==(margin>0) and b['no_tube_added']
 tl=F(b['S_lower'])+iv;assert tl>0 and tl**4<=32**4*(g/5)**5
 return {'UH':str(u),'I_upper':str(iv),'cost_upper':str(cv),'S_lower':b['S_lower'],'strict32':margin>0}
def hierarchy(initial,rows,cert,blocks,final):
 st=copy.deepcopy(initial);ns=st['nodes'];roots=st['roots'];st['path112']={'processed':[],'blocks':[]};done=set();trace=[]
 assert len(ns)==97972 and len(st['frontier'])==94580 and st['roots']==final['roots'] and st['frontier']==final['frontier']
 volume=sum((F(ns[k]['weight'],1<<(4*ns[k]['level'])) for k in st['frontier']),F(0));assert volume==1
 E=F(read(B/'upper_adaptive103/registration.json')['E076_upper']);oldU=sum(ns[k]['effU'] for k in roots);oldL=sum(ns[k]['effL'] for k in roots)
 def dfs(k):
  n=ns[k]
  if n['children']:
   v=[dfs(c) for c in n['children']];l=max(n['selfL'],sum(a for a,b in v));u=min(n['selfU'],sum(b for a,b in v))
  else:l,u=n['selfL'],n['selfU']
  assert l==n['effL'] and l<=u;n['effU']=u;return l,u
 for j,block in enumerate(blocks):
  name=block['name'];reg=read(P/(name+'_registration.json'));receipt=read(P/(name+'_receipt.json'));ks=reg['keys'];br=block['rows']
  assert receipt['status']=='complete' and reg['count']==receipt['count']==len(ks)==len(br)
  assert [r['key'] for r in br]==ks and not done.intersection(ks)
  expected=sorted(st['frontier'] if j==0 else set(ns)-done,key=lambda k:(-(ns[k]['selfU']-ns[k]['selfL']),ns[k]['level'],ns[k]['indices']))[:256 if j==0 else 4096]
  assert ks==expected
  assert reg['input_state_sha256']==hashlib.sha256(gzip.compress(ser(st).encode(),mtime=0)).hexdigest()
  assert reg['factor_catalog_sha256']==sha(P/'factor_catalog.json') and reg['nodes_sha256']==sha(P/'nodes.json.gz')
  budget(receipt['before'],oldU,E)
  for r in br:
   k=r['key'];assert r==cert[k] and r['new_selfU']<=ns[k]['selfU'];ns[k]['selfU']=r['new_selfU'];done.add(k);st['path112']['processed'].append(k)
  # Different algorithm: recursively rebuild every root at each complete block.
  vals=[dfs(k) for k in roots];lo=sum(a for a,b in vals);up=sum(b for a,b in vals);assert lo==oldL and up<=oldU
  for k in set(ns)-done:assert ns[k]['selfU']==initial['nodes'][k]['selfU']
  st['path112']['blocks'].append(name)
  assert receipt['state_sha256']==hashlib.sha256(gzip.compress(ser(st).encode(),mtime=0)).hexdigest()
  b=budget(receipt['after'],up,E);assert receipt['processed_total']==len(done) and receipt['strict_global_improvement']==(up<oldU)
  assert receipt['results_sha256']==block['sha256'] and receipt['range_L_U_unchanged'] and receipt['weighted_volume']=='1'
  trace.append({'name':name,'count':len(ks),'processed':len(done),'untouched_nodes_paid_old':len(ns)-len(done),'global_U_units':up,'volume':'1','budget':b});oldU=up
  write(H/'hierarchy_checkpoint.json',{'complete_blocks':j+1,'total':25,'processed':len(done)})
 assert st==final and len(done)==len(ns)==97972
 hierarchy=read(P/'numeric_hierarchy.json.gz')
 for k,n in hierarchy['nodes'].items():assert n==ns[k]
 assert hierarchy['frontier']==st['frontier'] and hierarchy['roots']==roots
 b=budget(hierarchy['budget'],oldU,E);cap=F(47619360,10**12);assert F(b['I_upper'])<=cap<F(51807888,10**12)
 oldfloor=read(H.parent/'audit_adaptive113/floor_results.json');TU=F(oldfloor['T_upper']);TL=F(oldfloor['T_lower']);IU=F(b['I_upper']);assert TL-IU<=TU-IU<0
 result={'status':'pass','nodes':len(ns),'physical_frontier':len(st['frontier']),'volume':'1','complete_prefixes':25,'all_lower_and_topology_unchanged':True,'all_unprocessed_self_paid_old':True,'LH':str(F(oldL,N)),'final_budget':b,'S_upper':str(TU-IU),'I_simple':str(cap),'cap_gain_over103':str(F(51807888,10**12)-cap),'trace':trace}
 write(H/'hierarchy_results.json',result);return result
def main():
 t=time.monotonic();assert not (H/'whole_results.json').exists();scope=read(H/'whole_registration_v2.json');assert sha(__file__)==scope['code_sha256']
 for r in scope['sources']:assert sha(r['path'])==r['sha256'],r['path']
 mf=read(P/'manifest.json');assert sha(P/'manifest.json')==scope['source112_manifest_sha256']
 for name,digest in mf['files'].items():assert sha(P/name)==digest,name
 assert read(H/'scalar_results.json')['status']=='pass' and read(H.parent/'audit_adaptive113/results.json')['status']=='pass'
 reg=read(P/'scalar_registration.json');initial=read(reg['source103_state_path']);assert sha(reg['source103_state_path'])==reg['source103_state_sha256'];final=read(S/'state.json.gz');assert sha(S/'state.json.gz')==scope['final_state_sha256']
 rows=read(P/'nodes.json.gz');by={r['key']:r for r in rows};cert={r['key']:r for r in read(P/'node_certificates.json.gz')};assert len(by)==len(cert)==97972 and set(by)==set(cert)
 linkage=bind(rows,initial);write(H/'linkage_results.json',linkage)
 ar=read(P/'author_verification.json');blocks=[];seen=set()
 for j,z in enumerate(ar['results']):
  assert sha(z['path'])==z['sha256'];rr=[json.loads(x) for x in Path(z['path']).read_text().splitlines()];assert len(rr)==z['count'];name='pilot' if j==0 else f'block{j:04d}'
  assert not seen.intersection(r['key'] for r in rr);seen.update(r['key'] for r in rr);assert all(r==cert[r['key']] for r in rr)
  blocks.append({'name':name,'rows':rr,'sha256':z['sha256']})
 assert len(blocks)==25 and seen==set(by)
 times=reg['time_endpoints'];assert times==read(B/'upper_adaptive103/registration.json')['time_endpoints'] and len(times)==73 and F(times[0])==F(1,512)
 factors=read(P/'factor_catalog.json');allout=[]
 with mp.get_context('spawn').Pool(8,initializer=init,initargs=(factors,times)) as pool:
  for block in blocks:
   name=block['name'];outpath=H/(name+'_checks.json')
   if outpath.exists():out=read(outpath)
   else:
    rr=list(pool.imap_unordered(one,[(by[r['key']],r) for r in block['rows']],chunksize=16));fail=[r for r in rr if r['status']!='pass']
    if fail:write(H/'failure.json',fail);raise AssertionError(fail[:1])
    out={'status':'pass','name':name,'count':len(rr),'keys_sha256':hashlib.sha256(ser(sorted(r['key'] for r in rr)).encode()).hexdigest(),'mean_exact_bins':sum(r['mean_exact_bin'] for r in rr),'root256_same_final_bin':sum(r['root256_same_final_bin'] for r in rr),'mean_improved':sum(r['mean_improved'] for r in rr),'self_improved':sum(r['self_improved'] for r in rr)};write(outpath,out)
   assert out['count']==len(block['rows']);allout.append(out);print('BAND_BLOCK_PASS',name,len(allout),flush=True)
 led=hierarchy(initial,by,cert,blocks,final)
 result={'status':'pass','pid':os.getpid(),'seconds':time.monotonic()-t,'nodes':97972,'positive_bands':97972*73,'blocks':25,'new_H_log_Arb_calls':0,'old_range_geometry_math_rerun':False,'all_full_offset_and_old_input_links':True,'mean_exact_bins':sum(z['mean_exact_bins'] for z in allout),'root256_same_final_bin':sum(z['root256_same_final_bin'] for z in allout),'mean_improved':sum(z['mean_improved'] for z in allout),'self_improved':sum(z['self_improved'] for z in allout),'final':{k:v for k,v in led.items() if k!='trace'},'code_sha256':sha(__file__),'registration_sha256':sha(H/'whole_registration_v2.json'),'shared_H_production':True,'new_C5_Cn_claim':False}
 write(H/'whole_results.json',result);print('WHOLE_PASS',ser(result),flush=True)
if __name__=='__main__':main()
