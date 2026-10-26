"""Finite zero-H path112 consumer. Frozen old sources stay read-only.
All new arithmetic is rational or integer-root outward arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,ceil
import json,gzip,hashlib,os,time,sys
import pathmean as P
import secant as Q
import ratio091 as R
H=Path(__file__).resolve().parent
S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_path112_astra')
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-alt-cont/cube-dim-5-bounds/upper_adaptive103')
N=1<<160
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ser(x):return json.dumps(x,sort_keys=True,separators=(',',':'),default=str)
def read(p):return json.loads(gzip.decompress(Path(p).read_bytes())) if str(p).endswith('.gz') else json.loads(Path(p).read_text())
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_bytes(gzip.compress(ser(x).encode(),mtime=0) if p.suffix=='.gz' else (ser(x)+'\n').encode());os.replace(t,p)
def root(q,k,bits=224,up=True):
 d=1<<bits;n=q.numerator*d**k//q.denominator
 if k==4:v=isqrt(isqrt(n))
 elif k==5:
  v=1<<((n.bit_length()+4)//5) if n else 0
  while v:
   z=(4*v+n//v**4)//5
   if z>=v:break
   v=z
  while (v+1)**5<=n:v+=1
  while v**5>n:v-=1
 else:raise ValueError(k)
 lo=F(v,d);assert lo**k<=q<(lo+F(1,d))**k
 return lo if not up or lo**k==q else lo+F(1,d)
def budget(up):
 reg=read(SRC/'registration.json');E=F(reg['E076_upper']);g=F(21,200000);u=F(up,N);B=root(u**4,5,384)+E;I=root(B**5,4,384);T=32*root((g/5)**5,4,384,False);C=root((5*B/g)**5,4,384);margin=32**4*(g/5)**5-I**4
 return {'UH':str(u),'I_upper':str(I),'I_display':float(I),'E076':str(E),'g105':str(g),'B_upper':str(B),'cost_upper':str(C),'cost_display':float(C),'S_lower':str(T-I),'S_lower_display':float(T-I),'exact_32_polynomial_margin':str(margin),'strict32':margin>0,'no_tube_added':True}
def scalars():
 reg=read(H/'scalar_registration.json');cat=read(H/'inherited_ratio_catalog.json');out=H/'ratio_catalog.json';journal=S/'new_ratio_rows.jsonl';S.mkdir(parents=True,exist_ok=True)
 assert sha(H/'inherited_ratio_catalog.json')==reg['inherited_catalog_sha256']
 existing={}
 if journal.exists():
  for line in journal.read_text().splitlines():
   z=json.loads(line);assert z['x'] in reg['missing_x'];existing[z['x']]=z
 if out.exists():assert len(existing)==len(reg['missing_x']);return
 write(H/'scalar_run_receipt.json',{'status':'running','pid':os.getpid(),'registration_sha256':sha(H/'scalar_registration.json'),'consumer_sha256':sha(__file__),'inherited_count':len(cat),'registered_new':len(reg['missing_x']),'resumed_successful':len(existing),'journal':str(journal)})
 t=time.monotonic()
 for x in reg['missing_x']:
  assert x not in cat
  if x not in existing:
   lo,hi=R.ratio_interval(F(x),192);z={'x':x,'lower':str(lo),'upper':str(hi),'origin':'path112_registered_new192'}
   with journal.open('a') as f:f.write(ser(z)+'\n');f.flush();os.fsync(f.fileno())
   existing[x]=z
  cat[x]={k:v for k,v in existing[x].items() if k!='x'}
 assert all(str(F(x)) in cat for x in reg['all_used_x']);write(out,cat)
 write(H/'scalar_run_receipt.json',{'status':'complete','pid':os.getpid(),'new_R':len(existing),'recomputed_old_R':0,'new_H':0,'seconds':time.monotonic()-t,'catalog_sha256':sha(out),'journal_sha256':sha(journal),'consumer_sha256':sha(__file__)})
 print('SCALARS_COMPLETE',len(existing),time.monotonic()-t,flush=True)
def factors():
 out=H/'factor_catalog.json'
 if out.exists():return
 reg=read(H/'scalar_registration.json');cat=read(H/'ratio_catalog.json');geos=read(H/'geometries.json.gz');plans=read(H/'plans.json.gz');factors={};t=time.monotonic()
 for row in plans:
  gid=row['geometry'];a=row['time'];p=P.darboux_plan(geos[gid]['offsets'],F(a),8);assert hashlib.sha256(ser(p).encode()).hexdigest()==row['plan_sha256'];e=Q.secant_envelope(p,cat);assert hashlib.sha256(ser(e['positive_endpoint_weights']).encode()).hexdigest()==row['alpha_beta_sha256'];f=e['factor_upper'];d=1<<192;rounded=F(ceil(f*d),d);assert rounded>=f;factors.setdefault(gid,{})[a]={'factor_upper':str(rounded),'exact_factor':str(f),'J_upper':str(e['J_upper']),'plan_sha256':row['plan_sha256']}
 write(out,factors);write(H/'factor_receipt.json',{'status':'complete','plans':len(plans),'seconds':time.monotonic()-t,'pid':os.getpid(),'all_plan_and_alpha_beta_hashes_match':True,'positive_weights_exact':True,'factor_outward_bits':192,'catalog_sha256':sha(out),'new_R':0,'new_H':0});print('FACTORS_COMPLETE',time.monotonic()-t,flush=True)
def node_value(row,fac,times):
 p=Path(row['packet_path']);assert sha(p)==row['packet_sha256'];ev=read(p);assert ev['point']==row['point'];ends=ev['endpoints'];assert [e['t'] for e in ends]==times
 mean=F(0)
 for i,e in enumerate(ends):
  nxt=ends[i+1]['H'] if i+1<len(ends) else {'lower':'0','upper':'0'};ml=F(e['H']['lower'])-F(nxt['upper']);mh=F(e['H']['upper'])-F(nxt['lower']);assert 0<ml<=mh;mean+=mh*fac[e['t']]
 mi=min(row['old_M_i'],row['old_U_i'],ceil(mean*N));li=row['old_L_i'];ui=row['old_U_i'];assert li<=mi<=ui
 l,u,m=F(li,N),F(ui,N),F(mi,N);lp=root(l**5,4);up=root(u**5,4);theta=F(mi-li,ui-li) if ui>li else F(0);chord=((1-theta)*lp+theta*up)/2**(4*row['level']);nu=min(row['old_selfU'],row['weight']*ceil(chord*N));assert row['old_selfL']<=nu
 return {'key':row['key'],'new_M_i':mi,'new_selfU':nu,'old_M_i':row['old_M_i'],'old_selfU':row['old_selfU'],'packet_sha256':row['packet_sha256']}
def propagate(st,r):
 ns=st['nodes'];k=r['key'];rt=k
 while ns[rt]['parent'] is not None:rt=ns[rt]['parent']
 old=ns[rt]['effU'];n=ns[k];assert r['new_selfU']<=n['selfU'];n['selfU']=r['new_selfU']
 while k is not None:
  n=ns[k];n['effU']=min(n['selfU'],sum(ns[c]['effU'] for c in n['children'])) if n['children'] else n['selfU'];assert n['effL']<=n['effU'];k=n['parent']
 return ns[rt]['effU']-old
def batch(kind):
 reg=read(H/'scalar_registration.json');rows=read(H/'nodes.json.gz');rm={r['key']:r for r in rows};times=reg['time_endpoints'];fac={g:{a:F(x['factor_upper']) for a,x in fs.items()} for g,fs in read(H/'factor_catalog.json').items()};sp=S/'state.json.gz'
 if sp.exists():st=read(sp)
 else:
  assert sha(reg['source103_state_path'])==reg['source103_state_sha256'];st=read(reg['source103_state_path']);st['path112']={'processed':[],'blocks':[]};write(sp,st)
 done=set(st['path112']['processed']);ns=st['nodes'];ml=max(ns[k]['level'] for k in st['frontier']);assert sum(ns[k]['weight']*2**(4*(ml-ns[k]['level'])) for k in st['frontier'])==2**(4*ml)
 if kind=='pilot':keys=reg['pilot_selected_physical_leaves'];assert keys==sorted(st['frontier'],key=lambda k:(-(ns[k]['selfU']-ns[k]['selfL']),ns[k]['level'],ns[k]['indices']))[:256] if not done else True;name='pilot'
 else:
  pilot=read(H/'pilot_receipt.json');assert pilot['status']=='complete' and pilot['strict_global_improvement'];keys=sorted((k for k in ns if k not in done),key=lambda k:(-(ns[k]['selfU']-ns[k]['selfL']),ns[k]['level'],ns[k]['indices']))[:4096];name='block%04d'%len(st['path112']['blocks'])
 if not keys:return
 rp=H/(name+'_registration.json');resultp=S/(name+'_results.jsonl');rcp=H/(name+'_receipt.json')
 if rp.exists():keys=read(rp)['keys']
 else:write(rp,{'keys':keys,'count':len(keys),'input_state_sha256':sha(sp),'source103_state_sha256':reg['source103_state_sha256'],'nodes_sha256':sha(H/'nodes.json.gz'),'factor_catalog_sha256':sha(H/'factor_catalog.json'),'consumer_sha256':sha(__file__),'range_L_U_unchanged':True,'new_H':0,'new_R':0,'pid':os.getpid(),'registration_unix':time.time()})
 before=sum(ns[k]['effU'] for k in st['roots']);write(rcp,{'status':'running','pid':os.getpid(),'count':len(keys),'before':budget(before),'registration_sha256':sha(rp)})
 saved={}
 if resultp.exists():saved={r['key']:r for r in map(json.loads,resultp.read_text().splitlines())}
 t=time.monotonic();up=before
 for ix,k in enumerate(keys):
  if k in done:continue
  r=saved.get(k)
  if r is None:
   r=node_value(rm[k],fac[rm[k]['geometry']],times)
   with resultp.open('a') as f:f.write(ser(r)+'\n');f.flush()
  delta=propagate(st,r);assert delta<=0;up+=delta;done.add(k);st['path112']['processed'].append(k)
  b=budget(up)
  if b['strict32'] and not (H/'first_strict32.json').exists():
   write(sp,st);snap=S/'first_strict32_state.json.gz';write(snap,st);write(H/'first_strict32.json',{'status':'pending_main_independent','budget':b,'state_path':str(snap),'state_sha256':sha(snap),'processed':len(done),'volume':'1'});print('FIRST_STRICT32',ser(b),flush=True);return
  if (ix+1)%256==0:print(name,ix+1,'cost',b['cost_display'],'seconds',time.monotonic()-t,flush=True)
 st['path112']['blocks'].append(name);write(sp,st);after=sum(ns[k]['effU'] for k in st['roots']);assert after==up and after<=before
 write(rcp,{'status':'complete','pid':os.getpid(),'count':len(keys),'processed_total':len(done),'before':budget(before),'after':budget(after),'strict_global_improvement':after<before,'state_sha256':sha(sp),'results_path':str(resultp),'results_sha256':sha(resultp),'seconds':time.monotonic()-t,'new_H':0,'new_R':0,'weighted_volume':'1','range_L_U_unchanged':True,'consumer_sha256':sha(__file__),'registration_sha256':sha(rp)})
 print('BATCH_COMPLETE',name,len(done),budget(after)['cost_display'],flush=True)
if __name__=='__main__':
 mode=sys.argv[1]
 if mode=='scalars':scalars()
 elif mode=='factors':factors()
 elif mode in ('pilot','block'):batch(mode)
