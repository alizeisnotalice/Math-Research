#!/usr/bin/env python3
"""Second audit178: independent full accepted-H reader. No old math imports.
prepare hashes finite sources only. run requires reviewed plan SHA and new output.
H enclosures remain shared input; no DD/log/H kernel is evaluated.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from itertools import combinations_with_replacement,product
from math import factorial,prod,isqrt
import argparse,gzip,hashlib,json,os,multiprocessing as mp,time,traceback,struct,sys
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parent
PACK=PROJECT/'upper_integral167/packet'
N=1<<160
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else Path(p).open()) as f:return json.load(f)
def stream(p):
 with gzip.open(p,'rt') as f:
  for line in f:yield json.loads(line)
def save(p,x):
 p=Path(p);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n');os.replace(t,p)
def resolve(n):
 p=PACK/n
 if p.exists():return p
 refs=read(PACK/'repository_reuse.json')['files'];return (PACK/refs[n]['relative_from_repository_packet']).resolve()
def prepare():
 assert not (HERE/'plan.json').exists()
 done=read(PACK/'completion.json');files={str(resolve(n)):{'sha256':h,'bytes':resolve(n).stat().st_size} for n,h in done['files'].items()}
 for n in ('completion.json','repository_reuse.json','stage_packet.py','verify_packet.py'):
  p=PACK/n;files[str(p)]={'sha256':digest(p),'bytes':p.stat().st_size}
 for p in (HERE/'reader178.py',PROJECT/'primary/finite_source_transfer.md',PROJECT/'round018/primary_integral156_acceptance.json'):
  files[str(p)]={'sha256':digest(p),'bytes':p.stat().st_size}
 for p,r in files.items():assert digest(p)==r['sha256']
 save(HERE/'input_lock.json',{'schema':'audit178-source-lock-v1','files':files})
 plan={'schema':'audit178-reader-plan-v1','status':'prepared_not_run','code_sha256':digest(HERE/'reader178.py'),'input_lock_sha256':digest(HERE/'input_lock.json'),'packet':str(PACK),'exact_shards':15,'exact_nodes':97971,'exact_bands':7151883,'time_endpoints':73,'roots':52360,'complete_tree_nodes':97972,'checkpoint_nodes':256,'output':str(HERE/'run001'),'workers_respawn':False,'wallcap':None,'old_H_interval_truth':'shared accepted enclosure dependency, not recomputed','new_math_scope':['independent binary parse and canonical projection','all scalar R upper certificates','independent piecewise-polynomial weight geometry, secant J,112 and150/156 factors','all real offsets and73 positive masses','range and mean from exact integer H arithmetic','Newton fourth-root power brackets and both old112/new156 chords','156 total gain once from112, full canonical root-child coverage and min tree','139 saved fourth-power brackets and128-panel sum, exact complete I156 and final-root interval comparisons'],'forbidden':['H/log/DD/gradient recalculation','sampling','lower-bound advancement','new grid or time point','rerun prior successful script','automatic resume']}
 save(HERE/'plan.json',plan)
 print(json.dumps({'status':'prepared_not_run','plan_sha256':digest(HERE/'plan.json'),'files':len(files),'bytes':sum(r['bytes'] for r in files.values()),'new_math_calls':0}))
def locks(plan):
 assert digest(__file__)==plan['code_sha256'] and digest(HERE/'input_lock.json')==plan['input_lock_sha256']
 for p,r in read(HERE/'input_lock.json')['files'].items():assert Path(p).stat().st_size==r['bytes'] and digest(p)==r['sha256'],p
def ceildiv(a,b):return -(-a//b)
def dyad(q,bits):
 q=Q(q);v=q*(1<<bits);assert v.denominator==1;return v.numerator
def root4floor(n):
 assert n>=0
 if n==0:return 0
 x=1<<ceildiv(n.bit_length(),4)
 while True:
  y=(3*x+n//(x*x*x))//4
  if y>=x:break
  x=y
 while x**4>n:x-=1
 while (x+1)**4<=n:x+=1
 assert x**4<=n<(x+1)**4
 return x
def power224(i):
 a=i**5<<96;r=root4floor(a)
 return r if r**4==a else r+1
def chord(row,mi):
 l,u=row['L_i'],row['U_i'];assert 0<l<=mi<=u
 a,b=power224(l),power224(u)
 if l==u:return row['weight']*ceildiv(a,1<<(64+4*row['level']))
 return row['weight']*ceildiv((u-mi)*a+(mi-l)*b,(u-l)<<(64+4*row['level']))
def polymul(a,b):
 c=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def inte(a,l,r,extra=0):return sum((x*(r**(j+extra+1)-l**(j+extra+1))/Q(j+extra+1) for j,x in enumerate(a)),Q(0))
def geometry(offsets):
 ab=[tuple(map(Q,x)) for x in offsets];knots=sorted({Q(0)}|{abs(t) for x in ab for t in x});pieces=[]
 for l,r in zip(knots,knots[1:]):
  mid=(l+r)/2;cdf=[Q(1)]
  for a,b in ab:
   assert a<b
   upper=[b,Q(0)] if b<mid else [Q(0),Q(1)]
   lower=[a,Q(0)] if a>-mid else [Q(0),Q(-1)]
   v=[(upper[j]-lower[j])/(b-a) for j in (0,1)]
   if v[0]+v[1]*mid<=0:cdf=[Q(0)];break
   cdf=polymul(cdf,v)
  surv=[-x for x in cdf];surv[0]+=1;pieces.append([l,r,surv,None])
 tail=Q(0)
 for p in reversed(pieces):
  l,r,s,_=p;w=[tail+sum((x*r**(j+1)/Q(j+1) for j,x in enumerate(s)),Q(0))]+[-x/Q(j+1) for j,x in enumerate(s)]
  p[3]=w;tail+=inte(s,l,r)
 j0=sum((inte(p[3],p[0],p[1]) for p in pieces),Q(0));assert j0==sum((inte(p[2],p[0],p[1],1) for p in pieces),Q(0))
 return ab,knots[-1],max(abs((a+b)/2) for a,b in ab),j0,pieces
def prove_R(cat):
 for xs,row in cat.items():
  x=Q(xs);assert x>=0;rs=[]
  for den in (16,256):
   q=1+x*x/den;D=1<<512;n=isqrt(q.numerator*D*D//q.denominator);u=Q(n+(n*n*q.denominator<q.numerator*D*D),D);assert u*u>=q;rs.append(u)
  assert (rs[0]+x/4)**2*(rs[1]+x/16)**8<=Q(row['upper'])
def factors(times):
 cat=read(resolve('ratios112.json'));old=read(resolve('factors112.json'));geos=read(resolve('geometries.json.gz'));new=read(resolve('mean_factors.json'))['geometry_factors'];prove_R(cat)
 result={};shapes={}
 for g,o in geos.items():
  ab,h,mu,j0,pieces=geometry(o['offsets']);shapes[g]=ab
  assert (h,mu,j0,2*j0)==tuple(Q(o['profile'][k]) for k in ('h','mu_infinity','J0','M2'))
  def R(z):return Q(1) if z==0 else Q(cat[str(z)]['upper'])
  js={}
  for t in times:
   z=1/Q(t);sec=Q(0);darb=Q(0)
   for l,r,s,w in pieces:
    for i in range(8):
     a=l+(r-l)*Q(i,8);b=l+(r-l)*Q(i+1,8);mass=inte(w,a,b);moment=inte(w,a,b,1);assert mass>=0
     left=(b*mass-moment)/(b-a);right=(moment-a*mass)/(b-a);assert left>=0 and right>=0 and left+right==mass
     ra,rb=R(a*z),R(b*z);sec+=left*min(ra,rb)+right*rb;darb+=mass*rb
   ju=min(sec,darb,j0*R(h*z));assert ju==Q(old[g][t]['J_upper']);js[t]=ju
   oldraw=min(R(h*z),1+mu*z+Q(25,24)*z*z*ju)
   assert Q(old[g][t]['factor_upper'])==Q(ceildiv(oldraw.numerator<<192,oldraw.denominator),1<<192)
  cs=[]
  for j,t in enumerate(times):
   a=1/Q(t);b=0 if j==72 else 1/Q(times[j+1]);theta=(Q(4,5)*(a**5-b**5)/(a**4-b**4)-b)/(a-b);assert 0<theta<1
   if j==72:assert theta==Q(4,5)
   ra,rb=R(h*a),R(h*b);ta=1+mu*a+Q(25,24)*a*a*js[t];tb=Q(1) if j==72 else 1+mu*b+Q(25,24)*b*b*js[times[j+1]]
   raw=min(Q(old[g][t]['factor_upper']),theta*ra+(1-theta)*min(ra,rb),theta*ta+(1-theta)*min(ta,tb))
   c=ceildiv(raw.numerator<<192,raw.denominator);assert Q(c,1<<192)==Q(new[g][j]['new_upper']) and Q(c,1<<192)<=Q(old[g][t]['factor_upper']);cs.append(c)
  result[g]=cs
 cats={g:read(resolve('range_'+g+'.json')) for g in ('092','095','103')}
 for cat in cats.values():prove_R(cat)
 return result,shapes,{g:{s:dyad(r['upper'],192) for s,r in cat.items()} for g,cat in cats.items()}
def take(f,n):
 b=f.read(n);assert len(b)==n;return b
def uint(f):
 n=int.from_bytes(take(f,2),'big');assert 0<n<65536;b=take(f,n);assert n==1 or b[0]!=0;return int.from_bytes(b,'big')
def worker(sid,rows,times,fac,shapes,cats,out):
 began=time.monotonic();dest=Path(out)/('shard%02d'%sid);dest.mkdir();count=0;bands=0;selfs={};chain='0'*64;chunk={};prefix_chunks=[]
 ident={'code_sha256':digest(__file__),'input_lock_sha256':digest(HERE/'input_lock.json'),'plan_sha256':digest(HERE/'plan.json')}
 try:
  p=resolve('endpoints%02d.bin.gz'%sid);before=digest(p)
  with gzip.open(p,'rb') as f:
   assert take(f,8)==b'I167H001' and int.from_bytes(take(f,4),'big')==sid and int.from_bytes(take(f,4),'big')==len(rows)
   for row in rows:
    assert row['id']%15==sid and row['shard']==sid and int.from_bytes(take(f,4),'big')==row['id'];vals=[]
    for t in times:
     e=int.from_bytes(take(f,4),'big',signed=True);assert 0<=e<=408;l,w=uint(f),uint(f);assert l>0;vals.append((l,l+w,e))
    proj,raw,ancestry=(take(f,32).hex() for _ in range(3));assert raw==row['packet_sha256']
    ends=[{'t':t,'H':{'lower':str(Q(l,1<<e)),'upper':str(Q(u,1<<e))}} for t,(l,u,e) in zip(times,vals)]
    projected={'point':[str(Q(s)) for s in row['point']],'endpoints':ends};assert hashlib.sha256(json.dumps(projected,sort_keys=True,separators=(',',':')).encode()).hexdigest()==proj
    assert row['_projection']=={'id':row['id'],'source_gzip_sha256':raw,'original_raw_sha256':ancestry,'normalized_projection_sha256':proj}
    lev=row['level'];idx=row['indices'];assert row['key']==str(lev)+':'+','.join(map(str,idx));assert idx==sorted(idx,reverse=True) and all(0<=i<(1<<lev) for i in idx)
    assert Q(row['point'][0])==1 and row['weight']==factorial(4)//prod(factorial(n) for n in Counter(idx).values())
    ab=[(Q(i,1<<lev)-Q(c),Q(i+1,1<<lev)-Q(c)) for i,c in zip(idx,row['point'][1:])];assert ab==shapes[row['geometry']] and max(abs(x) for pair in ab for x in pair)==Q(row['h_infinity'])
    E=max(v[2] for v in vals);L=[l<<(E-e) for l,u,e in vals]+[0];U=[u<<(E-e) for l,u,e in vals];upper=mean=lower=0
    for j,t in enumerate(times):
     ml,mh=L[j]-(U[j+1] if j<72 else 0),U[j]-L[j+1];assert 0<ml<=mh
     R=cats[row['range_source']][str(Q(row['h_infinity'])/Q(t))];assert R>=1<<192
     lower+=(ml<<(704-E))//R;upper+=mh*R;mean+=mh*fac[row['geometry']][j];bands+=1
    assert row['L_i']<<(512-160)<=lower
    assert row['U_i']<=3*N
    assert row['U_i']==3*N or upper<=row['U_i']<<(E+32)
    mi=ceildiv(mean,1<<(E+32));assert mi==row['new_mean_i']<row['m112_i'];assert row['L_i']<row['U_i'] and row['U_i']>=1<<144
    assert row['L_i']<=mi<=row['m112_i']<=row['U_i']
    assert chord(row,row['m112_i'])<=row['self112U'] and chord(row,mi)<=row['new_selfU']
    gain=row['weight']*(row['m112_i']-mi)//(1<<(4*lev+4));assert row['new_selfU']==row['self112U']-gain
    selfs[row['key']]=row['new_selfU'];chunk[row['key']]=row['new_selfU'];count+=1;chain=hashlib.sha256((chain+row['key']+str(mi)+str(row['new_selfU'])).encode()).hexdigest()
    if count%256==0:
     cp=dest/('prefix%05d.json'%count);save(cp,{'identity':ident,'verified_selfs':chunk});prefix_chunks.append({'name':cp.name,'sha256':digest(cp),'nodes':len(chunk)});chunk={}
     save(dest/'checkpoint.json',{'pid':os.getpid(),'nodes':count,'bands':bands,'chain':chain,'identity':ident,'verified_prefix_chunks':prefix_chunks})
   assert f.read(1)==b''
  assert before==digest(p)
  save(dest/'selfs.json',selfs);save(dest/'result.json',{'status':'complete','pid':os.getpid(),'sid':sid,'nodes':count,'bands':bands,'seconds':time.monotonic()-began,'chain':chain,'selfs_sha256':digest(dest/'selfs.json'),'new_H_calls':0})
 except BaseException:
  save(dest/'failure.json',{'pid':os.getpid(),'sid':sid,'nodes':count,'bands':bands,'traceback':traceback.format_exc(),'identity':ident,'verified_prefix_chunks':prefix_chunks});raise
def check139():
 inp=read(resolve('139/actual/inputs.json'));cert=read(resolve('139/actual/certificate.json'));res=read(resolve('139/actual/results.json'))
 e,c,m=map(Q,[inp['epsilon'],inp['Cbar_upper'],inp['M_upper']]);pl,pu=map(Q,inp['pi_interval']);assert e==Q(1,512) and c>0 and m==Q(read(resolve('certificate.json'))['M132_cap'])<3 and 0<pl<=pu
 d=c/(2*e*e);k=m-d;end=[];means=[];geo=[]
 for j in range(128):
  mass=Q((j+1)**3-j**3,3*128**3);moment=Q(3*((j+1)**7-j**7),7*128**4*((j+1)**3-j**3));theta=(moment*128**4-j**4)/((j+1)**4-j**4);assert mass>0 and 0<=theta<=1;geo.append((mass,moment,theta))
 assert len(cert['roots'])==257
 for n,r in enumerate(cert['roots']):
  idx=n if n<=128 else n-129;kind='endpoint' if n<=128 else 'mean';arg=d+k*(Q(idx**4,128**4) if kind=='endpoint' else geo[idx][1]);l,u=Q(r['lower']),Q(r['upper']);assert r['kind']==kind and r['index']==idx and Q(r['argument'])==arg and 0<=l<=u and l**4<=arg<=u**4 and u-l<=Q(1,1<<192)
  (end if kind=='endpoint' else means).append((l,u))
 lo=up=Q(0);assert len(cert['panels'])==128
 for j,r in enumerate(cert['panels']):
  mass,moment,theta=geo[j];a=mass*((1-theta)*end[j][0]+theta*end[j+1][0]);b=mass*means[j][1];assert 0<=a<=b;assert r['index']==j and list(map(Q,[r['mass'],r['mean_z4'],r['theta'],r['lower'],r['upper']]))==[mass,moment,theta,a,b];lo+=a;up+=b
 jl=5*e*e/(2*pu)*lo;ju=5*e*e/(2*pl)*up
 assert [jl,ju]==list(map(Q,res['auxiliary_J_interval'])) and ju==Q(read(resolve('certificate.json'))['qualified_J139_upper'])
 return ju

def treecheck(selfs):
 tree=read(resolve('tree.json.gz'));ns=tree['nodes'];roots=tree['roots'];special=read(resolve('certificate.json'))['unused_legacy_self_key'];assert len(roots)==len(set(roots))==52360 and len(ns)==97972 and set(ns)-set(selfs)=={special}
 expected={tuple(reversed(v)) for v in combinations_with_replacement(range(32),4)};assert {tuple(ns[k]['indices']) for k in roots}==expected
 incoming=Counter();eff={};leafcount=innercount=0
 for k,n in ns.items():
  i=n['indices'];l=n['level'];assert k==str(l)+':'+','.join(map(str,i)) and i==sorted(i,reverse=True) and all(0<=x<(1<<l) for x in i)
  assert n['weight']==factorial(4)//prod(factorial(x) for x in Counter(i).values())
  assert (n['parent'] is None)==(k in roots)
  if k in roots:assert l==5
  ch=n['children'] or [];assert len(ch)==len(set(ch))
  for child in ch:incoming[child]+=1;assert ns[child]['parent']==k and ns[child]['level']==l+1
  if ch:
   innercount+=1;children={tuple(sorted((2*x+b for x,b in zip(i,bits)),reverse=True)) for bits in product((0,1),repeat=4)}
   assert children=={tuple(ns[c]['indices']) for c in ch} and len(ch)==len(children) and sum(ns[c]['weight'] for c in ch)==16*n['weight']
  else:leafcount+=1;assert k!=special
 assert innercount==3392 and leafcount==94580
 assert all(incoming[k]==(0 if k in roots else 1) for k in ns)
 visited=set();stack=list(roots)
 while stack:
  k=stack.pop();assert k not in visited;visited.add(k);stack.extend(ns[k]['children'] or [])
 assert visited==set(ns)
 for k in sorted(ns,key=lambda k:-ns[k]['level']):
  child=ns[k]['children'] or []
  eff[k]=min(selfs[k],sum(eff[c] for c in child)) if child and k!=special else sum(eff[c] for c in child) if child else selfs[k]
 uh=Q(sum(eff[k] for k in roots),N);assert uh==Q(read(resolve('certificate.json'))['candidate_UH'])
 return uh

def run(args):
 assert digest(HERE/'plan.json')==args.reviewed_plan_sha
 plan=read(HERE/'plan.json');assert plan['status']=='prepared_not_run';locks(plan)
 out=Path(args.output).resolve();assert str(out)==plan['output'];out.mkdir(exist_ok=False);began=time.monotonic()
 save(out/'registration.json',{'status':'running','master_pid':os.getpid(),'plan_sha256':args.reviewed_plan_sha,'start_unix':time.time(),'H_enclosure_truth_shared':True})
 info=read(resolve('certificate.json'));times=info['times'];assert len(times)==73 and list(map(Q,times))==sorted(set(map(Q,times))) and Q(times[0])==Q(1,512) and Q(times[-1])==1
 pm={r['id']:r for r in stream(resolve('projection_map.jsonl.gz'))};assert set(pm)==set(range(97971))
 groups=[[] for _ in range(15)];ids=set();keys=set()
 for row in stream(resolve('node_index.jsonl.gz')):
  assert row['id'] not in ids and row['key'] not in keys;ids.add(row['id']);keys.add(row['key']);row['_projection']=pm[row['id']];groups[row['shard']].append(row)
 assert ids==set(range(97971)) and [len(g) for g in groups]==[6532]*6+[6531]*9
 fac,shapes,cats=factors(times);save(out/'factors_result.json',{'status':'complete','geometry_factors':365,'geometry_profiles':5,'all_R_catalog_uppers_requalified':True,'new_H_log_calls':0})
 ctx=mp.get_context('spawn');ps=[ctx.Process(target=worker,args=(sid,groups[sid],times,fac,shapes,cats,str(out))) for sid in range(15)]
 for p in ps:p.start()
 save(out/'workers.json',{'master_pid':os.getpid(),'worker_pids':[p.pid for p in ps],'fixed_shards':15})
 fault=False
 while any(p.is_alive() for p in ps):
  if any(p.exitcode not in (None,0) for p in ps):
   fault=True
   for p in ps:
    if p.is_alive():p.terminate()
   break
  for p in ps:p.join(.1)
 for p in ps:p.join()
 assert not fault and all(p.exitcode==0 for p in ps),'worker failure; see preserved shard failures'
 allself={};nodes=bands=0
 for sid in range(15):
  q=out/('shard%02d'%sid);r=read(q/'result.json');assert r['status']=='complete' and digest(q/'selfs.json')==r['selfs_sha256'];z=read(q/'selfs.json');assert not set(z)&set(allself);allself.update(z);nodes+=r['nodes'];bands+=r['bands']
 assert nodes==97971 and bands==7151883
 uh=treecheck(allself);ju=check139();I=uh+ju;accepted=read(PROJECT/'round018/primary_integral156_acceptance.json');g=Q(accepted['qualified_g']);power=3125*I**4/g**5
 assert uh==Q(accepted['new_UH']) and I==Q(accepted['complete_I_upper']) and power==Q(accepted['best_upper']['fourth_power'])
 a,b=map(Q,accepted['best_upper']['interval']);assert 0<a<=b and a**4<=power<=b**4 and power<Q(accepted['simple_rational_upper'])**4
 locks(plan)
 save(out/'results.json',{'status':'complete_full_arithmetic_from_accepted_H','pid':os.getpid(),'seconds':time.monotonic()-began,'worker_pids':[p.pid for p in ps],'worker_exitcodes':[p.exitcode for p in ps],'nodes':nodes,'bands':bands,'roots':52360,'tree_nodes':97972,'UH156':str(uh),'J139_upper':str(ju),'I156_upper':str(I),'final_U156_root_interval_verified':True,'H_interval_truth_shared_accepted_dependency':True,'new_H_DD_log_gradient_calls':0,'previous_successful_script_replays':0,'all_source_locks_before_and_after':True,'lower_changed':False})
 print(json.dumps({'status':'complete','nodes':nodes,'bands':bands,'seconds':time.monotonic()-began,'results_sha256':digest(out/'results.json')}))
if __name__=='__main__':
 a=argparse.ArgumentParser();s=a.add_subparsers(dest='mode',required=True);s.add_parser('prepare');r=s.add_parser('run');r.add_argument('--reviewed-plan-sha',required=True);r.add_argument('--output',required=True);args=a.parse_args()
 if args.mode=='prepare':prepare()
 else:run(args)
