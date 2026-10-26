#!/usr/bin/env python3
"""One finite <=4096-cell invocation; checkpointed complete-space ledger."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod,lcm,ceil,floor
import argparse,gzip,hashlib,json,multiprocessing as mp,os,resource,signal,sys,time,traceback
import flint
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_global080');SCALE=1<<160;TABLE=None;TS=None;FP=None

def sha(p):
 v=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):v.update(b)
 return v.hexdigest()
def serial(o):return json.dumps(o,sort_keys=True,separators=(',',':'))
def atomic(p,o):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_text(serial(o)+'\n');os.replace(t,p)
def gzatomic(p,o):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_bytes(gzip.compress(serial(o).encode(),mtime=0));os.replace(t,p)
def rq(x):
 a=F(str(x));return fmpq(a.numerator,a.denominator)
def ep(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def quant(x,up):return ceil(x*SCALE) if up else floor(x*SCALE)
def alarm(sig,frame):raise TimeoutError('590 second program alarm')
def evidencepath(rank):return S/'evidence'/f'{rank//1024:02d}'/f'{rank:05d}.json.gz'
def checkpointpath(rank):return S/'partial'/f'{rank:05d}.json.gz'
def init_worker(table,ts,fp):
 global TABLE,TS,FP
 ctx.prec=384;TABLE=table;TS=ts;FP=fp
def qform(u):return sum(x*x for x in u)/4-sum(u)**2/36

def geometry(row):
 c=list(map(F,row['point']))[1:];lo=[F(i,32)-x for i,x in zip(row['indices'],c)];hi=[F(i+1,32)-x for i,x in zip(row['indices'],c)];assert all(l<0<h for l,h in zip(lo,hi));mu=[(l+h)/2 for l,h in zip(lo,hi)];var=[(h-l)**2/12 for l,h in zip(lo,hi)]
 return {'lo':list(map(str,lo)),'hi':list(map(str,hi)),'mu':list(map(str,mu)),'eta':str(qform(mu)+F(2,9)*sum(var)),'B2':str(max(qform(u) for u in product(*zip(lo,hi)))),'h_inf':str(max(map(abs,lo+hi))),'lambda':str(max(map(abs,mu))),'volume':'1/1048576'}
def bounds(rows,g):
 B2,h,eta,lam=map(F,(g['B2'],g['h_inf'],g['eta'],g['lambda']));L=arb(0);U=arb(0);M=arb(0);positive=0;taylor=0;expn=0;minmass=None
 for i,row in enumerate(rows):
  a=F(row['t']);cur=row['H'];nxt=rows[i+1]['H'] if i+1<len(rows) else {'lower':'0','upper':'0'};ml=F(cur['lower'])-F(nxt['upper']);mu=F(cur['upper'])-F(nxt['lower']);assert 0<ml<=mu,('mass_not_positive',i,str(ml),str(mu));minmass=ml if minmass is None else min(minmass,ml)
  ex=ep(arb(rq(h/a)).exp(),True);s=arb(rq(B2/(20*a*a))).sqrt();old=ep(((1+s*s).sqrt()+s)**10,True);R=min(ex,old);C=min(R,1+lam/a+F(21,8)*eta*R/(a*a));assert R>=1 and C>=1
  L+=arb(rq(ml/R));U+=arb(rq(mu*R));M+=arb(rq(mu*C));positive+=1;taylor+=C<R;expn+=ex<=old
 l=ep(L);rawu=ep(U,True);rawm=ep(M,True);u=min(F(3),rawu);m=min(u,rawm);assert 0<l<=m<=u
 # Quantize range/mean outward first; the chord then uses these reproducible endpoints.
 li=quant(l,False);ui=min(3*SCALE,quant(u,True));mi=min(ui,quant(m,True));assert 0<li<=mi<=ui
 l=F(li,SCALE);u=F(ui,SCALE);m=F(mi,SCALE);q=arb(5)/4;v=arb(rq(g['volume']));lp=arb(rq(l))**q;up=arb(rq(u))**q;chord=v*(lp+(up-lp)*arb(rq((m-l)/(u-l))));baseline=v*arb(rq(rows[0]['H']['lower']))**q
 return {'L_i':li,'U_i':ui,'M_i':mi,'integral_lower_i':quant(ep(v*lp),False),'integral_upper_i':quant(ep(chord,True),True),'relative_baseline_upper_i':quant(ep(chord/baseline,True),True),'positive_bands':positive,'taylor_bands':taylor,'exp_bands':expn,'global_cap_used':rawu>3,'minimum_band_mass_lower':str(minmass),'raw_H_upper':str(rawu),'raw_mean_upper':str(rawm)}

def worker(row):
 began=time.monotonic();rank=row['rank'];out=evidencepath(rank);cp=checkpointpath(rank)
 if out.exists():
  r=json.loads(gzip.decompress(out.read_bytes()));assert r['fingerprint']==FP and r['point']==row['point'] and r['rootid']==row['rootid'];return {'rank':rank,'rootid':row['rootid'],'status':'reused_completed','weight':row['weight'],'bounds':r['bounds'],'evidence_sha256':sha(out),'runtime':r['runtime']}
 state={'fingerprint':FP,'rank':rank,'rootid':row['rootid'],'point':row['point'],'next_cone':0,'weights':{},'endpoints':[]};runtime={'pid':os.getpid(),'began_unix':time.time(),'new_H_count':0,'reused_H_count':0,'new_jet_count':0}
 try:
  if cp.exists():state=json.loads(gzip.decompress(cp.read_bytes()));assert state['fingerprint']==FP and state['point']==row['point']
  c=list(map(F,row['point']));D=lcm(*(x.denominator for x in c));nodes=[sum(int(x*D)*z for x,z in zip(c,v)) for v in TABLE['vertices']];assert len(set(nodes))==280;weights={int(k):rq(v) for k,v in state['weights'].items()};stage=time.monotonic();start_cone=state['next_cone']
  for ci in range(start_cone,3504):
   cone=TABLE['simplices'][ci];ks=[nodes[k] for k in cone['vertices']];det=int(cone['det_numerator'])
   for j,k in enumerate(ks):weights[k*k]=weights.get(k*k,fmpq(0))-fmpq(det,2*prod(k-l for n,l in enumerate(ks) if n!=j))
   state['next_cone']=ci+1
   if (ci+1)%512==0 or ci==3503:
    state['weights']={str(k):str(v) for k,v in weights.items() if v};gzatomic(cp,state)
  assert len(weights)==140 and sum(weights.values(),fmpq(0))==0;runtime['stencil_seconds']=time.monotonic()-stage;stage=time.monotonic();factor=arb(D)**4/(144*(2*arb.pi())**5);oldn=len(state['endpoints'])
  for i in range(oldn,len(TS)):
   ts=TS[i];t=F(ts)
   if t==F(1,64):ball=row['cached_H1over64'];origin='cached066_H1over64';runtime['reused_H_count']+=1
   else:
    val=factor*sum((arb(w)*arb(rq((144*D*t)**2+k2)).log() for k2,w in sorted(weights.items())),arb(0));assert val.is_finite() and val>0,('invalid_H',ts,str(val));ball={'lower':str(ep(val)),'upper':str(ep(val,True))};origin='new080_scalar';runtime['new_H_count']+=1
   state['endpoints'].append({'t':ts,'H':ball,'origin':origin})
   if (i+1)%8==0 or i+1==len(TS):
    gzatomic(cp,state)
  runtime['endpoint_seconds']=time.monotonic()-stage;stage=time.monotonic();g=geometry(row);b=bounds(state['endpoints'],g);runtime.update({'bounds_seconds':time.monotonic()-stage,'seconds':time.monotonic()-began,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'new_cone_count':3504-start_cone})
  r={'fingerprint':FP,'rank':rank,'rootid':row['rootid'],'indices':row['indices'],'point':row['point'],'weight':row['weight'],'D':D,'geometry':g,'weights':{str(k):str(v) for k,v in sorted(weights.items())},'endpoints':state['endpoints'],'bounds':b,'runtime':runtime};gzatomic(out,r);cp.unlink();return {'rank':rank,'rootid':row['rootid'],'status':'complete','weight':row['weight'],'bounds':b,'evidence_sha256':sha(out),'runtime':runtime}
 except BaseException as e:
  state['diagnostic']={'error':repr(e),'traceback':traceback.format_exc(),'runtime':runtime,'elapsed_seconds':time.monotonic()-began};gzatomic(cp,state);return {'rank':rank,'rootid':row['rootid'],'status':'failed','error':repr(e),'pid':os.getpid()}

def space_ledger(state,cap):
 num=state['sum_upper_weighted_units']+(1048576-state['computed_weight'])*cap;UH=F(num,SCALE);E=F(json.loads((H/'cap_scalar.json').read_text())['E076_norm']['upper']);I=(arb(rq(UH))**(arb(4)/5)+arb(rq(E)))**(arb(5)/4)
 goal=arb(rq('1/25000'));budget=(goal**(arb(4)/5)-arb(rq(E)))**(arb(5)/4)
 return {'completed_cells':state['next_rank'],'computed_weight':state['computed_weight'],'computed_weighted_volume':str(F(state['computed_weight'],1048576)),'uncomputed_weighted_volume':str(F(1048576-state['computed_weight'],1048576)),'computed_UH_upper':str(F(state['sum_upper_weighted_units'],SCALE)),'uncomputed_UH_upper':str(F((1048576-state['computed_weight'])*cap,SCALE)),'full_face_UH_upper':str(UH),'full_face_UH_upper_display':float(UH),'full_I_upper':str(ep(I,True)),'full_I_upper_display':float(I),'complete_E076_norm_upper':str(E),'UH_sufficient_budget_for_I_le_1over25000_lower':str(ep(budget)),'UH_sufficient_budget_display':float(budget),'no_tube_added':True,'new_global_g':False,'new_Cn_claim':False}

def main():
 began=time.monotonic();ctx.prec=384;ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');ap.add_argument('--max-cells',type=int,default=4096);args=ap.parse_args();assert 1<=args.max_cells<=4096
 reg=json.loads((H/'registration.json').read_text());policy=json.loads((H/'execution_policy_v2.json').read_text());assert sha(__file__)==policy['new_engine_sha256'];assert sha(H/'engine.py')==policy['old_engine_sha256'];assert sha(H/'roots.json.gz')==reg['roots_sha256'];assert sha(H/'geometry017.json.gz')==reg['table_sha256'];assert sha(H/'cap_scalar.json')==reg['cap_scalar_sha256'];assert sha(H/'mechanism.md')==reg['mechanism_sha256']
 fp=hashlib.sha256((sha(__file__)+sha(H/'registration.json')+sha(H/'execution_policy_v2.json')).encode()).hexdigest();roots=json.loads(gzip.decompress((H/'roots.json.gz').read_bytes()));table=json.loads(gzip.decompress((H/'geometry017.json.gz').read_bytes()));cap=quant(ep(arb(3)**(arb(5)/4)/arb(1048576),True),True)
 stpath=S/'state.json';state=json.loads(stpath.read_text()) if stpath.exists() else {'fingerprint':fp,'next_rank':0,'sum_lower_weighted_units':0,'sum_upper_weighted_units':0,'computed_weight':0,'blocks':[]};assert state['fingerprint'] in (fp,policy['old_fingerprint']);assert state['next_rank']<52360,'initial layer complete; no adaptive continuation'
 if state['fingerprint']==policy['old_fingerprint']:
  assert state['next_rank']==policy['transition_rank'] and sha(stpath)==policy['parent_state_sha256'],'v2 bridge must start at frozen exact parent state'
 blockid=len(state['blocks']);start=state['next_rank'];end=min(52360,start+args.max_cells);selected=roots[start:end];bid=f'{blockid:04d}';bp=S/'blocks'/bid;bp.mkdir(parents=True,exist_ok=True);bregpath=H/('block'+bid+'_registration.json');bcp=bp/'checkpoint.json';selection_sha=hashlib.sha256(serial(selected).encode()).hexdigest()
 if bregpath.exists():
  assert args.resume,'existing block registration requires --resume';breg=json.loads(bregpath.read_text());assert breg['start']==start and breg['end']==end and breg['selection_sha256']==selection_sha and breg['fingerprint']==fp
 else:
  assert not args.resume,'no prior block to resume';breg={'id':'upper_global080.block'+bid,'fingerprint':fp,'engine_sha256':sha(__file__),'registration_sha256':sha(H/'registration.json'),'parent_state_sha256':sha(stpath) if stpath.exists() else None,'parent_block_receipt_sha256':state['blocks'][-1]['receipt_sha256'] if state['blocks'] else None,'start':start,'end':end,'count':end-start,'selection_sha256':selection_sha,'rootids':[r['rootid'] for r in selected],'worker_count':12,'execution_policy_sha256':sha(H/'execution_policy_v2.json'),'wall_deadline':None,'workload_checkpoint':'512 cones and8 endpoints','registered_unix':time.time()};atomic(bregpath,breg)
 receipt={'blockid':blockid,'fingerprint':fp,'registration_sha256':sha(bregpath),'pid':os.getpid(),'began_unix':time.time(),'status':'running','execution_policy_sha256':sha(H/'execution_policy_v2.json'),'engine_sha256':sha(__file__),'start':start,'end':end,'jobs':[],'python':sys.executable,'flint':flint.__version__};atomic(bp/'receipt.json',receipt);pool=mp.get_context('spawn').Pool(12,initializer=init_worker,initargs=(table,reg['time_endpoints'],fp));allrows=[];failure=None
 try:
  for r in pool.imap_unordered(worker,selected,chunksize=1):
   allrows.append(r)
   if r['status']=='failed':failure=r;pool.terminate();break
   if len(allrows)%64==0:atomic(bcp,{'fingerprint':fp,'block_registration_sha256':sha(bregpath),'rows':allrows});print('block',bid,'completed',len(allrows),'wall',round(time.monotonic()-began,3),flush=True)
  if failure is None:pool.close()
  pool.join();receipt['workers_joined']=True
 except BaseException as e:
  pool.terminate();pool.join();receipt.update({'workers_joined':True,'exception':repr(e)})
 finally:atomic(bcp,{'fingerprint':fp,'block_registration_sha256':sha(bregpath),'rows':allrows})
 complete=len(allrows)==len(selected) and all(r['status'] in ('complete','reused_completed') for r in allrows)
 receipt.update({'status':'complete' if complete else 'failed' if failure else 'checkpointed','failure':failure,'wall_seconds':time.monotonic()-began,'completed_count':sum(r['status'] in ('complete','reused_completed') for r in allrows),'worker_pids':sorted({r.get('runtime',{}).get('pid',r.get('pid')) for r in allrows if r.get('runtime',{}).get('pid',r.get('pid')) is not None}),'peak_parent_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'children_user_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime,'children_system_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime})
 if complete:
  allrows.sort(key=lambda r:r['rank']);assert [r['rank'] for r in allrows]==list(range(start,end));gzatomic(bp/'ledger.json.gz',allrows)
  for r in allrows:state['sum_lower_weighted_units']+=r['weight']*r['bounds']['integral_lower_i'];state['sum_upper_weighted_units']+=r['weight']*r['bounds']['integral_upper_i'];state['computed_weight']+=r['weight']
  state['next_rank']=end;state['fingerprint']=fp;receipt['ledger_sha256']=sha(bp/'ledger.json.gz');receipt['space_ledger']=space_ledger(state,cap);receipt['stage_seconds_sum']={k:sum(r['runtime'][k] for r in allrows) for k in ('stencil_seconds','endpoint_seconds','bounds_seconds','seconds')};receipt['new_H_count']=sum(r['runtime']['new_H_count'] for r in allrows);receipt['reused_H_count']=sum(r['runtime']['reused_H_count'] for r in allrows)
  atomic(bp/'receipt.json',receipt);atomic(H/('block'+bid+'_receipt.json'),receipt);state['blocks'].append({'id':blockid,'start':start,'end':end,'registration_sha256':sha(bregpath),'receipt_sha256':sha(bp/'receipt.json'),'ledger_sha256':sha(bp/'ledger.json.gz')});atomic(stpath,state);atomic(H/'progress.json',{'state_sha256':sha(stpath),'fingerprint':fp,'blocks_completed':len(state['blocks']),'initial_layer_complete':end==52360,'space_ledger':receipt['space_ledger']})
 else:atomic(bp/'receipt.json',receipt);atomic(H/('block'+bid+'_receipt.json'),receipt)
 print(serial({k:v for k,v in receipt.items() if k!='jobs'}),flush=True)
if __name__=='__main__':main()
