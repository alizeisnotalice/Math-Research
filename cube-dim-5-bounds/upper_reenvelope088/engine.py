#!/usr/bin/env python3
"""One <=4096-node pure mean/chord reenvelope invocation; zero H kernels."""
from pathlib import Path
from fractions import Fraction as F
from math import ceil
import argparse,gzip,json,hashlib,multiprocessing as mp,os,time,resource,sys
from flint import arb,ctx,fmpq
import moment085 as M
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_reenvelope088');N=1<<160;FP=None;TS=None

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def atomic(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_text(serial(x)+'\n');os.replace(t,p)
def gz(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_bytes(gzip.compress(serial(x).encode(),mtime=0));os.replace(t,p)
def rq(x):
 x=F(str(x));return fmpq(x.numerator,x.denominator)
def ep(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def init(fp,ts):
 global FP,TS
 ctx.prec=384;FP=fp;TS=ts
def outpath(rank):return S/'nodes'/f'{rank//1024:02d}'/f'{rank:05d}.json.gz'
def worker(row):
 began=time.monotonic();out=outpath(row['rank'])
 if out.exists():
  r=json.loads(gzip.decompress(out.read_bytes()));assert r['fingerprint']==FP and r['key']==row['key'] and r['source_evidence_sha256']==row['source_evidence_sha256'];return r
 try:
  path=Path(row['source_evidence_path']);assert sha(path)==row['source_evidence_sha256'];old=json.loads(gzip.decompress(path.read_bytes()));assert old['indices']==row['indices'];point=list(map(F,old['point'] if 'point' in old else old['anchor']));level=row['level'];bounds=[(F(i,2**level),F(i+1,2**level)) for i in row['indices']];stats=M.box_statistics(bounds,point[1:]);geom=old['geometry'];assert stats['eta_old_Q']==F(geom['eta']) and stats['h_infinity']==F(geom['h_inf']);ends=old['endpoints'];assert [x['t'] for x in ends]==TS and len(ends)==73;B2=F(geom['B2']);hinf=stats['h_infinity'];mean=arb(0);new_bands=0
  for i,x in enumerate(ends):
   a=F(x['t']);nxt=ends[i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'};ml=F(x['H']['lower'])-F(nxt['upper']);mh=F(x['H']['upper'])-F(nxt['lower']);assert 0<ml<=mh;ex=ep(arb(rq(hinf/a)).exp(),True);z=arb(rq(B2/(20*a*a))).sqrt();oldR=ep(((1+z*z).sqrt()+z)**10,True);R=min(ex,oldR);factor=M.band_mean_factor(a,R,stats);mean+=arb(rq(mh*factor));new_bands+=1
  b=old['bounds'];l=F(b['L_i'],N);u=F(b['U_i'],N);oldm=F(b['M_i'],N);mi=min(b['M_i'],ceil(ep(mean,True)*N));m=F(mi,N);assert 0<l<=m<=oldm<=u<=3;q=arb(5)/4;v=arb(rq(stats['volume']));lp=arb(rq(l))**q;up=arb(rq(u))**q;chord=v*(lp+(up-lp)*arb(rq((m-l)/(u-l))));integralU=ceil(ep(chord,True)*N);newU=min(row['old_self_U'],row['weight']*integralU);assert row['old_self_L']<=newU
  r={'status':'complete','fingerprint':FP,'rank':row['rank'],'key':row['key'],'level':level,'indices':row['indices'],'weight':row['weight'],'old_self_L':row['old_self_L'],'old_self_U':row['old_self_U'],'new_self_U':newU,'rangeL_i':b['L_i'],'rangeU_i':b['U_i'],'old_mean_i':b['M_i'],'new_mean_i':mi,'new_chord_integral_upper_i':integralU,'M2':str(stats['E_max_abs_squared']),'eta_old_Q':str(stats['eta_old_Q']),'mu_infinity':str(stats['mu_infinity']),'h_infinity':str(hinf),'positive_bands':new_bands,'source_evidence_sha256':row['source_evidence_sha256'],'new_H_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-began};gz(out,r);return r
 except BaseException as e:
  import traceback
  r={'status':'failed','key':row['key'],'rank':row['rank'],'error':repr(e),'traceback':traceback.format_exc(),'pid':os.getpid()};gz(S/'diagnostics'/f"{row['rank']:05d}.json.gz",r);return r

def recompute(state):
 nodes=state['nodes']
 for k in sorted(nodes,key=lambda k:nodes[k]['level'],reverse=True):
  n=nodes[k]
  if n['children']:n['effU']=min(n['selfU'],sum(nodes[c]['effU'] for c in n['children']));n['effL']=max(n['selfL'],sum(nodes[c]['effL'] for c in n['children']))
  else:n['effU'],n['effL']=n['selfU'],n['selfL']
  assert n['effL']<=n['effU']
def budget(state,reg):
 lo=F(sum(state['nodes'][k]['effL'] for k in state['roots']),N);up=F(sum(state['nodes'][k]['effU'] for k in state['roots']),N);E=F(reg['E076_upper']);g=F(reg['valid_g']);B=ep(arb(rq(up))**(arb(4)/5)+arb(rq(E)),True);I=ep(arb(rq(B))**(arb(5)/4),True);cost=ep((5*arb(rq(B))/arb(rq(g)))**(arb(5)/4),True);assert lo==F(reg['initial_LH']);return {'LH_truncated_integral':str(lo),'UH_truncated_integral':str(up),'UH_display':float(up),'fixed_E076_upper':str(E),'B_upper':str(B),'I_upper':str(I),'I_upper_display':float(I),'valid_g':str(g),'g_acceptance_source_sha256':reg['g_acceptance_sha256'],'cost_upper':str(cost),'cost_upper_display':float(cost),'cost_strictly_below32':5**5*B**5<32**4*g**5,'I_strictly_below_old_complete_cap':I<F(reg['old_complete_I_cap']),'old_complete_I_cap':reg['old_complete_I_cap'],'physical_frontier_count':len(state['frontier']),'physical_frontier_unchanged':True,'no_tube_added':True,'new_H_calls':0}
def main():
 began=time.monotonic();ctx.prec=384;ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');args=ap.parse_args();reg=json.loads((H/'registration.json').read_text());assert sha(__file__)==reg['engine_sha256'] and sha(H/'moment085.py')==reg['moment085_sha256'] and sha(H/'nodes.json.gz')==reg['nodes_sha256'];fp=hashlib.sha256((sha(__file__)+sha(H/'registration.json')).encode()).hexdigest();nodes=json.loads(gzip.decompress((H/'nodes.json.gz').read_bytes()));sp=S/'state.json.gz'
 if sp.exists():state=json.loads(gzip.decompress(sp.read_bytes()))
 else:
  state=json.loads(gzip.decompress((S/'initial_state.json.gz').read_bytes()));assert sha(S/'initial_state.json.gz')==reg['initial_state_sha256'];state['reenvelope088']={'fingerprint':fp,'next_rank':0,'blocks':[]};gz(sp,state)
 meta=state['reenvelope088'];assert meta['fingerprint']==fp and meta['next_rank']<56139,'zero-H phase complete; no extra dispatch';bid=len(meta['blocks']);tag=f'{bid:04d}';start=meta['next_rank'];end=min(56139,start+4096);selected=nodes[start:end];brp=H/('block'+tag+'_registration.json');selection=hashlib.sha256(serial(selected).encode()).hexdigest();before=budget(state,reg)
 if brp.exists():assert args.resume;br=json.loads(brp.read_text());assert br['parent_state_sha256']==sha(sp) and br['selection_sha256']==selection
 else:
  assert not args.resume;br={'id':'upper_reenvelope088.block'+tag,'fingerprint':fp,'registered_unix':time.time(),'start':start,'end':end,'count':end-start,'selection_sha256':selection,'keys':[r['key'] for r in selected],'parent_state_sha256':sha(sp),'parent_receipt_sha256':meta['blocks'][-1]['receipt_sha256'] if meta['blocks'] else None,'source_engine_sha256':sha(__file__),'source_registration_sha256':sha(H/'registration.json'),'worker_count':12,'wallcap':None,'new_H_registered':0,'before':before};atomic(brp,br)
 print('REGISTERED',tag,start,end,flush=True);receipt={'block':bid,'pid':os.getpid(),'fingerprint':fp,'registration_sha256':sha(brp),'status':'running','before':before};rp=H/('block'+tag+'_receipt.json');atomic(rp,receipt);pool=mp.get_context('spawn').Pool(12,initializer=init,initargs=(fp,reg['time_endpoints']));rows=[];fault=None
 try:
  for r in pool.imap_unordered(worker,selected):
   rows.append(r)
   if r['status']=='failed':fault=r;pool.terminate();break
   if len(rows)%128==0:gz(S/('block'+tag+'_checkpoint.json.gz'),rows)
  if fault is None:pool.close()
  pool.join();receipt['workers_joined']=True
 except BaseException as e:pool.terminate();pool.join();receipt.update({'workers_joined':True,'exception':repr(e)})
 complete=len(rows)==len(selected) and all(r['status']=='complete' for r in rows);receipt.update({'status':'complete' if complete else 'failed','fault':fault,'count':len(rows),'worker_pids':sorted({r['pid'] for r in rows}),'wall_seconds':time.monotonic()-began,'new_H_calls':0})
 if not complete:atomic(rp,receipt);print(serial(receipt));return
 rows.sort(key=lambda r:r['rank']);assert [r['rank'] for r in rows]==list(range(start,end));gz(H/('block'+tag+'_cells.json.gz'),rows)
 for r in rows:
  n=state['nodes'][r['key']];assert n['selfL']==r['old_self_L'] and n['selfU']==r['old_self_U'];n['selfU']=r['new_self_U']
 recompute(state);after=budget(state,reg);assert F(after['UH_truncated_integral'])<=F(before['UH_truncated_integral']);assert state['frontier']==json.loads(gzip.decompress((S/'initial_state.json.gz').read_bytes()))['frontier'];receipt.update({'after':after,'improved_nodes':sum(r['new_self_U']<r['old_self_U'] for r in rows),'certificates_sha256':sha(H/('block'+tag+'_cells.json.gz')),'mean_seconds_sum':sum(r['seconds'] for r in rows)});atomic(rp,receipt);meta['next_rank']=end;meta['blocks'].append({'id':bid,'start':start,'end':end,'registration_sha256':sha(brp),'receipt_sha256':sha(rp),'certificates_sha256':receipt['certificates_sha256']});gz(sp,state);atomic(H/'progress.json',{'processed':end,'total':56139,'all_complete':end==56139,'state_sha256':sha(sp),'budget':after});print(serial(receipt),flush=True)
if __name__=='__main__':main()
