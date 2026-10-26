#!/usr/bin/env python3
"""One <=4096-node pure mean/chord reenvelope invocation; zero H kernels."""
from pathlib import Path
from fractions import Fraction as F
from math import ceil
import argparse,gzip,json,hashlib,multiprocessing as mp,os,time,resource,sys
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_range092');N=1<<160;FP=None;TS=None;CAT=None

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def atomic(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_text(serial(x)+'\n');os.replace(t,p)
def gz(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_bytes(gzip.compress(serial(x).encode(),mtime=0));os.replace(t,p)
def rq(x):
 x=F(str(x));return fmpq(x.numerator,x.denominator)
def ep(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def init(fp,ts,catalog):
 global FP,TS,CAT
 ctx.prec=384;FP=fp;TS=ts;CAT=catalog
def outpath(rank):return S/'nodes'/f'{rank//1024:02d}'/f'{rank:05d}.json.gz'
def worker(row):
 began=time.monotonic();out=outpath(row['rank'])
 if out.exists():
  r=json.loads(gzip.decompress(out.read_bytes()));assert r['fingerprint']==FP and r['key']==row['key'] and r['packet_sha256']==row['packet_sha256'];return r
 try:
  assert sha(row['packet_path'])==row['packet_sha256'];packet=json.loads(gzip.decompress(Path(row['packet_path']).read_bytes()));assert packet['key']==row['key'] and packet['original_raw_sha256']==row['original_raw_sha256'];point=list(map(F,packet['point']));assert point[0]==1;L=row['level'];offsets=[(F(i,2**L)-x,F(i+1,2**L)-x) for i,x in zip(row['indices'],point[1:])];h=max(abs(x) for span in offsets for x in span);assert h==F(row['h_infinity']);mu=F(row['mu_infinity']);M2=F(row['M2']);ends=packet['endpoints'];assert [e['t'] for e in ends]==TS;lo=arb(0);up=arb(0);mean=arb(0)
  for i,e in enumerate(ends):
   a=F(e['t']);nxt=ends[i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'};ml=F(e['H']['lower'])-F(nxt['upper']);mh=F(e['H']['upper'])-F(nxt['lower']);assert 0<ml<=mh;R=F(CAT[str(h/a)]['upper']);assert R>=1;C=min(R,1+mu/a+F(25,48)*M2*R/a**2);lo+=arb(rq(ml/R));up+=arb(rq(mh*R));mean+=arb(rq(mh*C))
  li=max(row['old_rangeL_i'],int(ep(lo)*N//1));ui=min(row['old_rangeU_i'],3*N,ceil(ep(up,True)*N));mi=min(row['old_M_i'],ui,ceil(ep(mean,True)*N));assert 0<li<=mi<=ui,('range_mean_intersection_failed',li,mi,ui);l,u,m=F(li,N),F(ui,N),F(mi,N);q=arb(5)/4;v=arb(1)/2**(4*L);lp=arb(rq(l))**q;upp=arb(rq(u))**q;chord=v*lp if ui==li else v*(lp+(upp-lp)*arb(rq((m-l)/(u-l))));il=int(ep(v*lp)*N//1);iu=ceil(ep(chord,True)*N);selfL=max(row['old_self_L'],row['weight']*il);selfU=min(row['old_self_U'],row['weight']*iu);assert selfL<=selfU,('self_intersection_failed',selfL,selfU)
  r={'status':'complete','fingerprint':FP,'rank':row['rank'],'key':row['key'],'level':L,'indices':row['indices'],'weight':row['weight'],'old_self_L':row['old_self_L'],'old_self_U':row['old_self_U'],'new_self_L':selfL,'new_self_U':selfU,'old_rangeL_i':row['old_rangeL_i'],'old_rangeU_i':row['old_rangeU_i'],'old_M_i':row['old_M_i'],'rangeL_i':li,'rangeU_i':ui,'mean_i':mi,'new_power_lower_i':il,'new_chord_upper_i':iu,'M2':row['M2'],'mu_infinity':row['mu_infinity'],'h_infinity':row['h_infinity'],'positive_bands':73,'packet_sha256':row['packet_sha256'],'original_raw_sha256':row['original_raw_sha256'],'new_H_calls':0,'old_exp_asinh_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-began};gz(out,r);return r
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
 lo=F(sum(state['nodes'][k]['effL'] for k in state['roots']),N);up=F(sum(state['nodes'][k]['effU'] for k in state['roots']),N);E=F(reg['E076_upper']);g=F(reg['valid_g']);B=ep(arb(rq(up))**(arb(4)/5)+arb(rq(E)),True);I=ep(arb(rq(B))**(arb(5)/4),True);cost=ep((5*arb(rq(B))/arb(rq(g)))**(arb(5)/4),True);assert lo>=F(reg['initial_LH']);return {'LH_truncated_integral':str(lo),'UH_truncated_integral':str(up),'UH_display':float(up),'fixed_E076_upper':str(E),'B_upper':str(B),'I_upper':str(I),'I_upper_display':float(I),'valid_g':str(g),'g_acceptance_source_sha256':reg['g_acceptance_sha256'],'cost_upper':str(cost),'cost_upper_display':float(cost),'cost_strictly_below32':5**5*B**5<32**4*g**5,'I_strictly_below_old_complete_cap':I<F(reg['old_complete_I_cap']),'old_complete_I_cap':reg['old_complete_I_cap'],'physical_frontier_count':len(state['frontier']),'physical_frontier_unchanged':True,'no_tube_added':True,'new_H_calls':0}
def main():
 began=time.monotonic();ctx.prec=384;ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');args=ap.parse_args();reg=json.loads((H/'registration.json').read_text());assert sha(__file__)==reg['engine_sha256'] and sha(H/'ratio091.py')==reg['ratio091_sha256'] and sha(H/'ratio_catalog.json')==reg['catalog_sha256'] and sha(H/'nodes.json.gz')==reg['nodes_sha256'];fp=hashlib.sha256((sha(__file__)+sha(H/'registration.json')).encode()).hexdigest();nodes=json.loads(gzip.decompress((H/'nodes.json.gz').read_bytes()));sp=S/'state.json.gz'
 if sp.exists():state=json.loads(gzip.decompress(sp.read_bytes()))
 else:
  state=json.loads(gzip.decompress((S/'initial_state.json.gz').read_bytes()));assert sha(S/'initial_state.json.gz')==reg['initial_state_sha256'];state['range092']={'fingerprint':fp,'next_rank':0,'blocks':[]};gz(sp,state)
 meta=state['range092'];assert meta['fingerprint']==fp and meta['next_rank']<56139,'zero-H phase complete; no extra dispatch';bid=len(meta['blocks']);tag=f'{bid:04d}';start=meta['next_rank'];end=min(56139,start+4096);selected=nodes[start:end];brp=H/('block'+tag+'_registration.json');selection=hashlib.sha256(serial(selected).encode()).hexdigest();before=budget(state,reg)
 if brp.exists():assert args.resume;br=json.loads(brp.read_text());assert br['parent_state_sha256']==sha(sp) and br['selection_sha256']==selection
 else:
  assert not args.resume;br={'id':'upper_range092.block'+tag,'fingerprint':fp,'registered_unix':time.time(),'start':start,'end':end,'count':end-start,'selection_sha256':selection,'keys':[r['key'] for r in selected],'parent_state_sha256':sha(sp),'parent_receipt_sha256':meta['blocks'][-1]['receipt_sha256'] if meta['blocks'] else None,'source_engine_sha256':sha(__file__),'source_registration_sha256':sha(H/'registration.json'),'worker_count':12,'wallcap':None,'new_H_registered':0,'before':before};atomic(brp,br)
 print('REGISTERED',tag,start,end,flush=True);receipt={'block':bid,'pid':os.getpid(),'fingerprint':fp,'registration_sha256':sha(brp),'status':'running','before':before};rp=H/('block'+tag+'_receipt.json');atomic(rp,receipt);pool=mp.get_context('spawn').Pool(12,initializer=init,initargs=(fp,reg['time_endpoints'],json.loads((H/'ratio_catalog.json').read_text())));rows=[];fault=None
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
  n=state['nodes'][r['key']];assert n['selfL']==r['old_self_L'] and n['selfU']==r['old_self_U'];n['selfL']=r['new_self_L'];n['selfU']=r['new_self_U']
 recompute(state);after=budget(state,reg);assert F(after['UH_truncated_integral'])<=F(before['UH_truncated_integral']);assert state['frontier']==json.loads(gzip.decompress((S/'initial_state.json.gz').read_bytes()))['frontier'];receipt.update({'after':after,'improved_upper_nodes':sum(r['new_self_U']<r['old_self_U'] for r in rows),'improved_lower_nodes':sum(r['new_self_L']>r['old_self_L'] for r in rows),'certificates_sha256':sha(H/('block'+tag+'_cells.json.gz')),'mean_seconds_sum':sum(r['seconds'] for r in rows)});atomic(rp,receipt);meta['next_rank']=end;meta['blocks'].append({'id':bid,'start':start,'end':end,'registration_sha256':sha(brp),'receipt_sha256':sha(rp),'certificates_sha256':receipt['certificates_sha256']});gz(sp,state);atomic(H/'progress.json',{'processed':end,'total':56139,'all_complete':end==56139,'state_sha256':sha(sp),'budget':after});print(serial(receipt),flush=True)
if __name__=='__main__':main()
