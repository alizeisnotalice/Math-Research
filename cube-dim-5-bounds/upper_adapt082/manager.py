#!/usr/bin/env python3
"""One registered finite adaptive batch, then exit with all workers joined."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product
from math import prod,factorial,lcm
import argparse,gzip,json,hashlib,multiprocessing as mp,os,time,resource,sys
from flint import arb,ctx
import kernel as K
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;OLD=H.parent/'upper_global080';S=K.S;N=1<<160

def key(level,i):return str(level)+':'+','.join(map(str,i))
def weight(i):return factorial(4)//prod(factorial(v) for v in Counter(i).values())
def children(level,i):return [(level+1,c) for c in sorted({tuple(sorted((2*a+b for a,b in zip(i,bits)),reverse=True)) for bits in product((0,1),repeat=4)})]
def initial(fp):
 p=json.loads(gzip.decompress((OLD/'initial_cells.json.gz').read_bytes()));roots=json.loads(gzip.decompress((OLD/'roots.json.gz').read_bytes()));nodes={}
 for xs in p['rows']:
  r=dict(zip(p['columns'],xs));i=roots[r['rank']]['indices'];k=key(5,i);lo=r['weight']*r['integral_lower_i'];up=r['weight']*r['integral_upper_i'];nodes[k]={'key':k,'level':5,'indices':i,'weight':r['weight'],'selfL':lo,'selfU':up,'effL':lo,'effU':up,'children':None,'parent':None,'rootid':r['rootid'],'source':'080'}
 return {'fingerprint':fp,'nodes':nodes,'frontier':sorted(nodes),'batches':[],'roots':sorted(nodes)}
def qualify(level,indices,table):
 mid=[F(2*i+1,2**(level+1)) for i in indices];eta=F(1,2**level*4096);cands=[mid]+[[x+eta*b**j for j,x in enumerate(mid)] for b in (F(1,8),F(1,16),F(1,32))];jits=[(-1,-2,-4,-8),(-1,-3,-7,-11),(-2,-5,-11,-19),(-3,-7,-17,-29),(-2,-11,-23,-43),(-5,-13,-31,-47),(-7,-19,-37,-49)]
 cands += [[x+F(v,3351328) for x,v in zip(mid,jit)] for jit in jits];diag=[]
 for ci,tail in enumerate(cands):
  p=[F(1)]+tail;D=lcm(*(x.denominator for x in p));nodes=[sum(int(x*D)*z for x,z in zip(p,v)) for v in table['vertices']];inside=all(F(i,2**level)<x<F(i+1,2**level) for i,x in zip(indices,tail));distinct=len(set(nodes));diag.append({'candidate':ci,'inside':inside,'distinct':distinct,'D':D,'minimum_absolute_projection':min(map(abs,nodes))})
  if inside and distinct==280:return list(map(str,p)),ci,diag
 raise ValueError(('finite_generic_candidates_exhausted',level,indices,diag))
def budget(state):
 ctx.prec=384;roots=[state['nodes'][k] for k in state['roots']];lh=F(sum(n['effL'] for n in roots),N);uh=F(sum(n['effU'] for n in roots),N);E=F(json.loads((OLD/'cap_scalar.json').read_text())['E076_norm']['upper']);B=arb(K.rq(uh))**(arb(4)/5)+arb(K.rq(E));bu=K.ep(B,True);I=arb(K.rq(bu))**(arb(5)/4);g=F(1,10000);cost=(5*arb(K.rq(bu))/arb(K.rq(g)))**(arb(5)/4);bl=K.ep(arb(K.rq(lh))**(arb(4)/5)+arb(K.rq(E)));gate=5**5*bu**5<32**4*g**5;impossible=5**5*bl**5>=32**4*g**5
 frontier=[state['nodes'][k] for k in state['frontier']];vol=sum((F(n['weight'],2**(4*n['level'])) for n in frontier),F(0));assert vol==1;gap=max((n['selfU']-n['selfL'],n['key']) for n in frontier)
 return {'LH_truncated_integral':str(lh),'UH_truncated_integral':str(uh),'LH_display':float(lh),'UH_display':float(uh),'physical_frontier_count':len(frontier),'physical_weighted_volume':str(vol),'physical_raw_U':str(F(sum(n['selfU'] for n in frontier),N)),'physical_raw_L':str(F(sum(n['selfL'] for n in frontier),N)),'fixed_E076_upper':str(E),'B_upper':str(bu),'I_upper':str(K.ep(I,True)),'I_upper_display':float(I),'conditional_g':str(g),'g_independently_qualified_here':False,'conditional_cost_upper':str(K.ep(cost,True)),'conditional_cost_display':float(cost),'conditional_cost_strictly_below32_exact':gate,'gate_exact_form':'5^5*B_upper^5 < 32^4*g^5','fixed_E_budget_lower_diagnostic_B':str(bl),'fixed_E_budget_route_impossible_below32':impossible,'diagnostic_scope':'only fixed E_upper recipe; not actual I or C5 lower','max_unexpanded_gap':str(F(gap[0],N)),'max_unexpanded_gap_key':gap[1],'no_tube_added':True}
def recompute(state):
 def visit(k):
  n=state['nodes'][k]
  if n['children']:
   for c in n['children']:visit(c)
   n['effL']=max(n['selfL'],sum(state['nodes'][c]['effL'] for c in n['children']));n['effU']=min(n['selfU'],sum(state['nodes'][c]['effU'] for c in n['children']))
  else:n['effL'],n['effU']=n['selfL'],n['selfU']
  assert n['effL']<=n['effU'],('hierarchy_inconsistent',k)
 for k in state['roots']:visit(k)
def main():
 began=time.monotonic();ctx.prec=384;ap=argparse.ArgumentParser();ap.add_argument('--parents',type=int,required=True);ap.add_argument('--resume',action='store_true');args=ap.parse_args();assert args.parents in (64,256)
 reg=json.loads((H/'registration.json').read_text());assert K.sha(__file__)==reg['manager_sha256'] and K.sha(H/'kernel.py')==reg['kernel_sha256'];assert K.sha(OLD/'manifest.json')==reg['source080_manifest_sha256'];assert K.sha(OLD/'initial_cells.json.gz')==reg['source080_cells_sha256'];assert K.sha(OLD/'geometry017.json.gz')==reg['table_sha256'];fp=hashlib.sha256((K.sha(__file__)+K.sha(H/'kernel.py')+K.sha(H/'registration.json')).encode()).hexdigest();st=S/'state.json.gz'
 if not st.exists():K.gzatomic(st,initial(fp))
 state=json.loads(gzip.decompress(st.read_bytes()));assert state['fingerprint']==fp;bid=len(state['batches']);assert bid<2,'two finite batches only, report before further work';assert args.parents==(64 if bid==0 else 256);tag=f'{bid:04d}';brp=H/('batch'+tag+'_registration.json');jobsfile=H/('batch'+tag+'_jobs.json.gz');table=json.loads(gzip.decompress((OLD/'geometry017.json.gz').read_bytes()));before=budget(state)
 if brp.exists():
  assert args.resume;br=json.loads(brp.read_text());assert br['parent_state_sha256']==K.sha(st);jobs=json.loads(gzip.decompress(jobsfile.read_bytes()));assert K.sha(jobsfile)==br['jobs_sha256']
 else:
  assert not args.resume;selected=sorted((state['nodes'][k] for k in state['frontier']),key=lambda n:(-(n['selfU']-n['selfL']),n['level'],n['indices']))[:args.parents];assert len(selected)==args.parents;mapping={};jobs={}
  for p in selected:
   assert p['children'] is None;cs=children(p['level'],p['indices']);assert sum(weight(i) for L,i in cs)==16*p['weight'];keys=[]
   for L,i in cs:
    k=key(L,i);assert k not in jobs and k not in state['nodes'];assert tuple(a//2 for a in i)==tuple(p['indices']);point,ci,diag=qualify(L,i,table);jobs[k]={'key':k,'level':L,'indices':list(i),'weight':weight(i),'point':point,'parent':p['key'],'candidate':ci,'qualification':diag,'cached_H1over64':None};keys.append(k)
   mapping[p['key']]=keys
  assert len(jobs)<=16*args.parents
  src=Path(reg['source066_path']);assert K.sha(src)==reg['source066_sha256'];old=json.loads(src.read_text());assert old['epsilon']=='1/64'
  for row in jobs.values():
   pointkey=','.join(row['point']);rec=old['cache'].get(pointkey)
   if rec is not None:assert rec['point']==row['point'];row['cached_H1over64']={'lower':rec['values'][0][0],'upper':rec['values'][0][1]}
  del old;jobs=[jobs[k] for k in sorted(jobs)];K.gzatomic(jobsfile,jobs);br={'id':'upper_adapt082.batch'+tag,'fingerprint':fp,'global_registration_sha256':K.sha(H/'registration.json'),'manager_sha256':K.sha(__file__),'kernel_sha256':K.sha(H/'kernel.py'),'parent_state_sha256':K.sha(st),'parent_receipt_sha256':state['batches'][-1]['receipt_sha256'] if state['batches'] else None,'registered_unix':time.time(),'selected_parents':[{'key':p['key'],'gap_units':p['selfU']-p['selfL']} for p in selected],'children_by_parent':mapping,'parent_count':len(selected),'child_count':len(jobs),'jobs_sha256':K.sha(jobsfile),'new_H_registered':sum(73-(r['cached_H1over64'] is not None) for r in jobs),'cached_H1over64_registered':sum(r['cached_H1over64'] is not None for r in jobs),'selection_rule':'descending exact self weighted integral gap, then level/indices','worker_count':12,'wall_limit':None,'all_four_axes_synchronous':True,'before':before};K.atomic(brp,br)
 print('REGISTERED',tag,'parents',br['parent_count'],'children',br['child_count'],'newH',br['new_H_registered'],'reuse',br['cached_H1over64_registered'],flush=True)
 receipt={'batch':bid,'fingerprint':fp,'pid':os.getpid(),'began_unix':time.time(),'registration_sha256':K.sha(brp),'status':'running','before':before};rp=H/('batch'+tag+'_receipt.json');K.atomic(rp,receipt);pool=mp.get_context('spawn').Pool(12,initializer=K.init_worker,initargs=(table,reg['time_endpoints'],fp));rows=[];fault=None
 try:
  for r in pool.imap_unordered(K.worker,jobs):
   rows.append(r)
   if r['status']=='failed':fault=r;pool.terminate();break
   if len(rows)%64==0:K.gzatomic(S/('batch'+tag+'_partial.json.gz'),rows);print(tag,len(rows),round(time.monotonic()-began,3),flush=True)
  if fault is None:pool.close()
  pool.join();receipt['workers_joined']=True
 except BaseException as e:pool.terminate();pool.join();receipt.update({'workers_joined':True,'exception':repr(e)})
 complete=len(rows)==len(jobs) and all(r['status'] in ('complete','reused_completed') for r in rows);receipt.update({'status':'complete' if complete else 'failed','fault':fault,'completed_children':sum(r['status'] in ('complete','reused_completed') for r in rows),'wall_seconds':time.monotonic()-began,'worker_pids':sorted({r.get('runtime',{}).get('pid',r.get('pid')) for r in rows if r.get('runtime',{}).get('pid',r.get('pid')) is not None}),'peak_parent_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
 if not complete:K.gzatomic(S/('batch'+tag+'_partial.json.gz'),rows);K.atomic(rp,receipt);print(K.serial(receipt));return
 bykey={r['rank']:r for r in rows};compact=[];jobmap={r['key']:r for r in jobs};degrading=[]
 for k,r in sorted(bykey.items()):
  j=jobmap[k];b=r['bounds'];lo=j['weight']*b['integral_lower_i'];up=j['weight']*b['integral_upper_i'];state['nodes'][k]={'key':k,'level':j['level'],'indices':j['indices'],'weight':j['weight'],'selfL':lo,'selfU':up,'effL':lo,'effU':up,'children':None,'parent':j['parent'],'source':'082.batch'+tag};compact.append({'key':k,'level':j['level'],'indices':j['indices'],'weight':j['weight'],'bounds':b,'evidence_sha256':r['evidence_sha256']})
 for pk,cs in br['children_by_parent'].items():
  assert all(c in bykey for c in cs);p=state['nodes'][pk];assert sum(state['nodes'][c]['weight'] for c in cs)==16*p['weight'];raw=sum(state['nodes'][c]['selfU'] for c in cs)
  if raw>p['selfU']:degrading.append(pk)
  p['children']=cs;state['frontier'].remove(pk);state['frontier']+=cs
 state['frontier']=sorted(state['frontier']);recompute(state);after=budget(state);assert F(after['UH_truncated_integral'])<=F(before['UH_truncated_integral']);assert F(after['LH_truncated_integral'])>=F(before['LH_truncated_integral']);K.gzatomic(H/('batch'+tag+'_cells.json.gz'),compact);receipt.update({'after':after,'raw_childsum_worse_parents':degrading,'new_H_count':sum(r['runtime']['new_H_count'] for r in rows),'reused_H1over64_count':sum(r['runtime']['reused_H_count'] for r in rows),'cell_certificates_sha256':K.sha(H/('batch'+tag+'_cells.json.gz')),'stage_seconds_sum':{n:sum(r['runtime'][n] for r in rows) for n in ('stencil_seconds','endpoint_seconds','bounds_seconds','seconds')}});K.atomic(rp,receipt);state['batches'].append({'batch':bid,'registration_sha256':K.sha(brp),'receipt_sha256':K.sha(rp),'certificates_sha256':receipt['cell_certificates_sha256']});K.gzatomic(st,state);K.gzatomic(H/('batch'+tag+'_frontier.json.gz'),state['frontier']);K.atomic(H/'progress.json',{'batch':bid,'state_sha256':K.sha(st),'budget':after,'hierarchy_nodes':len(state['nodes']),'physical_frontier_count':len(state['frontier'])});print(K.serial(receipt),flush=True)
if __name__=='__main__':main()
