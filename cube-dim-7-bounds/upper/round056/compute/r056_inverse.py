#!/usr/bin/env python3
"""R056 inverse-source finite-law moments and exact target-column sampling."""
import os,sys,time,json,hashlib,resource,concurrent.futures,threading,subprocess
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
import mpmath as mp
D=Path(__file__).resolve().parent;PROOF=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md');sha=lambda b:hashlib.sha256(b).hexdigest()
DS=('1/10','1/9','1/8','1/7','1/6','3/16','1/5','7/32','1/4','17/64','9/32','19/64','3/10','5/16','1/3');RATIOS=(4,8,16,64,256,2**20,2**40,2**80);THETAS=(Q(0),Q(1,3),Q(2,3),Q(1))
def reg():
 out={'batch':'R056','objective':'inverse-source weighted face expectation/variance and target-column factor-two diagnostics; no Cn proof','code_sha256':sha(Path(__file__).read_bytes()),'R051_source_sha256':sha(PROOF.read_bytes()),'d_exact':list(DS),'q_rule':'q=1-3d/4+d*i/254, i=0..127','rho_over_R':list(RATIOS),'h_rule':'h=d/1024','R_rule':'R=q+d-h','rho_rule':'rho=(rho/R)*R','n_rule':'ceil(64*rho/q), t=rho/q, zeta=n^-3','g_corrected':'g=R-(kappa-1)rho=d-h+zeta*rho; kappa*rho=rho+q-zeta*rho','r_rule':'r=q+d+theta*h, theta=0,1/3,2/3,1','combinations':15*128*8*4,'unique_source_integrals_expected':15*128*8*2,'tasks':'60 dynamic chunks, each d and 32 q indices; 1024 target combos and 512 quadratures per chunk','workers':15,'precision_dps':80,'stop':'all registered finite tuples, or first detected mismatch/exception stops only affected chunk with exact parameters; no grid extension','monitor':'0.5 sec then 5 sec','registered_unix_time':time.time()};(D/'registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'registered','sha256':sha((D/'registration.json').read_bytes()),'tasks':60,'combinations':out['combinations'],'integrals':out['unique_source_integrals_expected']}),flush=True)
def task(work):
 mp.mp.dps=80;start=time.monotonic();cpu=time.process_time();ds,q0,q1=work;d=Q(ds);h=d/1024;rows=[];samples=0;source_count=0;integrals=0;fail=None;maxerr=mp.mpf(0);maxvar=mp.mpf(0);columnmax=Q(0)
 for i in range(q0,q1):
  q=1-3*d/4+d*Q(i,254);R=q+d-h
  for ratio in RATIOS:
   rho=ratio*R;t=rho/q;n=-(-(64*rho/q).numerator//(64*rho/q).denominator);zeta=Q(1,n**3);g=d-h+zeta*rho;kaprho=rho+q-zeta*rho
   try:
    assert 0<g<R/2 and rho>=4*R and kaprho>rho
    def mpq(v):return mp.mpf(v.numerator)/v.denominator
    Rm=mpq(R);rm=mpq(rho);gm=mpq(g);ep=1-gm/(2*Rm);lp=-mp.log(ep)
    def f(y,p):
     e=1-(gm-y)/(2*Rm);v=-mp.log1p(-(gm-y)/(2*Rm));return v**p/e
    m1=mp.quad(lambda y:f(y,1),[0,gm])/rm;m2=mp.quad(lambda y:f(y,2),[0,gm])/rm;integrals+=2
    exact1=Rm/rm*lp**2;exact2=2*Rm/(3*rm)*lp**3;weight=1+(-2*Rm*mp.log(ep)-gm)/rm
    variance=m2/weight-(m1/weight)**2;err=max(abs(m1-exact1),abs(m2-exact2));maxerr=max(maxerr,err);maxvar=max(maxvar,variance)
    assert err<mp.mpf('1e-60') and variance>=-mp.mpf('1e-65')
    target=[]
    for theta in THETAS:
     r=q+d+theta*h;assert kaprho>r
     threshold=kaprho-r
     zs={Q(0),rho,-rho,min(rho,threshold),-min(rho,threshold)}
     if threshold<rho:
      for shift in (-h/16,h/16):zs.add(max(-rho,min(rho,threshold+shift)));zs.add(max(-rho,min(rho,-threshold+shift)))
     maximum=Q(0)
     for z in zs:
      ell=min(Q(1),(r+kaprho-abs(z))/(2*r));N=1+int(abs(z)<threshold);assert ell>=Q(1,2);ratio_col=Q(N,2)/ell;assert ratio_col<=1;maximum=max(maximum,ratio_col);samples+=1
     target.append({'theta':str(theta),'r':str(r),'max_tested_N_over_2ell':str(maximum) });columnmax=max(columnmax,maximum)
    rows.append({'i':i,'q':str(q),'R':str(R),'rho_over_R':ratio,'rho':str(rho),'n':n,'g':str(g),'kappa_rho':str(kaprho),'pin_ell':mp.nstr(ep,32),'inverse_weighted_mean_phi':mp.nstr(m1/weight,32),'inverse_weighted_variance_phi':mp.nstr(variance,32),'direct_formula_max_abs_error':mp.nstr(err,8),'targets':target});source_count+=1
   except Exception as ex:
    fail={'d':str(d),'i':i,'q':str(q),'R':str(R),'rho':str(rho),'n':n,'g':str(g),'error':repr(ex)};break
  if fail:break
 out={'d':str(d),'q_index_range':[q0,q1-1],'status':'exception_or_candidate_counterexample' if fail else 'completed','failure':fail,'source_states':source_count,'target_combinations':source_count*4,'quadratures':integrals,'column_exact_sample_checks':samples,'max_column_ratio':str(columnmax),'max_formula_abs_error':mp.nstr(maxerr,12),'max_inverse_weighted_variance':mp.nstr(maxvar,32),'results':rows,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'numerical moment diagnostic and exact sampled columns; not an operator theorem or universal integration certificate'}
 name='chunk_d'+str(DS.index(ds))+'_q'+str(q0)+'.json';(D/name).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'chunk','d':ds,'qrange':[q0,q1-1],'status':out['status'],'states':source_count,'integrals':integrals,'pid':out['pid'],'wall':out['wall_seconds']}),flush=True);return {'file':name,'sha256':sha((D/name).read_bytes()),'status':out['status'],'source_states':source_count,'target_combinations':out['target_combinations'],'quadratures':integrals,'column_checks':samples,'pid':out['pid'],'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes']}
def monitor(stop,parent):
 start=time.monotonic();delay=.5
 with (D/'monitor.jsonl').open('w') as f:
  while not stop.wait(delay):
   try:
    text=subprocess.check_output(['ps','-axo','pid=,ppid=,%cpu=,rss=,time=,command='],text=True);ws=[]
    for line in text.splitlines():
     z=line.strip().split(None,5)
     if len(z)==6 and int(z[1])==parent and 'spawn_main' in z[5]:ws.append({'pid':int(z[0]),'cpu_percent':float(z[2]),'rss_kib':int(z[3]),'cpu_time':z[4]})
    row={'wall_seconds':time.monotonic()-start,'worker_count':len(ws),'total_cpu_percent':sum(x['cpu_percent'] for x in ws),'total_rss_kib':sum(x['rss_kib'] for x in ws),'workers':ws}
   except Exception as e:row={'wall_seconds':time.monotonic()-start,'monitor_error':repr(e)}
   f.write(json.dumps(row,separators=(',',':'))+'\n');f.flush();delay=5
if __name__=='__main__':
 if sys.argv[1:]==['--register']:reg();sys.exit()
 assert sys.argv[1:]==['--run'];r=json.loads((D/'registration.json').read_bytes());assert r['code_sha256']==sha(Path(__file__).read_bytes()) and r['R051_source_sha256']==sha(PROOF.read_bytes());chunks=[(d,q,q+32) for d in DS for q in (0,32,64,96)];start=time.monotonic();stop=threading.Event();thread=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);thread.start();print(json.dumps({'event':'start','pid':os.getpid(),'workers':15,'chunks':len(chunks),'registration_sha256':sha((D/'registration.json').read_bytes()),'monitor':str(D/'monitor.jsonl')}),flush=True)
 try:
  with concurrent.futures.ProcessPoolExecutor(max_workers=15) as ex:results=list(ex.map(task,chunks))
 finally:stop.set();thread.join(timeout=10)
 out={'batch':'R056','registration_sha256':sha((D/'registration.json').read_bytes()),'wall_seconds':time.monotonic()-start,'results':results};(D/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','chunks':len(results),'wall_seconds':out['wall_seconds'],'states':sum(z['source_states'] for z in results),'integrals':sum(z['quadratures'] for z in results)}),flush=True)
