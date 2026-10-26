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
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adapt082');SCALE=1<<160;TABLE=None;TS=None;FP=None

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
def evidencepath(rank):return S/'evidence'/(rank.replace(':','_').replace(',','-')+'.json.gz')
def checkpointpath(rank):return S/'partial'/(rank.replace(':','_').replace(',','-')+'.json.gz')
def init_worker(table,ts,fp):
 global TABLE,TS,FP
 ctx.prec=384;TABLE=table;TS=ts;FP=fp
def qform(u):return sum(x*x for x in u)/4-sum(u)**2/36

def geometry(row):
 c=list(map(F,row['point']))[1:];lo=[F(i,2**row['level'])-x for i,x in zip(row['indices'],c)];hi=[F(i+1,2**row['level'])-x for i,x in zip(row['indices'],c)];assert all(l<0<h for l,h in zip(lo,hi));mu=[(l+h)/2 for l,h in zip(lo,hi)];var=[(h-l)**2/12 for l,h in zip(lo,hi)]
 return {'lo':list(map(str,lo)),'hi':list(map(str,hi)),'mu':list(map(str,mu)),'eta':str(qform(mu)+F(2,9)*sum(var)),'B2':str(max(qform(u) for u in product(*zip(lo,hi)))),'h_inf':str(max(map(abs,lo+hi))),'lambda':str(max(map(abs,mu))),'volume':str(F(1,2**(4*row['level'])))}
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
 began=time.monotonic();rank=row['key'];out=evidencepath(rank);cp=checkpointpath(rank)
 if out.exists():
  r=json.loads(gzip.decompress(out.read_bytes()));assert r['fingerprint']==FP and r['point']==row['point'] and r['rootid']==row['key'];return {'rank':rank,'rootid':row['key'],'status':'reused_completed','weight':row['weight'],'bounds':r['bounds'],'evidence_sha256':sha(out),'runtime':r['runtime']}
 state={'fingerprint':FP,'rank':rank,'rootid':row['key'],'point':row['point'],'next_cone':0,'weights':{},'endpoints':[]};runtime={'pid':os.getpid(),'began_unix':time.time(),'new_H_count':0,'reused_H_count':0,'new_jet_count':0}
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
   if t==F(1,64) and row['cached_H1over64'] is not None:ball=row['cached_H1over64'];origin='cached066_H1over64';runtime['reused_H_count']+=1
   else:
    val=factor*sum((arb(w)*arb(rq((144*D*t)**2+k2)).log() for k2,w in sorted(weights.items())),arb(0));assert val.is_finite() and val>0,('invalid_H',ts,str(val));ball={'lower':str(ep(val)),'upper':str(ep(val,True))};origin='new082_scalar';runtime['new_H_count']+=1
   state['endpoints'].append({'t':ts,'H':ball,'origin':origin})
   if (i+1)%8==0 or i+1==len(TS):
    gzatomic(cp,state)
  runtime['endpoint_seconds']=time.monotonic()-stage;stage=time.monotonic();g=geometry(row);b=bounds(state['endpoints'],g);runtime.update({'bounds_seconds':time.monotonic()-stage,'seconds':time.monotonic()-began,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'new_cone_count':3504-start_cone})
  r={'fingerprint':FP,'rank':rank,'rootid':row['key'],'indices':row['indices'],'level':row['level'],'point':row['point'],'weight':row['weight'],'D':D,'geometry':g,'weights':{str(k):str(v) for k,v in sorted(weights.items())},'endpoints':state['endpoints'],'bounds':b,'runtime':runtime};gzatomic(out,r);cp.unlink();return {'rank':rank,'rootid':row['key'],'status':'complete','weight':row['weight'],'bounds':b,'evidence_sha256':sha(out),'runtime':runtime}
 except BaseException as e:
  state['diagnostic']={'error':repr(e),'traceback':traceback.format_exc(),'runtime':runtime,'elapsed_seconds':time.monotonic()-began};gzatomic(cp,state);return {'rank':rank,'rootid':row['key'],'status':'failed','error':repr(e),'pid':os.getpid()}
