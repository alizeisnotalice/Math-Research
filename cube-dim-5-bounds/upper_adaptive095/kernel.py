"""New same-model point/time H values only, exact zeroth DD and091+085 envelopes."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm,prod,ceil,floor
import gzip,json,hashlib,os,time,resource,sys,traceback
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive095');N=1<<160;TABLE=None;TS=None;CAT=None;FP=None

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def atomic(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(serial(x)+'\n');os.replace(tmp,p)
def gz(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_bytes(gzip.compress(serial(x).encode(),mtime=0));os.replace(tmp,p)
def rq(x):
 x=F(str(x));return fmpq(x.numerator,x.denominator)
def ep(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def path(k,partial=False):return S/('partial' if partial else 'evidence')/(k.replace(':','_').replace(',','-')+'.json.gz')
def init(table,ts,catalog,fp):
 global TABLE,TS,CAT,FP
 ctx.prec=384;TABLE=table;TS=ts;CAT=catalog;FP=fp

def bounds(ends,row):
 g=row['geometry'];h=F(g['h_infinity']);mu=F(g['mu_infinity']);m2=F(g['M2']);lo=arb(0);up=arb(0);mean=arb(0);minmass=None
 for i,e in enumerate(ends):
  a=F(e['t']);nxt=ends[i+1]['H'] if i+1<len(ends) else {'lower':'0','upper':'0'};ml=F(e['H']['lower'])-F(nxt['upper']);mh=F(e['H']['upper'])-F(nxt['lower']);assert 0<ml<=mh,('nonpositive_band',i,str(ml),str(mh));minmass=ml if minmass is None else min(minmass,ml);R=F(CAT[str(h/a)]['upper']);C=min(R,1+mu/a+F(25,48)*m2*R/a**2);lo+=arb(rq(ml/R));up+=arb(rq(mh*R));mean+=arb(rq(mh*C))
 li=floor(ep(lo)*N);ui=min(3*N,ceil(ep(up,True)*N));mi=min(ui,ceil(ep(mean,True)*N));assert 0<li<=mi<=ui,('range_mean_not_separated',li,mi,ui);l,u,m=F(li,N),F(ui,N),F(mi,N);v=arb(rq(g['volume']));q=arb(5)/4;lp=arb(rq(l))**q;upow=arb(rq(u))**q;chord=v*lp if l==u else v*(lp+(upow-lp)*arb(rq((m-l)/(u-l))));return {'L_i':li,'U_i':ui,'M_i':mi,'integral_lower_i':floor(ep(v*lp)*N),'integral_upper_i':ceil(ep(chord,True)*N),'positive_bands':73,'minimum_band_mass_lower':str(minmass),'global_cap_used':ep(up,True)>3}

def worker(row):
 began=time.monotonic();k=row['key'];out=path(k);cp=path(k,True)
 if out.exists():
  r=json.loads(gzip.decompress(out.read_bytes()));assert r['fingerprint']==FP and r['point']==row['point'];return {'key':k,'status':'reused_completed','bounds':r['bounds'],'evidence_sha256':sha(out),'runtime':r['runtime']}
 state={'fingerprint':FP,'key':k,'point':row['point'],'next_cone':0,'weights':{},'endpoints':[],'segments':[]};runtime={'pid':os.getpid(),'began_unix':time.time(),'new_H_calls':0,'reused_H_values':0,'new_DD_cones':0,'new_jet_calls':0}
 try:
  if cp.exists():state=json.loads(gzip.decompress(cp.read_bytes()));assert state['fingerprint']==FP and state['point']==row['point']
  p=list(map(F,row['point']));D=lcm(*(x.denominator for x in p));nodes=[sum(int(x*D)*z for x,z in zip(p,v)) for v in TABLE['vertices']];assert len(set(nodes))==280;weights={int(x):rq(v) for x,v in state['weights'].items()};t=time.monotonic()
  for ci in range(state['next_cone'],3504):
   cone=TABLE['simplices'][ci];ks=[nodes[i] for i in cone['vertices']];det=cone['det_numerator']
   for j,x in enumerate(ks):weights[x*x]=weights.get(x*x,fmpq(0))-fmpq(det,2*prod(x-y for n,y in enumerate(ks) if n!=j))
   state['next_cone']=ci+1;runtime['new_DD_cones']+=1
   if (ci+1)%512==0 or ci==3503:state['weights']={str(x):str(v) for x,v in weights.items() if v};gz(cp,state)
  weights={x:v for x,v in weights.items() if v};assert len(weights)<=140 and sum(weights.values(),fmpq(0))==0;runtime['stencil_seconds']=time.monotonic()-t;t=time.monotonic();factor=arb(D)**4/(144*(2*arb.pi())**5)
  for i in range(len(state['endpoints']),73):
   ts=TS[i]
   if ts in row['reuse_H']:ball=row['reuse_H'][ts];origin='exact_registered_cache';runtime['reused_H_values']+=1
   else:
    a=F(ts);value=factor*sum((arb(v)*arb(rq((144*D*a)**2+x)).log() for x,v in sorted(weights.items())),arb(0));assert value.is_finite() and value>0,('invalid_H',ts,str(value));ball={'lower':str(ep(value)),'upper':str(ep(value,True))};origin='genuine_new_H';runtime['new_H_calls']+=1
   state['endpoints'].append({'t':ts,'H':ball,'origin':origin})
   state['current_runtime']=runtime;gz(cp,state)
  runtime['endpoint_seconds']=time.monotonic()-t;runtime['new_H_calls_this_invocation']=runtime['new_H_calls'];runtime['new_H_calls']=sum(e['origin']=='genuine_new_H' for e in state['endpoints']);runtime['reused_H_values']=sum(e['origin']=='exact_registered_cache' for e in state['endpoints']);t=time.monotonic();b=bounds(state['endpoints'],row);runtime.update({'envelope_seconds':time.monotonic()-t,'seconds':time.monotonic()-began,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss});ev={'fingerprint':FP,'key':k,'level':row['level'],'indices':row['indices'],'weight':row['weight'],'parent':row['parent'],'point':row['point'],'D':D,'weights':{str(x):str(v) for x,v in sorted(weights.items())},'geometry':row['geometry'],'endpoints':state['endpoints'],'bounds':b,'reuse_sources':row['reuse_sources'],'runtime':runtime};gz(out,ev);cp.unlink();return {'key':k,'status':'complete','bounds':b,'evidence_sha256':sha(out),'runtime':runtime}
 except BaseException as e:
  state['diagnostic']={'error':repr(e),'traceback':traceback.format_exc(),'runtime':runtime,'elapsed_seconds':time.monotonic()-began};gz(cp,state);return {'key':k,'status':'failed','error':repr(e),'pid':os.getpid()}
