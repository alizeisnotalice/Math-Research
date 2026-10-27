from fractions import Fraction as Q
from math import comb,log,sqrt,ceil
from pathlib import Path
from scipy.optimize import minimize_scalar
from flint import arb,fmpq,ctx
import json,time,datetime,sys,resource,hashlib,scipy,flint
D=Path(__file__).resolve().parent; ctx.prec=512
start=datetime.datetime.now(datetime.timezone.utc).isoformat(); w=time.perf_counter(); c=time.process_time(); trace=[]
def coeff(t):
 q=1+1/(2*t);h=(1+1/t)/2
 s=sum(Q(comb(7,j))*q**(7-j)*Q(-1,2)**j*(1-t**(-j))/j for j in range(1,8))
 return 1-h**8-4*s,8-4*q**7
def obj(x):
 x=float(x); t=Q.from_float(x); a,b=coeff(t); B=float(a)+float(b)*log(x)
 assert B>0
 y=((1+1/x)**4+sqrt(B))**2
 trace.append({'t_hex':x.hex(),'F_hex':y.hex()});return y
def arbq(x):return arb(fmpq(x.numerator,x.denominator))
def ep(x):
 def f(y):
  m,e=y.man_exp();return str(Q(int(m))*Q(2)**int(e))
 return {'lower':f(x.lower()),'upper':f(x.upper())}
try:
 for x in (6,7,8,9,10):obj(x)
 def bounded(x):
  assert len(trace)-5<80
  return obj(x)
 o=minimize_scalar(bounded,bounds=(2,32),method='bounded',options={'maxiter':80,'xatol':2**-24});assert o.success
 best=min(trace,key=lambda r:float.fromhex(r['F_hex'])); x=float.fromhex(best['t_hex']); k=round(Q.from_float(x)*2**20)
 rows=[]
 for v in (k-1,k,k+1):
  t=Q(v,2**20);a,b=coeff(t);B=arbq(a)+arbq(b)*arbq(t).log();assert B>0
  F=(arbq((1+1/t)**4)+B.sqrt())**2
  rows.append({'t':str(t),'numerator_on_grid':v,'B':ep(B),'F':ep(F)})
 winner=min(rows,key=lambda r:Q(r['F']['upper']));u=Q(winner['F']['upper']);upper=Q(u.numerator*10**10//u.denominator+1,10**10)
 result={'status':'completed','initial_count':5,'optimizer_success':bool(o.success),'optimizer_nfev':int(o.nfev),'optimizer_x_hex':float(o.x).hex(),'best_observed':best,'grid_center_integer':k,'grid_denominator':2**20,'rows':rows,'winner':winner,'strict_rational_upper':str(upper),'scope':'finite candidates only; no assertion of global optimum or C8'}
 (D/'trace.json').write_text(json.dumps(trace,indent=2)+'\n');(D/'result.json').write_text(json.dumps(result,indent=2)+'\n');status='completed';error=None;code=0
except BaseException as e:
 status='failed';error=repr(e);code=1
receipt={'status':status,'error':error,'exit_code':code,'start_utc':start,'finish_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.perf_counter()-w,'cpu_seconds':time.process_time()-c,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version,'flint':flint.__version__,'scipy':scipy.__version__,'precision_bits':512,'objective_calls':len(trace),'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (D/'registration.md',D/'run.py')}}
(D/'actualreceipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));sys.exit(code)
