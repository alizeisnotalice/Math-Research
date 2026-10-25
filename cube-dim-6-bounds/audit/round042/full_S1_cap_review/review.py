"""F042 independent four-configuration audit; never imports an author writer."""
from pathlib import Path
from fractions import Fraction as Q
import json, hashlib, os, time, platform
from flint import arb, fmpq, ctx
import sympy as S
P=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def B(q):
 q=Q(q); return arb(fmpq(q.numerator,q.denominator))
def fq(x): return Q(str(x.fmpq()))
def bounds(x): return fq(x.lower()),fq(x.upper())
def pack(x):
 a,b=bounds(x);return {'lower_exact':str(a),'upper_exact':str(b),'display':str(x)}
def saved_bounds(x):
 e=int(x['common_binary_exponent']);return Q(int(x['lower_mantissa']))*Q(2)**e,Q(int(x['upper_mantissa']))*Q(2)**e
def overlap(v,s):
 a,b=bounds(v);c,d=saved_bounds(s);assert max(a,c)<=min(b,d);return True

def main():
 if (P/'execution_receipt.json').exists(): raise SystemExit('refuse second execution')
 clock=json.loads((P/'input_manifest.json').read_text());start=time.monotonic()
 receipt={'status':'RUNNING','pid':os.getpid(),'started_unix':time.time(),'kernel_calls':0,'planned_configurations':4}
 (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 try:
  assert time.time()-clock['started_unix']<clock['budget_seconds'], 'audit deadline'
  for f,h in clock['input_sha256'].items(): assert sha(f)==h, f
  source=Path(clock['source_directory']);manifest=json.loads((source/'manifest.json').read_text())
  for f,h in manifest['whitelist_sha256'].items(): assert sha(source/f)==h,f
  for f,h in manifest['external_readonly_sha256'].items(): assert sha(f)==h,f
  saved=json.loads((source/'results.json').read_text());rows={x['id']:x for x in saved['configurations']}
  out=[]
  # Configuration 1: alternate atan branch for tail; formal derivative only.
  t,e=S.symbols('t e', positive=True)
  decomposition=S.factor(4/t**4-16/t**2+16/(t*t+S.Rational(1,4))-1/(t**4*(t*t+S.Rational(1,4))))
  R=4/(3*e**3)-16/e+16*S.pi-32*S.atan(2*e)
  residual=S.factor(S.diff(R,e)+1/(e**4*(e*e+S.Rational(1,4))))
  limit=S.limit(R,e,S.oo)
  assert decomposition==residual==limit==0
  out.append({'id':1,'name':'P_tail_alternate_branch_identity','partial_fraction_residual':str(decomposition),'derivative_residual':str(residual),'infinity_limit':str(limit),'positivity':'positive original integrand; proof_review section 2'})
  # Configuration 2: same single fixed geometry; epsilon is algebraically prescribed.
  ctx.prec=768;pi=arb.pi();rho=Q(1,64)
  D=B(Q(67320683,3251298631680))/(32*pi**6);cb=4*D
  cr=B(Q(44743079,2970))/pi**6
  V=B(Q(4,45))*B(6).sqrt()*pi**2*B(rho)**5
  eta=B(6).sqrt()*B(rho)
  cs=[B(6).sqrt()/2-eta,B(Q(10,3)).sqrt()/2-eta,1/B(2).sqrt()-eta]
  assert all(bounds(c)[0]>0 for c in cs)
  eps=(cb*pi*V)**B(Q(1,5));vp=V**B(Q(5,6))
  for name,v in [('D',D),('Cbig',cb),('Crem',cr),('volume',V),('epsilon_star',eps)]:overlap(v,rows[4][name])
  for v,s in zip(cs,rows[4]['c']):overlap(v,s)
  assert rows[4]['A']==[992,960,928] and rows[4]['rho']=='1/64'
  out.append({'id':2,'name':'fixed_cap_and_analytic_epsilon','rho':str(rho),'volume':pack(V),'D':pack(D),'Cbig':pack(cb),'Crem':pack(cr),'c':[pack(x) for x in cs],'epsilon':pack(eps),'unique_loose_optimizer':'proved analytically in proof_review section 3'})
  # Configuration 3: independently evaluate atan(2eps), not the author's reciprocal.
  short=B(Q(6,7))*cb**B(Q(1,6))*pi**B(Q(-5,6))*eps**B(Q(7,6))
  tailint=4/(3*eps**3)-16/eps+16*pi-32*(2*eps).atan()
  tail=vp*D*tailint;loose=vp*cb/(3*eps**3)
  assert bounds(tailint)[0]>0 and bounds(loose-tail)[0]>0
  for name,v in [('short_norm_upper',short),('tail_integral',tailint),('tail_norm_upper',tail),('loose_tail_norm_upper',loose),('whole_P_norm_upper',short+tail)]:overlap(v,rows[5][name])
  out.append({'id':3,'name':'complete_P_time_integral','short':pack(short),'tail_integral':pack(tailint),'tail':pack(tail),'whole_P':pack(short+tail),'loose_minus_exact_tail':pack(loose-tail)})
  # Configuration 4: positive expression in original c_i, sqrt(A_i), not scaled poles.
  a,b,c=cs;r,s,u=[B(n).sqrt() for n in (992,960,928)]
  jb=pi*(a*s*u+b*r*u+c*r*s)/(2*a*b*c*(a*s+b*r)*(a*u+c*r)*(b*u+c*s))
  bad=vp*cr*jb;total=short+tail+bad
  assert bounds(jb)[0]>0
  for name,v in [('B_density_integral_factor',jb),('whole_B_norm_upper',bad),('complete_cap_norm_upper',total),('cap_integral_Gq_upper',total**B(Q(6,5)))]:overlap(v,rows[6][name])
  target=Q(2118205441,10**15)
  term_upper=sum(bounds(v)[1] for v in (short,tail,bad))
  assert bounds(total)[1]<target and term_upper<target
  out.append({'id':4,'name':'B_full_time_and_cap_sum','JB_original_coordinates':pack(jb),'whole_B':pack(bad),'total_cap_norm_bound':pack(total),'cap_integral_bound':pack(total**B(Q(6,5))),'three_term_upper_sum_exact':str(term_upper),'strict_rational_cap_norm_bound':str(target),'strict_positive_margin':str(target-term_upper),'all_saved_intervals_overlap':True,'whole_face_certified':False})
  assert len(out)==4
  for f,h in clock['input_sha256'].items():assert sha(f)==h,('input drift',f)
  elapsed=time.time()-clock['started_unix'];assert elapsed<600
  result={'status':'PASS_FULL_S1_CAP_ONLY','configurations':out,'precision_bits':768,'source_sha256':sha(__file__),'source_manifest_sha256':sha(source/'manifest.json'),'pid':os.getpid(),'runtime_seconds':time.monotonic()-start,'audit_elapsed_seconds':elapsed,'python':platform.python_version(),'sympy':S.__version__,'kernel_calls':0,'native_calls':0,'matlab_calls':0,'grid_calls':0,'source_writer_calls':0,'whole_face_I_certified':False,'g_floor_certified':False,'new_C6_bound':False}
  (P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
  receipt.update(status='PASS',completed_unix=time.time(),runtime_seconds=time.monotonic()-start,configurations=4,audit_elapsed_seconds=elapsed)
  print(json.dumps({'status':result['status'],'pid':os.getpid(),'runtime_seconds':receipt['runtime_seconds'],'audit_elapsed_seconds':elapsed,'epsilon':str(eps),'P_short':str(short),'P_tail':str(tail),'B_full':str(bad),'cap_total':str(total),'strict_target':str(target)}))
 except BaseException as exc:
  receipt.update(status='FAILED_STOP_NO_RETRY',error=repr(exc),completed_unix=time.time(),runtime_seconds=time.monotonic()-start)
  raise
 finally:
  (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':main()
