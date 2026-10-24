"""Nine frozen scalar tail bounds; no full generator, LP or quadrature."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json,sys,time
import flint
from flint import arb,fmpq,ctx
import mpmath as mp
import sympy as sy
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'tail_result.json'
if OUT.exists():raise SystemExit('Existing evidence; no automatic repetition.')
ctx.prec=192;mp.mp.dps=70
start=time.monotonic()
result={'status':'started','arb_bits':192,'mpmath_dps':70,'counts':{'symbolic':0,'Arb_started':0,'Arb_completed':0,'mpmath_started':0,'mpmath_completed':0},'rows':[], 'scope':'Tail lower bound only; no full generator positivity or improved weak-type bound.'}
def save():OUT.write_text(json.dumps(result,indent=2)+'\n')
def A(q):
 q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def M(q):
 q=Q(q);return mp.mpf(q.numerator)/q.denominator
def pw(x,q):return (x.log()*A(q)).exp()
def guard():
 if time.monotonic()-start>=60:raise TimeoutError('Frozen 60-second budget.')
def ball(x):return x.str(65)
save()
try:
 z,p=sy.symbols('z p',positive=True)
 f=(p-1)-p*z+z**p
 assert sy.simplify(sy.diff(f,z,2)-p*(p-1)*z**(p-2))==0
 result['counts']['symbolic']+=1
 # Let X=r^-alpha and q=R^-alpha, then a=X-q and r^-beta=X^p.
 X,q=sy.symbols('X q',positive=True)
 D=-(X**p-q**p)+p*X**(p-1)*(X-q)
 assert sy.simplify(D-((p-1)*X**p-p*q*X**(p-1)+q**p))==0
 result['counts']['symbolic']+=1;save()
 alpha,R=Q(47,40),Q(33,20)
 for beta in [alpha,Q(6,5),Q(5,4)]:
  for r in [Q(1,4),Q(1),Q(8,5)]:
   guard();power=beta/alpha
   row={'beta':str(beta),'r':str(r),'status':'started'};result['rows'].append(row)
   result['counts']['Arb_started']+=1;save()
   if beta==alpha:T=arb(0);L=arb(0);difference=arb(0)
   else:
    za=pw(A(r/R),alpha)
    fz=A(power-1)-A(power)*za+pw(za,power)
    factor=A(4)/A(alpha)*pw(A(r),-beta)*pw(A(R-r),-alpha)
    T=factor*fz;L=factor*A(power*(power-1)/2)*(1-za)**2
    difference=T-L
   assert all(x.is_finite() for x in [T,L,difference])
   # Save outward displayed intervals; assess the actual saved balls.
   Tb,Lb,Db=map(lambda x:arb(ball(x)),[T,L,difference])
   if beta==alpha:passed=(Tb==0 and Lb==0)
   else:passed=(Tb.lower()>0 and Lb.lower()>0 and Db.lower()>0)
   row.update(T=ball(T),T_quadratic=ball(L),difference=ball(difference),strict_tail_checks=bool(passed))
   result['counts']['Arb_completed']+=1;save()
   result['counts']['mpmath_started']+=1;save()
   if beta==alpha:reference=mp.mpf(0)
   else:
    am,bm,rm,Rm,pm=map(M,[alpha,beta,r,R,power])
    # Separate direct D0 formula, not the z-form used by Arb.
    direct=(pm-1)*rm**(-bm)-pm*Rm**(-am)*rm**(am-bm)+Rm**(-bm)
    reference=4/am*direct*(Rm-rm)**(-am)
   # Floating comparison is diagnostic only, independent of rigorous signs above.
   # Arb midpoint string can be bracketed; use midpoint double solely for the diagnostic comparison.
   midpoint=mp.mpf(float(T.mid()))
   error=abs(midpoint-reference)
   row.update(mpmath_reference=mp.nstr(reference,65),float_midpoint_abs_difference=mp.nstr(error,10),
              numeric_agreement=bool(error<=mp.mpf('1e-13')*max(1,abs(reference))),status='completed')
   result['counts']['mpmath_completed']+=1;save()
 result['status']='nine_tail_bounds_and_two_identities_closed'
 result['all_strict_checks']=all(r['strict_tail_checks'] for r in result['rows'])
 result['all_numeric_agreements']=all(r['numeric_agreement'] for r in result['rows'])
except Exception as e:result.update(status='failed_or_interrupted',reason=f'{type(e).__name__}: {e}')
result.update(elapsed_seconds=time.monotonic()-start,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
 versions={'python':sys.version,'flint':flint.__version__,'mpmath':mp.__version__,'sympy':sy.__version__})
save();print(json.dumps({'status':result['status'],'counts':result['counts'],'all_strict_checks':result.get('all_strict_checks')},indent=2))
