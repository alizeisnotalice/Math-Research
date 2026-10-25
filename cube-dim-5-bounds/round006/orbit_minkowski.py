"""Primary finite arithmetic for an analytic orbit/Minkowski I bound."""
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import hashlib,json,time
p=Path(__file__).resolve().parent;start=time.monotonic();ctx.prec=256
beta=F(3,4);q=F(10,7);b=F(3,14)/q;a=5*b/16
assert b==F(3,20) and a==F(3,64) and 5*b-beta==0
assert 3840//160==24 and F(2)**4==16
J=2*(arb(7)/4).gamma()*((arb(17)/16).gamma()/(arb(67)/64).gamma())**16
c=3/(2**(arb(15)/2)*arb.pi()**(arb(5)/2))
I=c**(arb(10)/7)*arb.pi()**2/60*J**(arb(10)/7)
oldQ=arb(16)**(arb(5)/14)*(arb(7)/4).gamma()*(arb(4)/7).gamma()/(arb(19)/28).gamma()*((arb(61)/56).gamma()/(arb(239)/224).gamma())**16
ratio=J**(arb(10)/7)/oldQ
assert ratio<1
needed=arb(10)/3*(I/32)**(arb(7)/10)
assert I<arb(1)/60000
out={'status':'analytic_orbit_minkowski_candidate_primary_arithmetic_passed_pending_independent_audit','b':'3/20','a':'3/64','T_exponent':'0','J_upper_ball':J.str(70,more=True),'I_upper_ball':I.str(70,more=True),'new_over_old_I_upper_ball':ratio.str(70,more=True),'g_needed_for_cost32_ball':needed.str(70,more=True),'simple_rational_I_upper':'1/60000','new_C5_upper':False,'elapsed_seconds':time.monotonic()-start,'precision_bits':ctx.prec,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'registration_sha256':hashlib.sha256((p/'orbit_minkowski_registration.md').read_bytes()).hexdigest()}
(p/'orbit_minkowski.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
