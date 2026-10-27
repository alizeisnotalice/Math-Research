import sys,json,hashlib,time,os,resource
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
import mpmath as mp
root=Path(__file__).resolve().parents[3]; sys.path.insert(0,str(root/'revision_v3/numeric'))
from verify import certify,finite_certify
started=time.monotonic(); cpu=time.process_time(); mp.mp.dps=85
n=7
def g(t):
 y=mp.log(t); h=(1+1/t)/2; ell=lambda s: 1+1/(2*t)-mp.exp(-s)/2
 B=1-h**n+n*mp.quad(lambda s:1-ell(s)**(n-1)/2,[0,y]); J=mp.quad(lambda s:ell(s)**(n-2),[0,y]); bp=n/(2*t)+n*h**(n-1)/(2*t*t)+n*(n-1)*J/(4*t*t)
 return bp/mp.sqrt(B)-n/t**2*(1+1/t)**(mp.mpf(n)/2-1)
troot=mp.findroot(g,(mp.mpf('6.46'),mp.mpf('6.48'))); t=Q(808656059921,125000000000)
_,bl,bu,fl,fu=certify(n,t); finite=finite_certify(n,t,bu)[-1]; U=Q(12277480262982331,500000000000000)
r=json.loads((root/'upper/round028/cap_candidate/result.json').read_text()); m,e=r['candidate_interval']['lower']; old=Q(int(m))*Q(2)**int(e)
assert fl<=fu<U<old and finite<U
out={'raw_root_60_digits_numerical_not_certified':mp.nstr(troot,60),'root_residual':mp.nstr(g(troot),12),'rational_t':str(t),'new_F_lower':str(fl),'new_F_upper':str(fu),'strict_U':str(U),'round028_original_full_lower_dyadic':r['candidate_interval']['lower'],'new_upper_below_original_full_lower':True,'U_below_original_full_lower':True,'exact_gap_original_lower_minus_new_upper':str(old-fu),'finite_a':'2^256','finite_bound_below_U':True,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-started,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('n7_against_round028.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
