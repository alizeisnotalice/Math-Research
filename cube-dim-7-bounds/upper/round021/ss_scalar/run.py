#!/usr/bin/env python3
"""One 512-bit outward Arb scalar and exact rational elementary comparison."""
from flint import arb,fmpq,ctx
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,time,resource,sys,datetime,traceback
D=Path(__file__).resolve().parent
ctx.prec=512

def f(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Q) else arb(x)
def frac_endpoint(x):
 m,e=x.man_exp();return Q(int(m)*2**int(e)) if e>=0 else Q(int(m),2**(-int(e)))
def main():
 t=time.perf_counter();cpu=time.process_time();start=datetime.datetime.now(datetime.timezone.utc).isoformat()
 lock=json.loads((D/'sourcehash.json').read_text());p=Path(lock['proof_path']);h=lambda:hashlib.sha256(p.read_bytes()).hexdigest();assert h()==lock['proof_sha256']
 a=Q(33,32);b=Q(65,64);cubic=1+Q(7,10)+Q(7,10)**2/2+Q(7,10)**3/6;assert cubic>2
 rational=(1+a**7)*(1+b**7*(1+42*Q(7,10)))
 loose=Q(57,25)*(1+Q(64,57)*(1+42*Q(7,10)))
 assert 1+a**7 < Q(57,25) and b**7 < Q(64,57)
 assert rational<loose==Q(10013,125)<81<128
 val=(1+f(a)**7)*(1+f(b)**7*(1+42*arb(2).log()))
 lo=frac_endpoint(val.lower());hi=frac_endpoint(val.upper());assert lo<hi<rational
 grid=Q(-(-hi.numerator*10**12//hi.denominator),10**12);assert hi<grid and grid<=rational
 assert h()==lock['proof_sha256']
 return {'status':'candidate scalar computed','a':64,'dimension':7,'precision_bits':512,'expression':'[1+(33/32)^7][1+(65/64)^7(1+42 log 2)]','arb_lower_exact':str(lo),'arb_upper_exact':str(hi),'arb_display':val.str(45),'grid_upper_denominator':10**12,'grid_upper':str(grid),'strict_below_grid_upper':True,'cubic_exp_lower_exact':str(cubic),'cubic_exceeds_2':True,'log2_upper':'7/10','rational_upper_using_log_bound':str(rational),'loose_rational_upper':str(loose),'strict_chain':'actual < rational_upper < 10013/125 < 81 < 128','proof_source_unchanged':True,'proof_status':'candidate only, pending independent audit; no C7 assertion','started_utc':start,'wall_seconds':time.perf_counter()-t,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version,'executable':sys.executable}
if __name__=='__main__':
 try: print(json.dumps(main(),indent=2));sys.exit(0)
 except Exception as exc:traceback.print_exc(file=sys.stderr);print(json.dumps({'status':'FAIL','error':repr(exc)}));sys.exit(1)
