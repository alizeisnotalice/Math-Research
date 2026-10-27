"""Independently recheck stored exact rational certificates; no search."""
import sys,json,time,os,hashlib,resource
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
mp.mp.dps=400;p=Path(__file__).resolve();src=p.with_name('c34a_results.json');data=json.loads(src.read_text());results=[];start=time.monotonic();cpu=time.process_time()
for x in data['results']:
 c0=x['certificate']; Y=Q(x['Y']);a=Q(c0['a']);b=Q(c0['b']);c=Q(c0['c']);q=Q(c0['q']);
 assert Q(1,10)*b<=a<=Q(9,10)*b and Q(2,25)<=b<=Q(4,5) and Q(1,200)<=c<=Q(2,5) and (b+c)*Y>=1
 A=c/(4*(1+c/Y));f1=1-4*(b+c)/Y;f2=1-2*(b+c)**2/(1-b/Y);assert f1>0 and f2>0
 r=9*a*Y/32;terms=[Q(1)]
 for k in range(1,129):terms.append(terms[-1]*r/k)
 P=(b-a)*(f1*f2-Q(1,2)/sum(terms))-1/Y-Q(3,8*Y**3)
 assert P>q>0 and A>q and A==Q(c0['KA_lower']) and P==Q(c0['KP_lower'])
 M=lambda z:mp.mpf(z.numerator)/z.denominator
 native=(M(b-a)*(M(f1*f2)-mp.exp(-M(r))/2)-1/M(Y)-3/(8*M(Y)**3))
 assert native>M(P)
 results.append({'Y':int(Y),'q_lt_KA_and_KP_lower':True,'mp_400digits_exponential_diagnostic_pass':True,'q':str(q)})
out={'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'checks':results,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'algebraic inequalities only; no operator theorem claim'}
p.with_name('c34a_check.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',len(results),'exact rational certificates; 400-digit diagnostics')
