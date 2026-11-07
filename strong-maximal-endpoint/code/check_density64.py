"""Independent double-precision finite Gamma mixtures, not interval bounds."""
import json, math, re, sys
from pathlib import Path

def bisect(fun, lo, hi):
    for _ in range(90):
        mid=(lo+hi)/2
        if fun(mid)>0: hi=mid
        else: lo=mid
    return (lo+hi)/2

def number(s):
    return float(re.sub(r'`{1,2}[0-9.]*','',s).replace('*^','e'))

def lse(v):
    top=max(v)
    return top+math.log(math.fsum(math.exp(t-top) for t in v))

mu=bisect(lambda u:u-1-math.exp(-u),1.,2.)
def location(b):
    # The Gaussian-tail equation is decreasing in z.
    z=bisect(lambda z:math.sqrt(mu/b)*z-math.erfc(z/math.sqrt(2))/2,0.,math.sqrt(math.log(b))+2.)
    return b*mu-math.sqrt(mu*b)*z

def finite(b,c,s):
    p=1/(1+math.exp(c)); radius=math.sqrt(110*b)
    lo=max(0,math.floor(b*p-radius));hi=min(b,math.ceil(b*p+radius))
    vals=[]
    for k in range(lo,hi+1):
        g=s-c*k
        if g>0:
            vals.append(math.lgamma(b+1)-math.lgamma(k+1)-math.lgamma(b-k+1)-c*k-b*math.log1p(math.exp(-c))+(b-1)*math.log(g)-g-math.lgamma(b))
    return lse(vals)

inp=Path(sys.argv[1]) if len(sys.argv)>1 else Path('/tmp/strong-maximal-r064/affinity64_1.json')
out=Path(sys.argv[2]) if len(sys.argv)>2 else Path('/tmp/strong-maximal-r064/density_check64.json')
source=json.loads(inp.read_text()); rows=[]
for r in source['rows']:
    b=int(r['b'])
    if b>1048576: continue
    h=number(r['h']);c=mu-h/2;d=mu+h/2;s=location(b)
    fc=finite(b,c,s);fd=finite(b,d,s)
    rows.append({'b':b,'mode':r['mode'],'c_logdensity_error':fc-number(r['log_f_c']),'d_logdensity_error':fd-number(r['log_f_d'])})
# A genuinely independent three-shift multinomial sum at b=256.
r=next(r for r in source['rows'] if r['b']=='256' and r['mode']==0)
b=256;h=number(r['h']);c=mu-h/2;d=mu+h/2;s=location(b);alpha=math.sqrt(2)-1;beta=2-math.sqrt(2);T=1+alpha*math.exp(-c)+beta*math.exp(-d);vals=[]
for k in range(b+1):
    for ell in range(b-k+1):
        g=s-c*k-d*ell
        if g>0:
            vals.append(math.lgamma(b+1)-math.lgamma(k+1)-math.lgamma(ell+1)-math.lgamma(b-k-ell+1)+k*(math.log(alpha)-c)+ell*(math.log(beta)-d)-b*math.log(T)+(b-1)*math.log(g)-g-math.lgamma(b))
fh=lse(vals);ref=number(r['log_density_ratio_factor'])+(number(r['log_f_c'])+number(r['log_f_d']))/2
result={'status':'complete','rows':rows,'max_logdensity_abs_error':max(abs(r[k]) for r in rows for k in ('c_logdensity_error','d_logdensity_error')),'three_shift_b256_logdensity_error':fh-ref,'binomial_mass_cutoff':'2 exp(-220), mass only','scope':'Independent finite sums in Python double precision, not a certificate'}
out.write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
