"""121 independent finite scalar/weight audit. Does not import producer math."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,ceil
import json,gzip,hashlib,time,os,sys
H=Path(__file__).resolve().parent
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-path112-astra/cube-dim-5-bounds/upper_path112')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ser(x):return json.dumps(x,sort_keys=True,separators=(',',':'),default=str)
def read(p):return json.loads(gzip.decompress(Path(p).read_bytes())) if str(p).endswith('.gz') else json.loads(Path(p).read_text())
def write(p,x):Path(p).write_text(ser(x)+'\n')
def integ(p,l,r):return sum((v*(r**(i+1)-l**(i+1))/(i+1) for i,v in enumerate(p)),F(0))
def poly(p,x):return sum((v*x**i for i,v in enumerate(p)),F(0))
def mul_linear(p,b0,b1):
 out=[F(0)]*(len(p)+1)
 for i,v in enumerate(p):out[i]+=b0*v;out[i+1]+=b1*v
 return out
def ratio224(x):
 n,d=x.numerator,x.denominator;D=1<<224
 # Root the unscaled integer sum-of-squares, then divide by the denominator.
 ends=[]
 for k in (4,16):
  m=((k*d)**2+n*n)*D*D;q=isqrt(m)
  assert q*q<=m<(q+1)**2
  ends.append((F(q+n*D,k*d*D),F(q+1+n*D,k*d*D)))
 return ends[0][0]**2*ends[1][0]**8,ends[0][1]**2*ends[1][1]**8
def main():
 t=time.monotonic();assert not (H/'scalar_results.json').exists()
 audit=read(H/'scalar_registration.json')
 for r in audit['sources']:assert sha(r['path'])==r['sha256'],r['path']
 reg=read(P/'scalar_registration.json');cat=read(P/'ratio_catalog.json');old=read(P/'inherited_ratio_catalog.json')
 assert len(reg['missing_x'])==5772 and len(set(reg['missing_x']))==5772
 assert set(cat)==set(old)|set(reg['missing_x']) and not set(old)&set(reg['missing_x'])
 assert set(reg['all_used_x'])-set(old)==set(reg['missing_x'])
 assert len(set(reg['all_used_x'])&set(old))==184
 for x,r in old.items():
  assert cat[x]==r
  o=r['origin']
  if 'identity' in o:assert x=='0' and r['lower']==r['upper']=='1';continue
  assert sha(o['path'])==o['sha256']
  source=read(o['path'])
  rr=next(z['R_new'] for z in source['ratios'] if z['x']==x) if o.get('kind') else source[x]
  assert F(rr['lower'])==F(r['lower']) and F(rr['upper'])==F(r['upper'])
 digest=hashlib.sha256()
 for i,x in enumerate(reg['missing_x']):
  lo,up=ratio224(F(x));r=cat[x]
  assert F(r['lower'])<=lo<=up<=F(r['upper']),x
  digest.update(ser([x,lo,up]).encode())
  if (i+1)%512==0:write(H/'scalar_checkpoint.json',{'new_roots_done':i+1,'total':5772})
 geos=read(P/'geometries.json.gz');plans=read(P/'plans.json.gz');factors=read(P/'factor_catalog.json')
 assert len(geos)==5 and len(plans)==365
 profiles={};cdfpieces=0
 for gid,g in geos.items():
  off=[tuple(map(F,s)) for s in g['offsets']];pr=g['profile'];knots=sorted({F(0)}|{abs(v) for ab in off for v in ab})
  assert [(F(z['lo']),F(z['hi'])) for z in pr['pieces']]==list(zip(knots,knots[1:]))
  assert F(pr['h'])==knots[-1] and F(pr['mu_infinity'])==max(abs((l+r)/2) for l,r in off)
  pieces=[];lastW=None;M2=F(0);J0=F(0)
  for z in pr['pieces']:
   l,r=F(z['lo']),F(z['hi']);cdf=list(map(F,z['cdf']));W=list(map(F,z['W']));sur=list(map(F,z['survival']))
   assert len(cdf)<=5 and len(W)<=6
   # Five exact interior values determine a polynomial of degree at most four.
   for k in range(1,6):
    x=l+(r-l)*F(k,6);actual=F(1)
    for a,b in off:
     actual*=F(int(x>=abs(a))) if a==b else max(F(0),min(b,x)-max(a,-x))/(b-a)
    assert poly(cdf,x)==actual
   assert sur==[1-cdf[0]]+[-v for v in cdf[1:]]
   assert [F(i)*W[i] for i in range(1,len(W))]==[-v for v in sur]
   if lastW is not None:assert lastW==poly(W,l)
   lastW=poly(W,r);assert poly(W,l)>=0 and lastW>=0
   w=integ(W,l,r);assert w==F(z['weight_integral']) and w>=0
   J0+=w;M2+=2*integ([F(0),*sur],l,r)
   pieces.append({'lo':l,'hi':r,'cdf':cdf,'survival':sur,'W':W,'weight_integral':w})
  assert lastW==0 and J0==F(pr['J0'])==M2/2 and M2==F(pr['M2'])
  profiles[gid]={'offsets':off,'h':knots[-1],'pieces':pieces,'J0':J0,'M2':M2,'mu_infinity':F(pr['mu_infinity'])};cdfpieces+=len(pieces)
 used=set();allpanels=0;factor_digest=hashlib.sha256();seen=set()
 def R(x,which):return F(1) if not x else F(cat[str(x)][which])
 for k,z in enumerate(plans):
  gid=z['geometry'];a=F(z['time']);pr=profiles[gid];assert a>0 and (gid,a) not in seen;seen.add((gid,a))
  cells=[];ab=[];lower=upper=sec=F(0)
  for piece in pr['pieces']:
   lo,hi=piece['lo'],piece['hi'];W=piece['W']
   for j in range(8):
    l=lo+(hi-lo)*F(j,8);r=lo+(hi-lo)*F(j+1,8);w=integ(W,l,r)
    assert w>=0
    if not w:continue
    # Direct polynomial multiplication, independent of producer's first moment subtraction.
    alpha=integ(mul_linear(W,r/(r-l),-1/(r-l)),l,r)
    beta=integ(mul_linear(W,-l/(r-l),1/(r-l)),l,r)
    assert alpha>=0 and beta>=0 and alpha+beta==w
    ab.append({'left':l,'right':r,'alpha':alpha,'beta':beta})
    cells.append({'left':l,'right':r,'x_left':l/a,'x_right':r/a,'weight':w})
    lower+=w*R(l/a,'lower');upper+=w*R(r/a,'upper')
    sec+=alpha*min(R(l/a,'upper'),R(r/a,'upper'))+beta*R(r/a,'upper')
  nodes=sorted({F(0),pr['h']/a}|{c[q] for c in cells for q in ('x_left','x_right')});used.update(map(str,nodes))
  plan={'profile':pr,'a':a,'cells':cells,'ratio_nodes':nodes}
  assert hashlib.sha256(ser(plan).encode()).hexdigest()==z['plan_sha256']
  assert hashlib.sha256(ser(ab).encode()).hexdigest()==z['alpha_beta_sha256']
  assert len(cells)==z['panels'] and sum((c['weight'] for c in cells),F(0))==pr['J0']==F(z['J0'])
  ru=R(pr['h']/a,'upper');oldJ=pr['J0']*ru;J=min(sec,upper,oldJ);lower=max(lower,pr['J0']);assert lower<=J
  f=min(ru,1+pr['mu_infinity']/a+F(25,24)*J/a**2)
  out=factors[gid][z['time']];assert F(out['J_upper'])==J and F(out['exact_factor'])==f and out['plan_sha256']==z['plan_sha256']
  assert F(out['factor_upper'])==F(ceil(f*(1<<192)),1<<192)
  factor_digest.update(ser([gid,a,J,f]).encode());allpanels+=len(cells)
  if (k+1)%32==0:write(H/'plan_checkpoint.json',{'plans_done':k+1,'total':365})
 assert seen=={(g,F(a)) for g in geos for a in reg['time_endpoints']}
 assert used==set(reg['all_used_x']) and allpanels==19272
 result={'status':'pass','pid':os.getpid(),'seconds':time.monotonic()-t,'new_ratio_224_checks':5772,'old_used_ratios_consumed_without_recalculation':184,'all_inherited_catalog_rows':len(old),'new_ratio_enclosure_digest':digest.hexdigest(),'geometries':5,'CDF_pieces':cdfpieces,'geometry_time_plans':365,'positive_panels':allpanels,'factor_digest':factor_digest.hexdigest(),'all_weights_nonnegative':True,'sum_weights_M2_half':True,'all_factor_rounding_exact':True,'new_H_log_exp_calls':0,'code_sha256':sha(__file__),'registration_sha256':sha(H/'scalar_registration.json')}
 write(H/'scalar_results.json',result);print(ser(result),flush=True)
if __name__=='__main__':main()
