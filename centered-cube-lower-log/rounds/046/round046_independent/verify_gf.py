"""Independent NumPy broadcast field arithmetic and full target-basis verification."""
from pathlib import Path
import numpy as np,json,hashlib,math,time,argparse
H=Path(__file__).resolve().parent
def clrem(x,p):
 while x.bit_length()>=p.bit_length():x^=p<<(x.bit_length()-p.bit_length())
 return x
def cmul(a,b,p):
 z=0
 for i in range(b.bit_length()):
  if b>>i&1:z^=a<<i
 return clrem(z,p)
def cpower(a,e,p):
 z=1
 for bit in bin(e)[2:]:
  z=cmul(z,z,p)
  if bit=='1':z=cmul(z,a,p)
 return z
def gcd(a,b):
 while b:a,b=b,clrem(a,b)
 return a
def vmul(a,b,d,p):
 a=np.asarray(a,dtype=np.uint64).copy();b=np.asarray(b,dtype=np.uint64).copy();z=np.zeros(np.broadcast_shapes(a.shape,b.shape),dtype=np.uint64);top=np.uint64(1<<d);lo=np.uint64(p^(1<<d));mask=np.uint64((1<<d)-1)
 for _ in range(d):
  z^=a*(b&np.uint64(1));carry=(a>>(d-1));a=((a<<np.uint64(1))&mask)^(carry*lo);b>>=np.uint64(1)
 return z
def vpower(a,e,d,p):
 z=np.ones_like(a,dtype=np.uint64)
 for bit in bin(e)[2:]:
  z=vmul(z,z,d,p)
  if bit=='1':z=vmul(z,a,d,p)
 return z
def digest(x):return hashlib.sha256(x.astype('<u8',copy=False).tobytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,default=H.parent/'round046_gf/MANIFEST.json');ap.add_argument('--output',type=Path,default=H/'GF_VERIFIED.json');args=ap.parse_args();assert not args.output.exists();start=time.time();m=json.loads(args.manifest.read_text());d,b,p,g=m['d'],m['b'],m['modulus'],m['generator'];K=(1<<d)-1;B=1<<b
 assert d==47 and b==11 and p==(1<<47)+(1<<5)+1
 factors=[2351,4513,13264529];assert math.prod(factors)==K
 assert all(all(z%a for a in range(2,math.isqrt(z)+1))for z in factors+[d])
 x=2
 for _ in range(d):x=cmul(x,x,p)
 assert x==2 and gcd(cmul(2,2,p)^2,p)==1
 assert cpower(g,K,p)==1 and all(cpower(g,K//z,p)!=1 for z in factors)
 W=[0]
 for v in m['basis']:W+=[u^v for u in W]
 assert len(W)==len(set(W))==B
 V=np.array(W[1:],dtype=np.uint64);iv=vpower(V,K-1,d,p);assert np.all(vmul(V,iv,d,p)==1)
 # Full outer product independently, then retain unordered multiset indices.
 prods=vmul(V[:,None],V[None,:],d,p);ii,jj=np.triu_indices(B-1);products=np.sort(prods[ii,jj]);del prods,ii,jj
 assert len(products)==B*(B-1)//2 and np.all(products[1:]!=products[:-1]);product_sha=digest(products);del products
 ratios=vmul(V[:,None],iv[None,:],d,p);np.fill_diagonal(ratios,0);ratios=np.sort(ratios[ratios!=0]);assert len(ratios)==(B-1)*(B-2)and np.all(ratios[1:]!=ratios[:-1]);ratio_sha=digest(ratios)
 shifted=np.sort(vmul(g,ratios,d,p));assert len(np.intersect1d(ratios,shifted,assume_unique=True))==0
 gi=cpower(g,K-1,p);short=[g,gi,cmul(g,g,p),cmul(gi,gi,p)];assert not np.isin(short,ratios).any()
 rng=np.random.default_rng(460901);pairs=rng.integers(0,K+1,size=(64,2),dtype=np.uint64);got=vmul(pairs[:,0],pairs[:,1],d,p);expected=np.array([cmul(int(a),int(bb),p)for a,bb in pairs],dtype=np.uint64);assert np.array_equal(got,expected)
 out=dict(status='PASS',d=d,b=b,B=B,modulus=p,generator=g,span_size=len(W),full_product_count=B*(B-1)//2,full_ratio_count=(B-1)*(B-2),product_duplicates=0,ratio_duplicates=0,adjacent_ratio_intersections=0,short_ratio_hits=0,irreducible_exact=True,primitive_exact=True,prime_factors_exact_trial_division=factors,independent_scalar_vector_cases=64,products_sorted_SHA256=product_sha,ratios_sorted_SHA256=ratio_sha,manifest_SHA256=hashlib.sha256(args.manifest.read_bytes()).hexdigest(),script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.time()-start,numpy=np.__version__,scope='Exact integer full basis checks via independent NumPy broadcasts; no q-sized labels. No weak-type log gain claimed.')
 args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
