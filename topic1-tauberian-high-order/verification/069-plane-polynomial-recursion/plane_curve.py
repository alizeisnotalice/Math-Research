"""Exact rational polynomial-curve recursion for the specified 1D min-cap plane model."""
from fractions import Fraction as F
from math import comb
import argparse,json,sys,time,resource
sys.set_int_max_str_digits(0)
def trim(p):
 while len(p)>1 and p[-1]==0:p.pop()
 return p
def add(*ps):
 o=[F(0)]*max(map(len,ps))
 for p in ps:
  for i,c in enumerate(p):o[i]+=c
 return trim(o)
def scale(p,c):return trim([c*x for x in p])
def mul(a,b):
 o=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:o[i+j]+=x*y
 return trim(o)
def diff(p):return trim([i*p[i] for i in range(1,len(p))] or [F(0)])
def integ(p):return [F(0)]+[x/F(i+1) for i,x in enumerate(p)]
def at(p,x):
 a=F(0)
 for c in reversed(p):a=a*x+c
 return a
def reverse(p):
 out=[F(0)]*len(p)
 for i,c in enumerate(p):
  for j in range(i+1):out[j]+=c*comb(i,j)*(-1)**j
 return trim(out)
def integral(p):return at(integ(p),F(1))
def bernstein_nonnegative(p):
 n=len(p)-1
 return all(sum(p[i]*F(comb(k,i),comb(n,i)) for i in range(k+1))>=0 for k in range(n+1))
def audit(pieces,delta):
 return {'r_derivative_nonnegative':all(bernstein_nonnegative(diff(p['r'])) for p in pieces),
         'h_derivative_nonnegative':all(bernstein_nonnegative(diff(p['h'])) for p in pieces),
         'height_nonnegative':all(bernstein_nonnegative(p['h']) for p in pieces),
         'height_le_delta_r':all(bernstein_nonnegative(add(scale(p['r'],delta),scale(p['h'],-1))) for p in pieces)}
def stats(piece):return {'r_degree' :len(piece['r'])-1,'h_degree':len(piece['h'])-1,'P_degree':len(piece['P'])-1}
def moment(pieces):
 total=F(0)
 for piece in pieces:
  r,h=piece['r'],piece['h'];t=add(mul(add([F(1)],r),h),scale(mul(h,h),-1))
  total+=2*integral(mul(t,diff(r)))
 return total
def recurse(pieces,delta):
 total=at(pieces[-1]['P'],F(1)); lows=[]; highs=[]
 for piece in pieces:
  r,h,P=(piece[k] for k in ('r','h','P'))
  q=h
  M=add(scale(mul(q,q),F(1,2)),mul(q,add([F(1)],scale(r,-1))),P)
  newh=scale(add(q,scale(M,-1)),delta)
  lows.append({'r':q,'h':newh})
  q=add([F(1)],scale(r,-1),h)
  M=add([2*total],scale(P,-1),scale(mul(h,h),F(1,2)))
  newh=scale(add(q,scale(M,-1)),delta)
  highs.append({'r':reverse(q),'h':reverse(newh)})
 new=lows+list(reversed(highs)); prefix=F(0);gaps=[]
 for i,p in enumerate(new):
  prim=integ(mul(p['h'],diff(p['r'])));prim[0]=prefix;p['P']=prim
  if i:gaps.append(at(new[i-1]['r'],F(1))==at(p['r'],F(0)) and at(new[i-1]['h'],F(1))==at(p['h'],F(0)))
  prefix=at(prim,F(1))
 assert all(gaps)
 assert at(new[0]['r'],F(0))==0 and at(new[-1]['r'],F(1))==1
 return new
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--delta',default='1/20');ap.add_argument('--depth',type=int,default=3);ap.add_argument('--save-curves',action='store_true');args=ap.parse_args();delta=F(args.delta)
 assert 0<delta<F(1,10)
 start=time.monotonic();piece={'r':[F(0),F(1)],'h':[F(0),delta],'P':[F(0),F(0),delta/2]};pieces=[piece];rows=[]
 for j in range(args.depth+1):
  t=time.monotonic();I=moment(pieces)
  rows.append({'j':j,'piece_count':len(pieces),'moment_I':str(I),'bernstein_audit':audit(pieces,delta),'degree_max':{key:max(stats(p)[key] for p in pieces) for key in ('r_degree','h_degree','P_degree')},'seconds':time.monotonic()-t})
  if j<args.depth:pieces=recurse(pieces,delta)
 output={'delta':str(delta),'depth':args.depth,'method':'exact Fraction polynomial parameter curves, no numerical fitting','pid':__import__('os').getpid(),'elapsed_seconds':time.monotonic()-start,'maxrss_native':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'rows':rows}
 if args.save_curves:output['final_curves']=[{k:[str(x) for x in p[k]] for k in ('r','h','P')} for p in pieces]
 print(json.dumps(output,indent=2))
if __name__=='__main__':main()
