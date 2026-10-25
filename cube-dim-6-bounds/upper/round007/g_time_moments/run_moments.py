#!/usr/bin/env python3
"""One-center Arb384 H-time nodes, real fan gradients, and box certificates."""
from fractions import Fraction as Q
from functools import cmp_to_key
from pathlib import Path
from time import perf_counter
import hashlib,itertools,json,math,os
from flint import arb,acb,ctx

HERE=Path(__file__).resolve().parent
MAIN=HERE.parents[2]
FAN=MAIN/'upper/round003/direct_kernel/fan_integer_input.txt'
CENTER=MAIN/'upper/round003/direct_kernel/local_minimum/exact_local_reference_results.json'
TAUS=[Q(2**k,1000)for k in range(8)]
ctx.prec=384
def A(x):return arb(x.numerator)/x.denominator
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cmp(a,b):return(a[0]*b[1]>b[0]*a[1])-(a[0]*b[1]<b[0]*a[1])
def eq(a,b):return a[0]*b[1]==b[0]*a[1]
def node(n,d):
 g=math.gcd(n,d);return(n//g,d//g)
def jet(x,e,k):
 if k==0:return(x/e).atan()
 ii=acb(0,1)
 return ((1 if k&1 else -1)*(ii**k/acb(e,x)**k-(-ii)**k/acb(e,-x)**k)/(2*k*ii)).real
def dd(z,e):
 z=sorted(z,key=cmp_to_key(cmp));v=[arb(n)/d for n,d in z]
 b=[jet(x,e,0)for x in v]
 for k in range(1,len(z)):
  b=[jet(v[j],e,k)if eq(z[j],z[j-k])else(b[j-k+1]-b[j-k])/(v[j]-v[j-k])for j in range(k,len(z))]
 return b[0]
def dd_grad(z,e,slopes):
 # z and slope rows have already been sorted together.
 n=len(z);v=[arb(a)/b for a,b in z];d=len(slopes[0])
 if any(eq(z[j-1],z[j])for j in range(1,n)):
  value=dd(z,e);out=[arb(0)for _ in range(d)]
  for j in range(n):
   sens=dd(list(z)+[z[j]],e)
   for i in range(d):out[i]+=sens*slopes[j][i]
  return value,out,True
 b=[(x/e).atan()for x in v]
 db=[[e/(e*e+v[j]*v[j])*slopes[j][i]for i in range(d)]for j in range(n)]
 for k in range(1,n):
  nb=[];nd=[]
  for j in range(k,n):
   denom=v[j]-v[j-k];value=(b[j-k+1]-b[j-k])/denom
   deriv=[(db[j-k+1][i]-db[j-k][i]-value*(slopes[j][i]-slopes[j-k][i]))/denom for i in range(d)]
   nb.append(value);nd.append(deriv)
  b,db=nb,nd
 return b[0],db[0],False
def anchor_check():
 out=[];e=arb(1)/250
 for label,z,ds in [
  ('distinct',[(-2,1),(-1,1),(0,1),(1,1),(2,1),(3,1)],[3,-1,2,0,-2,1]),
  ('confluent',[(-2,1),(-1,1),(0,1),(0,1),(2,1),(3,1)],[3,-1,2,-3,-2,1])]:
  value,g,repeated=dd_grad(z,e,[[arb(x)]for x in ds]);via=arb(0)
  for j in range(6):via+=ds[j]*dd(z+[z[j]],e)
  assert (g[0]-via).contains(0)
  # A separate real derivative of the barycentric formula at distinct nodes.
  direct=None
  if not repeated:
   direct=arb(0)
   xs=[Q(a,b)for a,b in z]
   for j,x in enumerate(xs):
    product=Q(1);logder=Q(0)
    for k,y in enumerate(xs):
     if k!=j:product*=x-y;logder+=Q(ds[j]-ds[k])/(x-y)
    direct+=(ds[j]*e/(e*e+A(x)*A(x))-(A(x)/e).atan()*A(logder))/A(product)
   assert (g[0]-direct).contains(0)
  out.append({'name':label,'DD':value.str(90,more=True),'AD_gradient':g[0].str(90,more=True),
   'copied_node_DD6_gradient':via.str(90,more=True),
   'barycentric_derivative':None if direct is None else direct.str(90,more=True),'PASS':True})
 (HERE/'anchor_calibration.json').write_text(json.dumps(out,indent=2)+'\n')
def fan():
 it=iter(FAN.read_text().split());out=[]
 for _ in range(int(next(it))):
  det=int(next(it));rs=[];alpha=Q(det)
  for _ in range(6):
   q=int(next(it));ray=tuple(int(next(it))for _ in range(6));rs.append((q,ray));alpha/=q
  out.append((A(2*alpha),rs))
 return out
def evaluate(x,den,e,cones,started,gradient=False):
 tic=perf_counter();total=arb(0);grad=[arb(0)for _ in range(5)];terms=0;events=0
 for perm in itertools.permutations(range(6)):
  for mask in range(32):
   signs=[-1 if mask>>j&1 else 1 for j in range(6)]
   y=[signs[j]*x[perm[j]]for j in range(6)]
   for alpha,rs in cones:
    ns=[];ss=[]
    for q,ray in rs:
     ns.append(node(sum(y[j]*ray[j]for j in range(6)),den*q))
     if gradient:
      slope=[arb(0)for _ in range(5)]
      for j in range(6):
       if perm[j]>0:slope[perm[j]-1]=arb(signs[j]*ray[j])/q
      ss.append(slope)
    if gradient:
     order=sorted(range(6),key=cmp_to_key(lambda i,j:cmp(ns[i],ns[j])))
     value,gg,repeated=dd_grad([ns[j]for j in order],e,[ss[j]for j in order])
     total+=alpha*value
     for i in range(5):grad[i]+=alpha*gg[i]
     events+=int(repeated)
    else:total+=alpha*dd(ns,e)
    terms+=2
   if terms%26880==0 and perf_counter()-started>1200:raise TimeoutError('registered 20 minute cap')
 assert terms==967680
 factor=(2*arb.pi())**6
 return {'H_arb384':(total/factor).str(100,more=True),
  'gradient_x2_to_x6':None if not gradient else[(g/factor).str(100,more=True)for g in grad],
  'complete_group_terms':terms,'computed_sign_paired_terms':terms//2,
  'confluent_gradient_terms':events,'elapsed_seconds':perf_counter()-tic}
def bounds(v):
 # lower()/upper() are outward Arb endpoints; retain nonnegative exact mass information.
 return (arb(0)if v.lower()<0 else v.lower(),v.upper())
def bands(rows):
 vals=[arb(row['H_arb384'])for row in rows];out=[]
 for k in range(7):
  delta=vals[k]-vals[k+1];lo,hi=bounds(delta);assert hi>0
  out.append((lo,hi))
 return vals,out
def U2(r,l,vals,masses,details=False):
 terms=[A(TAUS[k])**-2*(A(r)/A(TAUS[k])).exp()*masses[k][1]for k in range(l,7)]
 tail=A(TAUS[7])**-2*(A(r)/A(TAUS[7])).exp()*vals[7].upper()
 value=sum(terms,arb(0))+tail
 return (value,terms,tail)if details else value
def band_lower(r,l,vals,masses):
 return sum(((-A(r)/A(TAUS[k])).exp()*masses[k][0]for k in range(l,7)),arb(0))+(-A(r)/A(TAUS[7])).exp()*vals[7].lower()
def actual_box(x,den,r):
 cc=[Q(v,den)for v in x[1:]]
 hs=[(max(-r,-v),min(r,1-v))for v in cc]
 return cc,hs
def taylor_lower(r,l,x,den,vals,masses,gradient):
 cc,hs=actual_box(x,den,r);linear=arb(0);normsq=Q(0)
 for g,(lo,hi)in zip(gradient,hs):
  interval=arb(A((lo+hi)/2),A((hi-lo)/2));linear+=g*interval
  normsq+=max(abs(lo),abs(hi))**2
 curvature=A(normsq)/32*U2(r,l,vals,masses)
 return vals[l]+linear-curvature,linear,curvature
def largest(predicate):
 # Exact rational radius with 1e-6 resolution, capped at 0.1.
 lo,hi=0,100001
 while hi-lo>1:
  m=(lo+hi)//2
  if predicate(Q(m,1000000)):lo=m
  else:hi=m
 return Q(lo,1000000)
def finalize(out,x,den):
 vals,masses=bands(out['time_nodes']);target=arb(1)/60000
 out['time_bands']=[{'a':str(TAUS[k]),'b':str(TAUS[k+1]),'mass_lower':lo.str(95,more=True),
  'mass_upper':hi.str(95,more=True)}for k,(lo,hi)in enumerate(masses)]
 old=[]
 for k,v in enumerate(vals):
  r=largest(lambda r:v*(-A(r)/A(TAUS[k])).exp()>target)if v>target else Q(0)
  old.append({'epsilon':str(TAUS[k]),'radius':str(r)})
 out['same_center_original_multiplicative']=old
 oldbest=max(Q(d['radius'])for d in old);out['best_original_radius']=str(oldbest)
 out['band_multiplicative_radius']=str(largest(lambda r:band_lower(r,0,vals,masses)>target))
 results=[]
 for row in out.get('gradients',[]):
  l=row['time_index'];gradient=[arb(g)for g in row['gradient_x2_to_x6']]
  r=largest(lambda r:taylor_lower(r,l,x,den,vals,masses,gradient)[0]>target)
  lower,linear,loss=taylor_lower(r,l,x,den,vals,masses,gradient)
  moment,parts,tail=U2(r,l,vals,masses,True)
  center_lower=sum((A(TAUS[k+1])**-2*masses[k][0]for k in range(l,7)),arb(0))
  cc,hs=actual_box(x,den,r)
  results.append({'epsilon':str(TAUS[l]),'radius':str(r),'G_lower_arb384':lower.str(95,more=True),
   'strict_G_above_1_over_60000':bool(lower>target),
   'box_bounds_x2_to_x6':[[str(c+a),str(c+b)]for c,(a,b)in zip(cc,hs)],
   'linear_term_enclosure':linear.str(95,more=True),'curvature_loss_upper':loss.str(95,more=True),
   'J2_center_lower':center_lower.str(95,more=True),'J2_center_upper':U2(Q(0),l,vals,masses).str(95,more=True),
   'J2_whole_box_upper':moment.str(95,more=True),'J2_band_upper_terms':[p.str(95,more=True)for p in parts],
   'J2_tail_upper':tail.str(95,more=True),'Hii_whole_box_upper':(arb(4225)/4224*moment).str(95,more=True),
   'radius_over_best_original':str(r/oldbest)if oldbest else None})
 out['taylor_boxes']=results
 if results:
  out['best_taylor_box']=max(results,key=lambda r:Q(r['radius']))
 out['global_g_target_certified']=False;out['new_C6_bound']=False
 out['uncovered']='The complement of the accepted clipped positive-face boxes and their exact symmetry images; no global lower certificate.'
def main():
 started=perf_counter();anchor_check();center=json.loads(CENTER.read_text());x=center['point_numerators'];den=center['denominator']
 inputs=[FAN,CENTER,MAIN/'audit/round006/report.md',MAIN/'audit/round006/check_time_moments_results.json',
  HERE/'registration.md',HERE/'local_theorem.md',HERE/'../..'/'round006/g/sharp_integrand_ode.md']
 # The last relative proof path is resolved explicitly to avoid clone/main ambiguity.
 inputs[-1]=HERE.parents[1]/'round006/g/sharp_integrand_ode.md'
 out={'status':'RUNNING_REGISTERED_ONE_CENTER_BATCH','pid':os.getpid(),'center_numerators':x,'denominator':den,
  'precision_bits':384,'max_seconds':1200,'max_new_H_nodes':8,'max_same_center_gradients':5,
  'input_sha256':{str(p.resolve()):sha(p.resolve())for p in inputs},'source_sha256':sha(Path(__file__)),
  'time_nodes':[],'gradients':[]}
 (HERE/'input_manifest.json').write_text(json.dumps(out,indent=2)+'\n');cones=fan()
 try:
  for k,tau in enumerate(TAUS):
   row=evaluate(x,den,A(tau),cones,started);row.update({'tau':str(tau),'time_index':k});out['time_nodes'].append(row)
   out['elapsed_seconds']=perf_counter()-started;(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
   print('H',tau,row['H_arb384'][:45],'seconds',row['elapsed_seconds'],flush=True)
  for k in range(5):
   if not arb(out['time_nodes'][k]['H_arb384'])>arb(1)/60000:continue
   row=evaluate(x,den,A(TAUS[k]),cones,started,True);row.update({'tau':str(TAUS[k]),'time_index':k})
   assert(arb(row['H_arb384'])-arb(out['time_nodes'][k]['H_arb384'])).contains(0)
   out['gradients'].append(row);finalize(out,x,den);out['elapsed_seconds']=perf_counter()-started
   (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
   print('gradient',TAUS[k],'seconds',row['elapsed_seconds'],'box_r',out['taylor_boxes'][-1]['radius'],flush=True)
 except TimeoutError as ex:out['stop_reason']=str(ex)
 if len(out['time_nodes'])==8:finalize(out,x,den)
 out['status']='COMPLETE_LOCAL_ONLY'if len(out['time_nodes'])==8 and len(out['gradients'])==5 else'STOPPED_LOCAL_ONLY'
 out['elapsed_seconds']=perf_counter()-started;(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
