#!/usr/bin/env python3
"""Finite Arb384 Gauss3 functional, exact algebraic DD, relative error audit."""
from collections import Counter
from functools import cmp_to_key
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
from time import perf_counter
from math import comb,factorial
import hashlib,json,os
from flint import arb,acb,ctx
HERE=Path(__file__).resolve().parent
FAN=HERE/'../../round003/direct_kernel/fan_integer_input.txt'
ctx.prec=384

def sgn(a,b):
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if a*b>0:return (a>0)-(a<0)
    z=a*a-15*b*b
    assert z!=0
    return ((a>0)-(a<0))*((z>0)-(z<0))
def cmpn(x,y):return sgn(x[0]*y[2]-y[0]*x[2],x[1]*y[2]-y[1]*x[2])
def eq(x,y):return x[0]*y[2]==y[0]*x[2] and x[1]*y[2]==y[1]*x[2]
def av(x):return (arb(x[0])+x[1]*arb(15).sqrt())/x[2]
def jet(x,e,k):
    if not k:return (x/e).atan()
    i=acb(0,1)
    z=(i**k/acb(e,x)**k-(-i)**k/acb(e,-x)**k)/(2*k*i)
    return (1 if k&1 else -1)*z.real

def kernel(x,e,cones):
    orbit=Counter()
    for perm in set(permutations(x)):
      for mask in range(64):
        y=tuple(tuple((-1 if mask>>j&1 else 1)*v for v in perm[j]) for j in range(6))
        orbit[y]+=1
    # Unique-permutation reduction needs the stabilizer multiplier.
    stab=factorial(6)//len(set(permutations(x)))
    total=arb(0)
    for y,mult in orbit.items():
      for det,rays in cones:
        nodes=[];alpha=arb(det)
        for den,ray in rays:
          nodes.append((sum(y[j][0]*ray[j]for j in range(6)),sum(y[j][1]*ray[j]for j in range(6)),10*den));alpha/=den
        z=sorted(nodes,key=cmp_to_key(cmpn));v=list(map(av,z));a=[[arb(0)for _ in range(6)]for _ in range(6)]
        for j in range(6):a[j][0]=jet(v[j],e,0)
        for k in range(1,6):
          for j in range(k,6):a[j][k]=jet(v[j],e,k)if eq(z[j],z[j-k])else(a[j][k-1]-a[j-1][k-1])/(v[j]-v[j-k])
        total+=stab*mult*alpha*a[5][5]
    assert sum(orbit.values())*stab==46080
    return total/(2*arb.pi())**6,len(orbit)

def derivative_constants():
    q=F(6,5);a=[F(0)]+[F(comb(63+k,k),32**k)for k in range(1,7)]
    p=[F(1)]+[F(0)]*6;co=[F(0)]*7;b=F(1)
    for j in range(1,7):
      p=[sum(p[i]*a[k-i]for i in range(k+1))for k in range(7)];b*=(q-j+1)/j
      for k in range(7):co[k]+=abs(b)*p[k]
    out=[]
    for n in range(1,4):
      d=co[2*n]*factorial(2*n);c=F(1,2**(2*n)*factorial(2*n));rho=c*d
      assert rho*2**(2*n)>1
      out.append({'n':n,'D_2n':str(d),'Peano_sup_upper':str(c),'rho':str(rho),'theta_at_h_equals_2epsilon':str(rho*2**(2*n))})
    return out

def main():
    start=perf_counter();tok=iter(FAN.read_text().split());cones=[]
    for _ in range(int(next(tok))):
      det=int(next(tok));rays=[]
      for _ in range(6):rays.append((int(next(tok)),tuple(int(next(tok))for _ in range(6))))
      cones.append((det,rays))
    fanhash=hashlib.sha256(FAN.read_bytes()).hexdigest()
    assert fanhash=='7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'
    out={'status':'running','pid':os.getpid(),'epsilon':'1/1000','q':'6/5','geometric_nodes':243,'cells':1,'rows':[],'fan_sha256':fanhash,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'derivative_constants':derivative_constants()}
    e=arb(1)/1000;q=arb(6)/5;Q=arb(0);count=0
    for lo in range(6):
      for mid in range(6-lo):
        hi=5-lo-mid
        if perf_counter()-start>1150:out['status']='STOP_WALL_LIMIT';break
        x=((10,0),)+((5,-1),)*lo+((5,0),)*mid+((5,1),)*hi
        value,orb=kernel(x,e,cones);assert value>0
        multiplicity=factorial(5)//(factorial(lo)*factorial(mid)*factorial(hi))
        weight=arb(multiplicity)*(arb(5)/18)**(lo+hi)*(arb(4)/9)**mid
        Q+=weight*value**q;count+=multiplicity
        row={'counts_low_mid_high':[lo,mid,hi],'multiplicity':multiplicity,'weight':str(weight),'H_arb384':value.str(95,more=True),'orbit_images':orb,'elapsed_seconds':perf_counter()-start}
        out['rows'].append(row);out['quadrature_Q_interval']=Q.str(95,more=True);out['elapsed_seconds']=perf_counter()-start
        (HERE/'gauss3_relative_results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(row),flush=True)
      else:continue
      break
    if len(out['rows'])==21:
      assert count==243
      out['status']='PASS_INTERVAL_QUADRATURE_FUNCTIONAL_REMAINDER_NOT_CLOSED'
      out['formal_Q_norm_interval']=Q**(arb(5)/6)
      out['formal_Q_norm_interval']=str(out['formal_Q_norm_interval'])
      out['Q_exceeds_target_I']=bool(Q>(arb(1)/16000)**q)
      out['theta_global_Gauss3']=str((arb(73354040202609)/16384000000/46080)*1000**6)
      out['strict_relative_envelope_cell_lower_bound']='31250000000000'
      out['gap']='The relative Peano remainder has theta>=1 on the global cell; Q is not I. No certified I target or C6 improvement.'
    out['elapsed_seconds']=perf_counter()-start
    (HERE/'gauss3_relative_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k!='rows'},indent=2))
if __name__=='__main__':main()
