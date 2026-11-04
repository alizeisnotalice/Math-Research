import numpy as np
from scipy.integrate import quad

def gain(x,k):
 y=1-x;m=min(x,y);aa=[-2*m,2*m]
 for a in [2*(x-k),2*(y-k)]:
  if -2*m<a<2*m:aa.append(a)
 return max(a*a/4+a/2+x*y+2*k*k-(k+max(0,k-x+a/2))*(k+max(0,k-y+a/2)) for a in aa)
for k in np.linspace(0,2,17):
 q,e=quad(lambda x:gain(x,k),0,.5,epsabs=1e-12,points=[v for v in [k/2,k,1-k] if 0<v<.5]);print(k,8*q,e)
