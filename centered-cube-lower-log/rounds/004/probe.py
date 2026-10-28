from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent
s=1/9;a=10/27;alpha=729/29;gamma=27/29
def C(x,b):
 x=np.abs(x);return np.where(x<=b,2*b**3/3-b*x*x+x**3/2,np.where(x<=2*b,(2*b-x)**3/6,0.))
def F(L):
 x=np.mod(L,1);return 1+alpha*(2*a**3/3+2*s**3-C(x,a)-C(x-1,a)-3*C(x,s)-3*C(x-1,s))
x=np.linspace(1,2,100001);f=F(x);rat=f/(gamma*x*x)
m=np.arange(1,1025);b=4*alpha*((np.sin(np.pi*m*a)/(np.pi*m))**4+3*(np.sin(np.pi*m*s)/(np.pi*m))**4)
xx=np.linspace(.01,4,1001);fs=1+(np.sin(np.pi*xx[:,None]*m)**2)@b
out=dict(fourier_error=float(np.max(abs(fs-F(xx)))),plateau_error=float(np.max(abs(F(np.linspace(10/9,11/9,500))-gamma*np.linspace(10/9,11/9,500)**2))),points=[dict(L=float(t),ratio=float(F(t)/(gamma*t*t))) for t in [1,10/9,11/9,1.25,1.5,1.75,2,3]],near_plateau=[])
for v in [16,64,256,1024,4096]:
 ids=np.flatnonzero(abs(rat-1)<=1/v);out['near_plateau'].append(dict(v=v,left=float(x[ids[0]]),right=float(x[ids[-1]]),length=float(x[ids[-1]]-x[ids[0]])))
(P/'results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
