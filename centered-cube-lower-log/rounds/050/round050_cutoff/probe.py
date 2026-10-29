from pathlib import Path
import numpy as np,json,hashlib
H=Path(__file__).resolve().parent;Ls=[1.,7/6,11/9];dx=.0001
def F(x):
 k=np.floor(x);r=x-k;return k+1.5*np.minimum(r,.5)+.5*np.maximum(r-.5,0)
def weak(y,step):
 z=np.sort(y);return float(np.max(z*np.arange(len(z),0,-1)*step))
u=(np.arange(10000)+.5)*dx;per=np.max([(F(u+L/2)-F(u-L/2))/L for L in Ls],axis=0);wp=weak(per,dx);rows=[]
for N in [4,8,16]:
 x=-N-1+(np.arange(round((2*N+2)/dx))+.5)*dx
 z=np.max([(F(np.clip(x+L/2,-N,N))-F(np.clip(x-L/2,-N,N)))/L for L in Ls],axis=0)
 bd=float(z[np.abs(x)>N-1].sum()*dx/(2*N));wc=weak(z,dx)/(2*N)
 rows.append(dict(n=1,N=N,J=len(Ls),periodic_weak_grid=wp,normalized_cutoff_weak_grid=wc,normalized_boundary_L1_grid=bd,boundary_candidate_upper=2*len(Ls)/N,weak_difference=abs(wp-wc),difference_candidate_upper=3*len(Ls)/N))
out=dict(rows=rows,dx=dx,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Actual one-dimensional step-density window quadrature, no outward integral certificate; general cutoff estimates proved separately.')
(H/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
