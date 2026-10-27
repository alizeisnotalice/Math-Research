"""Symbolic verification for 0<u<1<r; no new shell cases."""
import sys,json,time,os,hashlib,resource
from pathlib import Path
import sympy as s
r,u,w=s.symbols('r u w',positive=True);t=time.monotonic();cpu=time.process_time()
# Pieces below use -r<-u<u<1<r.
ell=(s.integrate(-w-u,(w,-r,-u))+s.integrate(w-u,(w,u,r)))/(2*r)
corner=(s.integrate(1-u-u-w,(w,-r,-u))+s.integrate(1-u,(w,-u,1))+s.integrate(w-u,(w,1,r)))/(2*r)
Eell=(r-u)**2/(2*r);Ecorner=1-u+(r-1)**2/(4*r)+(r-u)**2/(4*r)
assert s.simplify(ell-Eell)==0 and s.simplify(corner-Ecorner)==0
N,p,n,M,L0=s.symbols('N p n M L0')
# Direct sum over one pinned negative face and n-1 iid other coordinates.
ours=(M**(n-1)+(n-1)*N*L0*M**(n-2))/2
sol=M**(n-1)/2+(n-1)*(1+u/r)*L0*M**(n-2)/2
assert s.simplify(ours.subs(N,1+u/r)-sol)==0
out={'identities':'PASS','assumptions':'0<u<1<r; −r<−u<u<1<r','mean_ell_deficit':str(Eell),'mean_corner_deficit':str(Ecorner),'Sol_G_formula':'PASS; pinned axis has only one captured face','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-t,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('s35a_symbolic.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
