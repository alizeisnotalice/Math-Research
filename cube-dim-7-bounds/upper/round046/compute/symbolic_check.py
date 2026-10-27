#!/usr/bin/env python3
"""Symbolic 1D primitives and R=1 row-mean cancellation; no new parameter scans."""
import os,time,json,hashlib,resource
from pathlib import Path
import sympy as s
start=time.monotonic();cpu=time.process_time();u,g=s.symbols('u g',positive=True)
E=-2*s.log(1-u/2)-u
assert s.simplify(s.diff(E,u)-u/(2-u))==0 and E.subs(u,0)==0
# integral_0^g E(u) du = (4-2g)log(1-g/2)+2g-g^2/2
F=(4-2*g)*s.log(1-g/2)+2*g-g**2/2
assert s.simplify(s.diff(F,g)-E.subs(u,g))==0 and F.subs(g,0)==0
# r=1, a=1-g: row integral of 2rho*(J-1) across w in [0,1]
row=2*(1-g)*E.subs(u,g)+g*E.subs(u,g)+F-g**2/2
assert s.simplify(row)==0
out={'status':'PASS','identities':['E(u)=integral_0^u v/(2-v)dv','r=1 face mean J=1 by exact algebra'],'domain':'0<g<1/2, rho>=4, R=1','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('symbolic_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
