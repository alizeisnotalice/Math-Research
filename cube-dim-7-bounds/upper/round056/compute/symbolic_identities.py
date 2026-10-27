#!/usr/bin/env python3
"""Symbolic integral primitives; factor-two branch algebra, no new grid."""
import os,time,json,hashlib,resource
from pathlib import Path
import sympy as s
start=time.monotonic();cpu=time.process_time();R,rho,z,pin=s.symbols('R rho z pin',positive=True)
# One edge: substitution z=-log ell gives dw=2R ell dz; two edges / (2rho).
f1=s.integrate(z,(z,0,pin))*2*R/rho;f2=s.integrate(z*z,(z,0,pin))*2*R/rho
assert s.simplify(f1-R*pin**2/rho)==0 and s.simplify(f2-2*R*pin**3/(3*rho))==0
# When |z|<kappa*rho-r, ell=1 and N=2; otherwise N=1 and ell>=1/2.
ell=s.symbols('ell',positive=True);assert s.simplify(2/(2*s.Integer(1))-1)==0 and s.simplify((1/(2*ell)-1)-(1-2*ell)/(2*ell))==0
out={'status':'PASS','identities':['mean_phi_over_ell = R/rho*(-log ellpin)^2','second_phi_over_ell = 2R/(3rho)*(-log ellpin)^3'],'column_argument':'N=2 forces ell=1; N=1 and ell>=1/2 implies N/(2ell)<=1','scope':'elementary symbolic algebra; geometric hypotheses require separate proof','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};Path(__file__).with_name('symbolic_identities.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
