#!/usr/bin/env python3
"""New panel identities only. No actual parameters, old toy or J evaluation."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,os,time
import sympy as s

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,x):
    with (HERE/name).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def main():
    assert not (HERE/'symbolic_results.json').exists()
    write('symbolic_registration.json',{'task':'upper_increment139_new_panel_symbolic_identities',
        'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),
        'finite_objects':['panel mass','conditional moment4','affine mean','integrated chord','telescope mass','telescope weighted mean'],
        'actual_J_cases':0,'old_scalar_norm_Gamma_or_toy_runs':0})
    began=time.monotonic();a,b,D,M,z,L,U=s.symbols('a b D M z L U',positive=True)
    mass=(b**3-a**3)/3;moment=3*(b**7-a**7)/(7*(b**3-a**3))
    theta=(moment-a**4)/(b**4-a**4);lam=(z**4-a**4)/(b**4-a**4)
    assert s.simplify(s.integrate(z**2,(z,a,b))-mass)==0
    assert s.simplify(s.integrate(z**6,(z,a,b))/mass-moment)==0
    mean=D+(M-D)*moment
    assert s.simplify(mean-((1-theta)*(D+(M-D)*a**4)+theta*(D+(M-D)*b**4)))==0
    assert s.simplify(s.integrate(z**2*((1-lam)*L+lam*U),(z,a,b))-mass*((1-theta)*L+theta*U))==0
    n=s.symbols('n',integer=True,positive=True);j=s.symbols('j',integer=True,nonnegative=True)
    assert s.simplify(s.summation(((j+1)**3-j**3)/(3*n**3),(j,0,n-1))-s.Rational(1,3))==0
    assert s.simplify(s.summation(((j+1)**7-j**7)/(7*n**7),(j,0,n-1))-s.Rational(1,7))==0
    out={'task':'upper_increment139','status':'PASS_new_panel_symbolic_identities_only','pid':os.getpid(),
        'seconds':time.monotonic()-began,'registration_sha256':sha(HERE/'symbolic_registration.json'),
        'identities':6,'actual_J_cases':0,'actual_root_atoms':0,'actual_panels':0,
        'old_scalar_norm_Gamma_or_toy_runs':0,'new_H_R_gradient_log_exp_Gamma_calls':0}
    write('symbolic_results.json',out);print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
