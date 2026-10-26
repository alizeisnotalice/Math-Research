#!/usr/bin/env python3
"""One preregistered new diagnostic,not159's five moment identities."""
import argparse,hashlib,json,os,time
from pathlib import Path
import sympy as s

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--registration',type=Path,required=True);ap.add_argument('--sha',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    assert sha(a.registration)==a.sha and not a.output.exists();r=json.loads(a.registration.read_text())
    assert r['case']=='capped_endpoint_upper_polynomial_first_second_derivatives_only'
    for p,h in {**r['source_locks'],**r['code_locks']}.items():assert sha(p)==h
    began=time.monotonic();b,u,B,V,mu,k=s.symbols('b u B V mu k',nonnegative=True);d=s.symbols('d',positive=True)
    z=b+d*u;upper=1+mu*z+k*z*z*(B+V*u)
    first=d*mu+k*(2*d*z*(B+V*u)+z*z*V)
    second=k*(2*d*d*(B+V*u)+4*d*z*V)
    e1=s.expand(s.diff(upper,u)-first);e2=s.expand(s.diff(upper,u,2)-second);assert e1==e2==0
    variables=(b,d,u,B,V,mu,k)
    assert all(c>=0 for c in s.Poly(s.expand(first),*variables).coeffs())
    assert all(c>=0 for c in s.Poly(s.expand(second),*variables).coeffs())
    result=dict(status='PASS_new_capped_upper_monotone_convex_derivative_identity',pid=os.getpid(),seconds=time.monotonic()-began,registration_sha256=a.sha,
        first_derivative=str(first),second_derivative=str(second),residuals=[str(e1),str(e2)],nonnegative_monomial_coefficients=True,
        real_scalar_factor_mean_corner_H_R_J_calls=0,old159_five_group_replays=0,old150_111_baseline_replays=0,
        scope='symbolic normalized-u upper polynomial only;no time moments or physical factor data')
    with a.output.open('x') as f:json.dump(result,f,sort_keys=True,indent=2);f.write('\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
