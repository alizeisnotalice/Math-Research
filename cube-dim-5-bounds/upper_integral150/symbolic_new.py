#!/usr/bin/env python3
"""New inverse-time weighted moments and safe-chord toy, no actual factors."""
from datetime import datetime,timezone
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,time
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,x):
    with (HERE/n).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def main():
    assert not (HERE/'symbolic_results.json').exists()
    write('symbolic_registration.json',{'task':'upper_integral150_new_symbolic_cases','pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),
        'finite_cases':['inverse-time weighted moment','z2R convexity derivative','tail theta4/5','Chebyshev density t^-3 toy','powerfree chord decrement toy'],
        'new_actual_H_R_factor_mean_or_UH_calls':0,'old_baseline_runs':0})
    began=time.monotonic();a,b,t,z,r=s.symbols('a b t z r',positive=True);R=s.Function('R')
    i0=(a**-4-b**-4)/4;i1=(a**-5-b**-5)/5
    assert s.simplify(i1/i0-s.Rational(4,5)*(a**-5-b**-5)/(a**-4-b**-4))==0
    assert s.simplify(s.diff(z*z*R(r*z),z,2)-(2*R(r*z)+4*r*z*s.Subs(s.diff(R(t),t),t,r*z)+r*r*z*z*s.Subs(s.diff(R(t),t,2),t,r*z)))==0
    za=s.symbols('za',positive=True);assert s.simplify((s.Rational(4,5)*za**5/za**4)/za)==s.Rational(4,5)
    weighted=F(7,24)/F(3,8);reference=F(31,160)/F(15,64);assert weighted==F(7,9)<reference==F(62,75)
    L,U=F(1,16)**4,F(1,8)**4;fL,fU=F(1,16)**5,F(1,8)**5
    oldm,newm=(L+U)/2,(3*L+U)/4;slope=(fU-fL)/(U-L)
    assert slope==F(31,240)>=F(1,16) and slope*(oldm-newm)>=F(1,16)*(oldm-newm)
    out={'task':'upper_integral150','status':'PASS_new_time_moment_and_chord_symbolic_cases','pid':os.getpid(),'seconds':time.monotonic()-began,'registration_sha256':sha(HERE/'symbolic_registration.json'),
        'inverse_time_tail_weight':'4/5','Chebyshev_toy':[str(weighted),str(reference)],'new_chord_toy_slope':str(slope),
        'new_actual_H_R_factor_mean_or_UH_calls':0,'old_baseline_runs':0,'global_Cn_claim':False}
    write('symbolic_results.json',out);print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
