from pathlib import Path
import sympy as S,json,os,time,hashlib
root=Path(__file__).parent
reg=json.loads((root/"symbolic175_registration.json").read_bytes())
assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==reg["code_sha256"]
assert hashlib.sha256((root/"primary_secant175.md").read_bytes()).hexdigest()==reg["proof_sha256"]
a=time.monotonic();t=S.symbols("t",nonnegative=True);D=1+t+t**2+t**3;R=1+t**4/D
res=S.cancel(S.diff(R,t)-t**3*(4+3*t+2*t**2+t**3)/D**2)
assert res==0 and R.subs(t,0)==1 and R.subs(t,1)==S.Rational(5,4)
r={"status":"exact_new_slope_monotonicity_identity_passed","pid":os.getpid(),"seconds":time.monotonic()-a,"derivative_residual":str(res),"endpoint0":str(R.subs(t,0)),"endpoint1":str(R.subs(t,1)),"new_true_range_root_or_gain_cases":0,"new_H_mean_old_baseline_calls":0,"sympy_version":S.__version__}
with (root/"symbolic175_result.json").open("x") as f:json.dump(r,f,indent=2)
print(json.dumps(r))
