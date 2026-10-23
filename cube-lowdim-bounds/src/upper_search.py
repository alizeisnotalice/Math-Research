"""Repeatable floating discovery of kernel parameters; NOT a certificate."""
import math,json
from pathlib import Path
from scipy.optimize import minimize_scalar
ROOT=Path(__file__).resolve().parents[1]
def J(s):return .5-2/s-4/s**2*math.log1p(-s/2)
def T(r):return 2*(2-r)*(1-math.log(2)+math.log(r))-2*(1-r)*math.log1p(-r)
def at_R(R):
    inner=minimize_scalar(lambda r:J(2*r/R)/T(r),bounds=(.6,.98),method='bounded',options={'xatol':1e-14})
    e=inner.fun/(R*(R-1))
    return {'R':float(R),'epsilon_candidate':float(e),'r_candidate':float(inner.x),'cost_candidate':float(R*R/(R-1)-e/6)}
if __name__=='__main__':
    opt=minimize_scalar(lambda R:at_R(R)['cost_candidate'],bounds=(1.75,2.05),method='bounded',options={'xatol':1e-13})
    out={'fixed_radius':at_R(3797/2000),'joint':at_R(opt.x),'status':'heuristic only; exact rounded parameters separately certified'}
    (ROOT/'results/upper_search.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
