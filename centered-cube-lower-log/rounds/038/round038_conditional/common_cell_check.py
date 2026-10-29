"""Deterministic common-D source mass calculation; no high-n tail simulation."""
from pathlib import Path
from fractions import Fraction as F
import math
import json
import hashlib
import sys
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    out = HERE/"COMMON_CELL.json"
    if out.exists():
        raise RuntimeError("Refusing replacement.")
    k, v, p = 255, 8, .5
    n = 116*v*(k//3)**2
    m = math.ceil(math.sqrt(v))
    freq = np.arange(1,m+1,dtype=float)
    b = 4*(729/29)*((np.sin(np.pi*freq*10/27)/(np.pi*freq))**4
                    +3*(np.sin(np.pi*freq/9)/(np.pi*freq))**4)
    c = np.sqrt(2*b/3)
    edges = np.arange(k+1)/k
    # rho_bg=1-p+(p/K) sum pi*m*c_m*cos(2*pi*m*y).
    primitive = (p/(2*k))*(np.sin(2*np.pi*edges[:,None]*freq)@c)
    cell = np.full(k,1/k)+np.diff(primitive)
    pchi = float(cell[:k//2].sum())
    analytic_pchi = k//2/k + p/(2*k)*float(np.sin(2*np.pi*(k//2/k)*freq)@c)
    assert abs(analytic_pchi-pchi)<1e-15
    bg_lower = float(F(1,2)-F(20,6*k)*sum((F(1,j)for j in range(1,m+1)),F(0)))
    assert bg_lower>0 and cell.min()>0
    data = {"scope":"common cell law only, no high-dimensional tail experiment",
        "K":k,"v":v,"n":n,"M":m,"p":p,"w":4,"delta":1/(4*k*v),
        "coarse_plus_cells":k//2,"common_background_lower_bound":bg_lower,
        "cell_mass_sum":float(cell.sum()),"cell_mass_min":float(cell.min()),
        "cell_mass_max":float(cell.max()),"actual_pchi":pchi,
        "pchi_minus_half":pchi-.5,
        "binomial_count_mean":n*pchi,"half_binomial_count_mean":n/2,
        "half_mean_shift_in_actual_sd":n*(pchi-.5)/math.sqrt(n*pchi*(1-pchi)),
        "v_ge_16logn":v>=16*math.log(n),"B_ge_v":4>=v,
        "original_matched_lower_bound_claimed":False,
        "cell_masses":cell.tolist(),"c_m":c.tolist(),
        "source_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python":sys.executable,"numpy":np.__version__,
        "floating_point_integral_not_interval_certified":True}
    out.write_text(json.dumps(data,separators=(",",":")))
    print(json.dumps({key:value for key,value in data.items()
                      if key not in ["cell_masses","c_m"]},indent=2))


if __name__=="__main__":
    main()
