"""Exact comparisons of the already executed public/subpacket certificates."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import math
import json

base = Path(__file__).resolve().parent
def value(name):
    cfg = json.loads((base/"exactinputs"/(name+".json")).read_text())
    result = json.loads((base/"exactinputs"/(name+".check2.json")).read_text())
    mass = sum(F(p[2]) for p in cfg["pts"])
    ratio = F(int(result["ratio_num"]), int(result["ratio_den"]))
    return cfg, mass, ratio, ratio*mass

values = {name:value(name) for name in
          ("frame","frame_line_A","frame_line_B","c2lower_1p70272")}
out = {"values":{}, "fill_increments":{}}
with localcontext() as ctx:
    ctx.prec = 50
    decimal = lambda x: str(Decimal(x.numerator)/Decimal(x.denominator))
    for name, (_, mass, ratio, area) in values.items():
        out["values"][name] = dict(mass=str(mass), area=str(area), ratio=str(ratio),
            ratio_decimal=decimal(ratio), area_decimal=decimal(area))
    for old,new in (("frame","frame_line_A"),("frame","frame_line_B"),
        ("frame_line_A","c2lower_1p70272"),("frame_line_B","c2lower_1p70272"),
        ("frame","c2lower_1p70272")):
        _,m0,r0,A0=values[old]; _,m1,r1,A1=values[new]
        cost=m1-m0; gain=A1-A0
        out["fill_increments"][old+" -> "+new] = dict(added_mass=str(cost),
            added_mass_decimal=decimal(cost), gained_area=str(gain),
            gained_area_decimal=decimal(gain), gain_per_added_mass=str(gain/cost),
            gain_per_added_mass_decimal=decimal(gain/cost),
            strict_ratio_increase=m0*gain>A0*cost)
    cfg,mass,ratio,area=values["c2lower_1p70272"]
    D=2**40
    Px,Py,reach=map(F,(cfg["Px"],cfg["Py"],cfg["reach"]))
    bound=(int(reach/Px)+2)*(int(reach/Py)+2)*mass*D*D/4
    maxk=math.isqrt(bound.numerator//bound.denominator)+1
    rho_max=F(maxk,D)
    q=math.ceil(rho_max/min(Px,Py))+1
    N=10**12
    finite=ratio*F(N-2*q,N)**2
    out["finite_strict_L1"] = dict(D=D, epsilon=str(F(1,2*D)),
        rho_max_rational=str(rho_max), rho_max_decimal=decimal(rho_max),
        q=q,N=N,ratio_lower=str(finite),ratio_lower_decimal=decimal(finite),
        strictly_exceeds_17027237335_over_1e10=finite>F(17027237335,10**10))
    budget=(Px*Py-area)/ratio
    out["fixed_period_fill_budget"] = dict(uncovered_certificate_area=str(Px*Py-area),
        uncovered_certificate_area_decimal=decimal(Px*Py-area),
        necessary_added_mass_upper_bound=str(budget),
        necessary_added_mass_upper_bound_decimal=decimal(budget))
assert all(v["strict_ratio_increase"] for v in out["fill_increments"].values())
assert out["finite_strict_L1"]["strictly_exceeds_17027237335_over_1e10"]
(base/"exact_comparisons.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out["fixed_period_fill_budget"],indent=2))
print("All exact comparison assertions passed.")
