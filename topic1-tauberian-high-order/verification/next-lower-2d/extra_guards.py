from pathlib import Path
from fractions import Fraction as F
import json
import argparse
parser=argparse.ArgumentParser();parser.add_argument("--input",type=Path,default=Path(__file__).with_name("finite_family.json"));src=parser.parse_args().input
r=json.loads(src.read_text());u=F(r["u"]);e=F(r["e"]);w=F(r["w_fraction"]);rho=F(r["rho_fraction"]);a=4*u*u;M=F(r["M"])
idx=int(r["hstar"]["lower_hex"],16);lo=(w-e)*idx/2**400;hi=(w-e)*(idx+1)/2**400
FF=lambda z:a*(2*rho*z+z*z)
JJ=lambda z:a*(rho*z*z+z**3/3)
Al=2*(e*lo-JJ(hi)/e);D=4*w*w-M;rhs_hi=D-2*Al
Hlo=(2*w-e)*lo-lo*lo/2+JJ(lo)/e
checks={"strong_double_guard_18":FF(w-e)>18*e*e,"slice_rhs_positive":D-2*(2*(e*hi-JJ(lo)/e))>0,"slice_rhs_below_H_hlow":rhs_hi<Hlo,"N_matches":int(r["N"])==u**(-9),"analytic_family_w":w==u/(1+2*u)+(F(128,27)-F(1280,27)*u+F(1033,9)*u*u)*u**6}
assert all(checks.values());print(json.dumps({"checks":checks,"all_pass":True,"scope":"cached-root exact guards; no new root search"}))
