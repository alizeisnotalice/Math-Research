"""Small exact audit: no full enumeration or source mutation."""
from fractions import Fraction as F
from pathlib import Path
from decimal import Decimal, localcontext
import hashlib, json, math

base=Path(__file__).resolve().parents[2]
cdir=base/"checks/lower"
src=base/"proofs/lower-baseline.md"
expected="7bc3fa85b86760b4c8193cec1ffee6cd86e242fe4a4ee1bc2d7e0e980561b9ae"
assert hashlib.sha256(src.read_bytes()).hexdigest()==expected
cfg=json.loads((cdir/"exactinputs/c2lower_1p70272.json").read_text())
Px,Py,reach=map(F,(cfg["Px"],cfg["Py"],cfg["reach"]))
pts=[tuple(map(F,p)) for p in cfg["pts"]]
mass=sum(p[2] for p in pts)
D=2**40
nb=(int(reach/Px)+2)*(int(reach/Py)+2)
z=nb*mass*D*D/4
rmax=F(math.isqrt(z.numerator//z.denominator)+1,D)
q=4
assert nb==30 and rmax==F(6915482572195,D)
assert q*min(Px,Py)>rmax
assert rmax*rmax>nb*mass/4
R=F(482167717377842991882283179477574801392641890014291050323,
     283174368139841919615520431313000857449005056000000000000)
N=10**12
finite=F(N-2*q,N)**2*R
assert finite>F(17027237335,10**10)
v=json.loads((cdir/"exact_comparisons.json").read_text())
assert F(v["values"]["c2lower_1p70272"]["mass"])==mass
assert F(v["values"]["c2lower_1p70272"]["ratio"])==R
assert F(v["finite_strict_L1"]["ratio_lower"])==finite
for name,ids in [("frame",range(2)),("frame_line_A",range(26)),
                 ("frame_line_B",list(range(2))+list(range(26,50)))]:
 sub=json.loads((cdir/f"exactinputs/{name}.json").read_text())
 assert [tuple(map(F,p)) for p in sub["pts"]]==[pts[i] for i in ids]
 for key in ("Px","Py","reach","lam"):
  assert F(sub[key])==F(cfg[key])
a=F(v["values"]["c2lower_1p70272"]["area"])
assert a==R*mass and a<Px*Py
assert F(v["fixed_period_fill_budget"]["necessary_added_mass_upper_bound"])==(Px*Py-a)/R
with localcontext() as ctx:
 ctx.prec=45
 dec=lambda x:str(Decimal(x.numerator)/Decimal(x.denominator))
 print("source SHA256:",expected)
 print("positive atoms:",len(pts),"; nb:",nb,"; qminP>rho_max:",True)
 print("rho_max:",str(rmax))
 print("uniform smoothing epsilon:",str(F(1,2*D)))
 print("ratio:",str(R))
 print("ratio decimal:",dec(R))
 print("finite strict L1 bound:",dec(finite))
 print("subpacket membership and fixed periods: True")
 print("portable exact input/output checks: True")
 print("fixed-period additional-mass necessary upper bound:",dec((Px*Py-a)/R))
 print("This audit does not re-enumerate full halos or certify any new tube patch.")

