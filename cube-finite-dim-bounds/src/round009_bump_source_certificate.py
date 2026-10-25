"""One exact necessary coefficient bound, fixed c=G(edge)/3, unit ball bump."""
from pathlib import Path
import json,hashlib
import sympy as s
from flint import arb,ctx
z,R=s.symbols('z R',real=True,nonzero=True)
q=3*z*z+2*z+1
F=(R-1)/z+2*R*s.log(s.Abs(z))-R*s.log(q)+R/s.sqrt(2)*s.atan((3*z+1)/s.sqrt(2))
assert s.simplify(s.diff(F,z)-(1-R/q)/z**2)==0
x,y,w=s.symbols('x y w',real=True);bs=[0,-(y+w)/2,(x-w)/2,(x-y)/2]
assert s.expand(sum((bs[i]-bs[j])**2 for i in range(4) for j in range(i+1,4))-(x*x+y*y+w*w))==0
ctx.prec=128;pi=arb.pi();T=4*arb(2).log()-1;r=15/(2*T);root=(3*r-2).sqrt();a=(-1-root)/3;b=(-1+root)/3
I=(r-1)*(1/a-1/b)+2*r*((-a)/b).log()+r/arb(2).sqrt()*(pi-2*(root/arb(2).sqrt()).atan())
c=T/(24*pi**3);mu_upper=4*c*I/pi;lb=32/(9*pi);eps=mu_upper/lb
assert eps<arb(3347)/1000000
candidate=arb(1)/200;source_upper=mu_upper-candidate*lb;assert source_upper<0
out={'status':'rigorous necessary coefficient upper bound and one rejected correction; no Cn bound','fixed_c':'(4 log2-1)/(24 pi^3)','b':'(1-|x|^2)_+^2','point':[1,0,0],'G_Jensen_lower':'5/(16 pi^3 |x|^2)','mu_density_upper_ball':str(mu_upper),'epsilon_necessary_upper_ball':str(eps),'simple_strict_necessary_condition':'epsilon < 3347/1000000','rejected_epsilon':'1/200','corrected_source_density_upper_ball':str(source_upper),'scope':'only fixed c and single unit-ball bump; no cost/g certification, no other c/b/frame exclusion','precision_bits':ctx.prec,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(Path(__file__).resolve().parents[1]/'results/round009_bump_source_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
