#!/usr/bin/env python3
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,datetime
O=Path(__file__).resolve().parent
j=json.loads((O/'input.json').read_text())
d=F(j['d']);a=F(j['a']);b=F(j['b']);theta=F(j['theta']);h=1-theta
m=F(j['spike_mass']);rho=F(j['background_density']);R=F(j['background_cube_radius']);big=F(j['large_component_mass_threshold']);coeff=F(j['L6_leakage_coefficient'])
r0=a-d;M=rho*(2*R)**6+m
out=[]
def check(name,predicate,**vals):
 assert predicate,name
 out.append({'name':name,'pass':True,**{k:str(v) for k,v in vals.items()}})
check('01_parameter_order',0<d<a<b<1 and r0>0,d=d,a=a,b=b,r0=r0)
check('02_total_mass',M==2080,background_mass=rho*(2*R)**6,spike_mass=m,total_mass=M)
check('03_unique_big',2*big>m,twice_big_threshold=2*big,total_spike=m)
S=rho+big/(64*r0**6)
check('04_small_bound',S<h,S=S,h_minus_S=h-S)
B=F(1,2)/r0**6+F(1,4)+(1-a)/(4*r0)
check('05_big_bound',B<h,B=B,h_minus_B=h-B)
V=(b-a)*(2*b)**5;allow=coeff*M
check('06_leakage',V>allow,slab_volume=V,L6_allowance=allow,strict_gap=V-allow)
check('07_support_margin',b+d==1 and 2*b+d<2,b_plus_d=b+d,two_b_plus_d=2*b+d)
check('08_constants',theta==F(43,15625)/256 and coeff==F(43,15625)/4 and rho<h,theta=theta,coefficient=coeff,h=h)
assert len(out)==8
r={'status':'PASS_ALL_8_EXACT_CHECKS','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'arithmetic':'Python fractions.Fraction; no floats','checks':out,'kernel_calls':0,'union_engine_calls':0,'matlab_calls':0,'spatial_grids':0,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((O/'input.json').read_bytes()).hexdigest()}
(O/'results.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
