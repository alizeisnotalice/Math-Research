"""Exact 96-cone, rate and spectral-volume certificate; no cost search."""
from pathlib import Path
from fractions import Fraction as F
import runpy,json
import sympy as s
OUT=Path(__file__).resolve().parent;ns=runpy.run_path(str(OUT/'cb11_frame_green.py'));tau=s.symbols('tau',nonnegative=True)
rows=[]
for cc in ns['CONES']:
    det=s.det(s.Matrix.hstack(*(s.Matrix(r) for r in cc)));assert abs(det)==1
    rates=[tuple(ns['rate_coeff'](r)) for r in cc]
    rows.append({'rays':list(map(list,cc)),'det':int(det),'rates_body_face_coefficients':list(map(list,rates))})
A=s.Matrix([1,0,0]);B=s.Matrix([1,1,0]);C=s.Matrix([2,1,1]);D=s.Matrix([1,1,1]);a,b,c=s.symbols('a b c',real=True)
assert s.Matrix.hstack(A,B,C).inv()*s.Matrix([a,b,c])==s.Matrix([a-b-c,b-c,c])
assert s.Matrix.hstack(B,C,D).inv()*s.Matrix([a,b,c])==s.Matrix([b-c,a-b,b+c-a])
AA=4+4*tau;BB=4+6*tau;CC=8+10*tau;DD=6+6*tau
vol=s.factor(48*(1/(AA*BB*CC)+1/(BB*CC*DD))/6)
peak=s.factor(6*vol/8);assert s.factor(peak-5/(8*(1+tau)*(2+3*tau)*(4+5*tau)))==0
rays={r for cc in ns['CONES'] for r in cc};assert len(rays)==50
rec={'status':'PASS exact 96-cone algebra','cones':rows,'signed_permutation_chambers':48,'base_cones':2,'distinct_rays':50,'spectral_unit_ball_volume':str(vol),'p0_times_pi3_t3':str(peak),'tail_G_T_infinity_bound':'p0_times_pi3_t3/(2*pi^3*T^2)','ABC_inverse':['a-b-c','b-c','c'],'BCD_inverse':['b-c','a-b','b+c-a'],'coverage':'explicit nonnegative inverse coordinates cover each sorted positive chamber; signed permutations partition full space except boundaries'}
(OUT/'cb11_cones.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps({k:v for k,v in rec.items() if k!='cones'},indent=2))
