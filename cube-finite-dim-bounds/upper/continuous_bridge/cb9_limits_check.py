"""Exact normalization/degenerate anchors at already registered spatial classes."""
from pathlib import Path
import sympy as s,json
OUT=Path(__file__).resolve().parent
u=s.symbols('u',positive=True);r,t=s.symbols('r t',positive=True)
X,Y,Z=s.symbols('X Y Z',nonnegative=True);S=X+Y+Z;S2=X*Y+Y*Z+Z*X
C=X*(Y-Z)**2+Y*(Z-X)**2+Z*(X-Y)**2
assert s.expand(C-(S*S2-9*X*Y*Z))==0
axis=4*(5*u+1)/(u+1)**3;edge=4*(5*u+2)/((u+1)**2*(u+4));line=(20*u**2+72*r*r*u+96*r**4)/(u*(u+4*r*r)**3)
axisI=s.integrate(axis,(u,0,s.oo));edgeI=s.integrate(edge,(u,0,s.oo));assert axisI==12 and s.simplify(edgeI-4*(4*s.log(2)-1))==0
line_pt=s.factor(4*t*line.subs(u,16*t*t)/(4*s.pi**3));line_limit=s.limit(t*line_pt,t,0,dir='+');assert line_limit==3/(32*s.pi**3*r*r)
V,Q,z=s.symbols('V Q z',positive=True);space=s.integrate((z*z+Q)**-4,(z,-s.oo,s.oo));tm=s.integrate(t**4*(t*t+V)**-s.Rational(7,2),(t,0,s.oo));assert space==5*s.pi/(16*Q**s.Rational(7,2)) and tm==1/(5*V)
rec={'status':'PASS exact anchors and Feynman normalizations','SOS_identity':'S1*S2-9XYZ = X(Y-Z)^2+Y(Z-X)^2+Z(X-Y)^2','axis_R_u':str(axis),'axis_integral_R':str(axisI),'G_axis_exact':'3/(8*pi^3)','edge_R_u':str(edge),'edge_partial_fractions':str(s.apart(edge,u)),'edge_integral_R':str(edgeI),'G_edge_exact':'(4*log(2)-1)/(8*pi^3)','body_line_R_u':str(line),'body_line_pt':str(line_pt),'limit_t_times_pt_on_bodyline':str(line_limit),'Feynman_space_integral':str(space),'Feynman_time_integral':str(tm),'Feynman_G_prefactor':'3/(32*pi^3)','simplex_volume':'1/6, coordinate Lebesgue measure, not normalized probability','scope':'no additional spatial points beyond registered degenerate classes; no g minimum or cost claim'}
(OUT/'cb9_normalization.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
