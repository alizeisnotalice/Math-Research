"""Exact G(a,1,1) log formula and independent conditional-density residues."""
from pathlib import Path
import sympy as s,json
OUT=Path(__file__).resolve().parent
a,u,t=s.symbols('a u t',positive=True)
p=(a-1)**2;q=(a+1)**2;S=a*a+2;T=a**4+14*a*a+9;C=2*(a*a-1)**2
N=20*u**3+24*S*u**2+4*T*u+4*C
R=N/((u+4)*(u+p)**2*(u+q)**2)
A=s.factor(N.subs(u,-4)/((p-4)**2*(q-4)**2))
Cp=s.factor(N.subs(u,-p)/((4-p)*(q-p)**2));Cq=s.factor(N.subs(u,-q)/((4-q)*(p-q)**2))
Bp=s.factor(s.diff(N/((u+4)*(u+q)**2),u).subs(u,-p));Bq=s.factor(s.diff(N/((u+4)*(u+p)**2),u).subs(u,-q))
assert s.factor(R-(A/(u+4)+Bp/(u+p)+Cp/(u+p)**2+Bq/(u+q)+Cq/(u+q)**2))==0
assert s.factor(A+Bp+Bq)==0
# M(a)=(a-1)^2 integral R; use log(p/4),log(q/4) for exact common-scale cancellation.
mp=s.factor(-(a-1)**2*Bp);mq=s.factor(-(a-1)**2*Bq);mr=s.factor((a-1)**2*(Cp/p+Cq/q))
# Independent conditional integral: b=(0,-1,d,d), d=(a-1)/2.
d=(a-1)/2;I=s.I;z0=I*t;z1=-1+I*t;zd=d+I*t
res0=1/(2*I*t*(1+2*I*t)*((-d)*(-d+2*I*t))**2)
res1=1/(2*I*t*(1-2*I*t)*((-1-d)*(-1-d+2*I*t))**2)
H=1/(d*(d+2*I*t)*(d+1)*(d+1+2*I*t)*(2*I*t)**2)
resd=-H*(2*zd/(d*(d+2*I*t))+2*(zd+1)/((d+1)*(d+1+2*I*t))+2/(2*I*t))
pt_res=s.factor(t**4/(4*s.pi**4)*(2*s.pi*I)*(res0+res1+resd))
pt_formula=s.factor((4*t)*R.subs(u,16*t*t)/(4*s.pi**3))
assert s.factor(pt_res-pt_formula)==0
rec={'status':'PASS exact curve identities','curve':'x=(a,1,1), a>0, a!=1; here a ranges [1/4,5/2] excluding1','R_u':'N3/((u+4)(u+(a-1)^2)^2(u+(a+1)^2)^2)','partial_fractions':{'A4':str(A),'Bp':str(Bp),'Cp':str(Cp),'Bq':str(Bq),'Cq':str(Cq)},'M_formula':'(a-1)^2 * 32*pi^3*G(a,1,1) = mp*log((a-1)^2/4)+mq*log((a+1)^2/4)+mr','mp':str(mp),'mq':str(mq),'mr':str(mr),'conditional_density_residue_identity':'PASS independent simple poles it,-1+it and double pole (a-1)/2+it','scope':'one curve for fixed e1 directional source only; no full-frame cost/source optimization'}
(OUT/'cb10_curve_formula.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
