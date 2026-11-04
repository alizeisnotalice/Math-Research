from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from pathlib import Path
import json,time

def dec(x,mode):
 with localcontext() as c:c.prec=45;c.rounding=mode;return str(Decimal(x.numerator)/Decimal(x.denominator))
def out(lo,hi):return [dec(lo,ROUND_FLOOR),dec(hi,ROUND_CEILING)]
start=time.perf_counter();u=F(1,1024);k=F(16,3)-F(64,3)*u;lam=F(128,27)-F(1280,27)*u+F(1024,9)*u**2;a=4*u*u;d0=u/(1+2*u);w=d0+lam*u**6+u**8;rho=1-2*w;e=k*u**3;hwidth=w-e;det=1-e**2;m0=(2*d0)**2;M=a*rho**2;mout=(2*w)**2;D=mout-M;N=2**90
FF=lambda h:a*(2*rho*h+h*h);JJ=lambda h:a*(rho*h*h+h**3/3)
guards={'e_positive':e>0,'e_lt_w_over_2':e<w/2,'rho_gt_2w':rho>2*w,'Fprime_gt_e':2*a*rho>e,'double_budget':FF(hwidth)>2*e**2,'det_gt_outer':det>mout,'M_gt_half_outer':2*M>mout}
assert all(guards.values())
lo=F(0);hi=hwidth
for _ in range(400):
 mid=(lo+hi)/2
 if FF(mid)<e*e:lo=mid
 else:hi=mid
Al=2*(e*lo-JJ(hi)/e);Au=2*(e*hi-JJ(lo)/e);budgetlow=D-2*Au;budgethigh=D-2*Al
radius=rho+w;defectupper=(M+2*e**2)/(2*radius)**2;strictmargin=u*u-defectupper;infnorm=1/(1-e);finite=F(2*N-5,2*N)**2*M/(det-M)-m0/(1-m0)
G=(M-det*m0)/((det-M)*(1-m0)*a*(2*u)**6);limit=F(4,27)
guards.update(root_bracket=FF(lo)<e*e<=FF(hi),h_lt_e=hi<e,A_positive=Al>0,budget_strict=budgetlow>0,center_strict=strictmargin>0,finite_cube_contained=infnorm*radius<N,finite_ratio_gt_Ln=finite>0)
idx=lo*2**400/hwidth;assert idx.denominator==1
result={'u':str(u),'kappa':str(k),'lambda_rational':str(lam),'family':'kappa=16/3-64u/3; lambda=128/27-1280u/27+1024u^2/9','w':'delta*+lambda_rational*u^6+u^8','e':str(e),'w_fraction':str(w),'rho_fraction':str(rho),'det_fraction':str(det),'M':str(M),'N':str(N),'geometric_input':{'T':'[[1,e],[e,1]]','G':'outer square excluding non-strict witness set','selected_A':'G intersect {x1<c} where c is chosen so |A|=M','E_N':'T([-N,N]^2) minus union over z in Z^2 of (Tz+A)'},'hstar':{'scale':'w-e','pow2':400,'lower_hex':hex(idx.numerator),'upper_hex':hex(idx.numerator+1)},'A_bounds':out(Al,Au),'D_minus_2A':out(budgetlow,budgethigh),'center':{'coordinates':['0','0'],'radius':str(radius),'defect_upper':str(defectupper),'strict_margin':str(strictmargin),'mean_bound':'actual mean >= 1-defect_upper >1-u^2; geometric containment is separately assumed'},'finite_halo_minus_Ln':{'exact_fraction':str(finite),'positive':finite>0},'G':out(G,G),'formal_limit':str(limit),'guards':guards,'all_pass':all(guards.values()),'seconds':time.perf_counter()-start,'scope':'finite conditional construction using proved geometry; exact rational family; not a certified global optimizer'}
Path('finite_family.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_pass':result['all_pass'],'guards':guards,'G':result['G'],'budget':result['D_minus_2A'],'seconds':result['seconds']}))
