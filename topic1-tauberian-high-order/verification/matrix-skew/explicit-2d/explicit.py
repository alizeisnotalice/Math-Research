from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from pathlib import Path
import json,time,argparse

def dec(x,mode):
 with localcontext() as c:c.prec=50;c.rounding=mode;return str(Decimal(x.numerator)/Decimal(x.denominator))
def bounds(lo,hi):return [dec(lo,ROUND_FLOOR),dec(hi,ROUND_CEILING)]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();start=time.perf_counter()
 old=json.loads(args.input.read_text());row=next(x['attempts'][-1] for x in old['cases'] if x['n']==2);assert row['status']=='certified';u=F(row['u']);k=F(row['kappa']);lam=F(row['lambda']);n=2;v=2*u;a=v**2;d0=u/(1+2*u);w=d0+lam*u**6;rho=1-2*w;e=k*u**3;det=1-e**2;M=a*rho**2;mout=(2*w)**2;D=mout-M;K=a;N=int(row['N']);hroot=row['roots']['hstar'];assert hroot['scale']=='w-e' and hroot['pow2']==400
 hidx=int(hroot['lower_hex'],16);assert int(hroot['upper_hex'],16)==hidx+1;hl=(w-e)*hidx/2**400;hh=(w-e)*(hidx+1)/2**400
 FF=lambda z:K*(2*rho*z+z*z)
 JJ=lambda z:K*(rho*z*z+z**3/3)
 H=lambda z:(2*w-e)*z-z*z/2+JJ(z)/e
 C=lambda z:JJ(z)/e-z*z/2
 Al=2*(e*hl-JJ(hh)/e);Au=2*(e*hh-JJ(hl)/e)
 rhslo=D-2*Au;rhshi=D-2*Al
 def bisect(target):
  assert 0<target<H(hl);lo=F(0);hi=hl
  for _ in range(400):
   mid=(lo+hi)/2
   if H(mid)<target:lo=mid
   else:hi=mid
  return lo,hi
 zl,_=bisect(rhslo);_,zh=bisect(rhshi);cl=w-zh;ch=w-zl;Cl=C(zl);Cu=C(zh)
 nbrlo=2*e**2-2*Au-Cu;nbrhi=2*e**2-2*Al-Cl;r=rho+w;deflo=(M+nbrlo)/(2*r)**2;defhi=(M+nbrhi)/(2*r)**2;avlo=1-defhi;avhi=1-deflo;marginlo=u**2-defhi;marginhi=u**2-deflo
 oldupper=F(row['center']['defect_upper']);Nrecord=int(row['N']);inv=F(row['center']['inverse_operator_norm'])
 checks={'prior_root_bracket':FF(hl)<e**2<=FF(hh),'Fprime0_gt_e':2*K*rho>e,'hstar_lt_e':hh<e,'e_lt_w':e<w,'rhs_positive':rhslo>0,'rhs_below_H_hlow':rhshi<H(hl),'z_bracket':0<zl<=zh<hl,'z_target_bracket':H(zl)<=rhslo and H(zh)>=rhshi,'C_increasing_guard':2*K*rho>e,'C_positive':Cl>0,'c_positive':0<cl<=ch<w,'A_positive':Al>0,'neighbor_positive':nbrlo>=0,'new_defect_below_prior_bound':defhi<oldupper,'actual_average_above_threshold':marginlo>0,'N_matches':N==Nrecord,'witness_inside_finite_box':inv*r<N,'M_gt_half_outer':2*M>mout}
 idxlo=zl*2**400/hl;idxhi=zh*2**400/hl;assert idxlo.denominator==idxhi.denominator==1
 result={'n':2,'u':str(u),'kappa':str(k),'lambda':str(lam),'N':str(N),'source_hstar':hroot,'z_bracket':{'scale':'hstar_lower','denominator_power2':400,'lower_hex':hex(idxlo.numerator),'upper_hex':hex(idxhi.numerator)},'A':bounds(Al,Au),'z':bounds(zl,zh),'c':bounds(cl,ch),'C':bounds(Cl,Cu),'neighbor_mass':bounds(nbrlo,nbrhi),'defect_fraction':bounds(deflo,defhi),'average':bounds(avlo,avhi),'strict_margin':bounds(marginlo,marginhi),'prior_defect_upper':bounds(oldupper,oldupper),'improvement_over_prior_average_lower_bound':bounds(oldupper-defhi,oldupper-deflo),'checks':checks,'all_pass':all(checks.values()),'seconds':time.perf_counter()-start,'scope':'conditional arithmetic validation of supplied exact geometry; no independent proof of bad-set/slab formulas'}
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_pass':result['all_pass'],'checks':checks,'z':result['z'],'c':result['c'],'margin':result['strict_margin'],'seconds':result['seconds']}))
if __name__=='__main__':main()
