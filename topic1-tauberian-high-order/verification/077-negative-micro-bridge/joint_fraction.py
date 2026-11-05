from fractions import Fraction as F
import json
rows=[]
for u in (F(1,64),F(1,256),F(1,1024)):
 g=1/(1+2*u);a=(2*u)**3;k=F(8);Bn=F(3);L0=Bn*k**3/(8*9*g)
 rho_old=g;rho_new=g+2*a*L0*u**6
 mo=a*rho_old**3;mn=a*rho_new**3;e=k*u**3
 dp=(1-e)**2*(1+2*e);dm=(1+e)**2*(1-2*e)
 ratio=lambda d,m:d/(d-m)
 scale=a*(2*u)**9
 det=(ratio(dm,mn)-ratio(dp,mn))/scale;mass=(ratio(dp,mn)-ratio(dp,mo))/scale
 rows.append({'u':str(u),'scalar_positivity_guard':dp>mn>0 and dm>mn and rho_new<1,'exact_det_positive':det>0,'exact_mass_positive':mass>0,'det_approx':float(det),'mass_approx':float(mass),'sum_approx':float(det+mass),'det_numerator_bits':det.numerator.bit_length(),'det_denominator_bits':det.denominator.bit_length()})
print(json.dumps({'scope':'fixed baseline lambda=0,c=0; exact arithmetic; decimal shown for convenience only','rows':rows},indent=2))
