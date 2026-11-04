from fractions import Fraction as F
import json
out=[]
for t in (F(1,2),F(1,4),F(1,8)):
 v=t**8
 for label,p in [('threshold',t**19),('strict',t**19*(1+t)),('saturated',v)]:
  A=p**3; H=min(p,A/v**3); R=(p*p*H/v**2)/A
  out.append({'t':str(t),'label':label,'A_gt_v57_8':A>t**57,'A_eq_v57_8':A==t**57,'truncation_saturated':p*p>=v**3,'R_ge_v_minus_1_4':R>=t**-2,'R_over_v_minus_1_4':str(R*t**2)})
for t in (F(1,2),F(1,4)):
 v=t**3;p=t**7;A=p**3;H=min(p,A/v**3);R=(p*p*H/v**2)/A
 out.append({'t':str(t),'label':'v7_example','A_eq_v7':A==v**7,'R_eq_v_minus_1_3':R==1/t,'R_over_1_over_v':str(R*v)})
print(json.dumps(out,indent=2))
