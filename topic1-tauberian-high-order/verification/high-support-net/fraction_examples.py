from fractions import Fraction as F
import json
out=[]
for q in (3,4):
 t=F(1,8);v=t**(q-1);s=t**(2*q-1)
 for scale in (F(1,2),F(1),F(2)):
  p=scale*s; A=p**q; P=p**(q-1); H=min(p,A/v**q); price=H*P/v**(q-1)
  out.append({'q':q,'scale':str(scale),'v':str(v),'p_each':str(p),'A':str(A),'price_over_A':str(price/A),'net_sign':(A>price)-(A<price),'p_le_v':p<=v,'saturated':H==p,'boundary_P':P==v**(2*q-1)})
print(json.dumps(out,indent=2))
