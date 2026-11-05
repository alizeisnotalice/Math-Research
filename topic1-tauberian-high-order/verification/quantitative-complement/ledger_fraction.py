from fractions import Fraction as F
import json
a=F(1,8);b=F(1,4)
e=[8+a,8+b-a,9-a-3*b,10-5*a,11-2*a,F(17,2)-3*a,9-b,9-a]
f=e.copy();f[5]=F(19,2)-3*a
out={'alpha':str(a),'beta':str(b),'base_exponents':[str(x) for x in e],'improved_exponents':[str(x) for x in f],'both_min':'65/8'==str(min(e))==str(min(f)),'delta':str(min(e)-8),'upper_bound_proof':'d<=alpha, d<=beta-alpha imply alpha>=d,beta>=2d; d<=1-alpha-3beta then 8d<=1'}
print(json.dumps(out,indent=2))
