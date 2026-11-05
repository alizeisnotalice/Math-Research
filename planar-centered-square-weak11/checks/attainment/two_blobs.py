from fractions import Fraction as F
from pathlib import Path
import os,sys,json
sep=F(19,10); eps=F(1,100)
s=F(49,50); rho=F(69,50)
rs=F(99,100); rd=F(139,100)
single=2*(2*s)**2-(2*s-sep)*(2*s)
extra=(2*rho-sep)*2*(rho-s)
total=single+extra; mass=F(8)
assert s+eps==rs and rho+eps==rd
assert F(4)>4*rs**2 and F(8)>4*rd**2
assert total/mass==F(10317,10000)>1
out={'python':sys.version,'pid':os.getpid(),'epsilon':str(eps),'single_union_area':str(single),'double_extra_area':str(extra),'receiver_area_lower':str(total),'input_mass':str(mass),'ratio_lower':str(total/mass),'strict_single':True,'strict_double':True}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(out,ensure_ascii=False,indent=2))
