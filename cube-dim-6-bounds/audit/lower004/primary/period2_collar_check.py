# Exact geometric check for the lower-bound period-two witness family.
from fractions import Fraction as Q
from pathlib import Path
import json,hashlib
u=Q(1122462,1000000);v=Q(1200936,1000000);a=Q(26,25)
assert u**6<2 and v**6<3
sp=[1+u**4/2,(u**3+u**4*v**3)/4]+[(1+u)/2]*4
assert all(x>a for x in sp)
assert a**6>Q(6,5)
out={'status':'PASS_exact_uniform_radius_to_spacing','delta_range':'0<=delta<=1/5','all_spacings_strictly_above':'26/25','r_max_strictly_below':'2*(26/25)','reason':'mass<=4^6*(1+delta); r=mass^(1/6)/2<2 min_spacing; half-cell source indices contained in {-1,0,1,2}, reflection gives {-1,0,1,2,3}; coarse collar2 sufficient','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'])
