from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,math
r=Path(__file__).resolve().parents[3];src=r/'upper/round002/fan/exact_cells_results.json';dst=r/'upper/round003/direct_kernel/fan_integer_input.txt'
cells=json.loads(src.read_text())['cells'];a=iter(map(int,dst.read_text().split()));assert next(a)==len(cells)
for c in cells:
 det=next(a);scales=[]
 for ray,rate in zip(c['rays_a0_equals_one'],c['rates']):
  q=next(a);v=[next(a) for _ in range(6)];x=list(map(Q,ray));scale=Q(v[0])/x[0];assert scale>0
  assert all(Q(z)==scale*y for z,y in zip(v,x));assert Q(q)==scale*Q(rate);scales.append(scale)
 assert Q(det)==Q(c['determinant'])*math.prod(scales)
assert next(a,None) is None
out={'status':'PASS_exact_ray_rate_determinant_scaling_bridge','cells':len(cells),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'integer_input_sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
