import json,sys,hashlib
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round133/fixed_corner_envelope_diagnostic')
source=A/'diagnostic.json';d=json.loads(source.read_text());out={}
for t in ['1/63500','1/65000']:
 strong=[];rectangle_only=[];bad_total=ordered_total=0
 for row in d['results']:
  bad=[c for c in row['corners']if Q(c['fixed_robust_family_upper'])<Q(t)];bad_total+=len(bad);ordered=[]
  for c in bad:
   point=[Q(1)]+list(map(Q,c['face_point']))
   if all(a>=b for a,b in zip(point,point[1:]))and point[-1]>=0:ordered.append(c)
  ordered_total+=len(ordered)
  if ordered:
   witness=min(ordered,key=lambda c:(Q(c['fixed_robust_family_upper']),c['corner_bits']))
   strong.append({'node_id':row['node_id'],'corner_bits':witness['corner_bits'],'ordered_body_point':['1']+witness['face_point'],'upper_exact':witness['fixed_robust_family_upper'],'target_minus_upper_exact':str(Q(t)-Q(witness['fixed_robust_family_upper']))})
  elif bad:rectangle_only.append(row['node_id'])
 out[t]={'strong_ordered_domain_obstruction_box_count':len(strong),'only_full_rectangle_obstruction_box_count':len(rectangle_only),'total_rectangle_obstruction_corners':bad_total,'ordered_domain_obstruction_corners':ordered_total,'only_full_rectangle_obstruction_node_ids':rectangle_only,'ordered_witnesses':strong}
result={'status':'PASS_ORDERED_DOMAIN_FIXED_FAMILY_SUPPLEMENT','diagnostic_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'targets':out,'scope':'each original box fixed selected centers and independent endpoint rectangle robust formula; not union325 expanded candidates; not upper bound on G','new_exp_native_checker_producer_calls':0}
(P/'ordered_domain_classification.json').write_text(json.dumps(result,indent=2)+'\n');print({t:{k:v for k,v in x.items()if k not in ['ordered_witnesses','only_full_rectangle_obstruction_node_ids']}for t,x in out.items()})
