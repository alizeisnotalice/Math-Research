"""Read full certificates; independent interval arithmetic, no union rerun."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time
p=Path(__file__).resolve().parents[1]; started=time.monotonic()
cpath=p/'lower_round003/a2_equals_a3_results.json'; C=json.loads(cpath.read_text()); opt=json.loads((p/'lower_opt013/results.json').read_text())
assert C['coefficients']==opt['C_coefficients']
den=int(opt['root_denominator']); lo=[Q(int(x),den) for x in opt['root_lower_numerators']];hi=[x+Q(1,den) for x in lo]
assert lo[0]**5<2<hi[0]**5 and lo[1]**5<3<hi[1]**5

def bounds(coeff):
 L=H=Q(0)
 for key,value in coeff.items():
  i,j=map(int,key.split(','));c=Q(value);a=lo[0]**i*lo[1]**j;b=hi[0]**i*hi[1]**j
  L+=c*(a if c>=0 else b);H+=c*(b if c>=0 else a)
 return L,H
rows=[]
for r in opt['candidates']:
 keys=set(C['coefficients'])|set(r['ratio_coefficients']);diff={k:str(Q(r['ratio_coefficients'].get(k,'0'))-Q(C['coefficients'].get(k,'0'))) for k in keys}
 assert all(Q(diff.get(k,'0'))==Q(r['difference_from_exact_C_coefficients'].get(k,'0')) for k in set(diff)|set(r['difference_from_exact_C_coefficients']))
 a,b=bounds(diff);assert b<0
 rows.append({'direction':r['direction'],'exact_difference_interval':[str(a),str(b)]})
clo=Q(C['ratio_interval'][0]); alt=[]
for name in ['weighted_first_tail','weighted_tail_tail','deleted_first_tail','deleted_tail_tail']:
 d=json.loads((p/f'lower_alt013/{name}_outer_results.json').read_text());ratio=Q(d['ratio']); gap=clo-ratio;assert gap>0
 alt.append({'candidate':name,'C_lower_minus_family_outer':str(gap)})
out={'status':'all_12_finite_candidate_comparisons_strictly_below_full_C','source_C_sha256':hashlib.sha256(cpath.read_bytes()).hexdigest(),'opt':rows,'alt':alt,'scope':'independent certificate comparison only; geometry and union rely on separately archived certificates','elapsed_seconds':time.monotonic()-started}
(p/'round005/compare_new_lower.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['elapsed_seconds'])
