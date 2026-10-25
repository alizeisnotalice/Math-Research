from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time
p=Path(__file__).resolve().parents[1];start=time.monotonic();sources=[]
def read(rel):
 a=p/rel;sources.append((rel,hashlib.sha256(a.read_bytes()).hexdigest()));return json.loads(a.read_text())
C=read('lower_round003/a2_equals_a3_results.json');clo=Q(C['ratio_interval'][0]);den=int(C['root_isolation_denominator']);lo=[Q(int(x),den) for x in C['root_lower_numerators']];hi=[x+Q(1,den) for x in lo]
assert lo[0]**5<2<hi[0]**5 and lo[1]**5<3<hi[1]**5
def bounds(coeff):
 L=H=Q(0)
 for key,val in coeff.items():
  i,j=map(int,key.split(','));c=Q(val);a=lo[0]**i*lo[1]**j;b=hi[0]**i*hi[1]**j
  L+=c*(a if c>=0 else b);H+=c*(b if c>=0 else a)
 return L,H
rows=[]
d=read('lower_opt014/exact_results.json');assert d['C_coefficients']==C['coefficients']
for idx,r in enumerate(d['candidates']):
 keys=set(r['ratio_coefficients'])|set(C['coefficients']);coeff={k:str(Q(r['ratio_coefficients'].get(k,'0'))-Q(C['coefficients'].get(k,'0'))) for k in keys}
 L,H=bounds(coeff);assert H<0;rows.append({'family':'opt014','candidate':idx,'ratio_minus_full_C_interval':[str(L),str(H)]})
for r in read('lower_round014/completion_receipt.json')['rows']:
 H=Q(r['exact_equality_radius_finite_family_ratio_enclosure'][1]);gap=clo-H;assert gap>0;rows.append({'family':'D5_finite_outer','candidate':r['candidate'],'C_lower_minus_outer':str(gap)})
rot=read('lower_round015/rotation_4_3_5_results.json');gap=clo-Q(rot['ratio_interval'][1]);assert gap>0;rows.append({'family':'rotation_finite_outer','C_lower_minus_outer':str(gap)})
for mode in ['fixed','codilation']:
 d=read(f'astra_lower001/exact_{mode}_results.json');gap=clo-Q(d['ratio_exact']);assert gap>0;rows.append({'family':'satellite_rational_inner_only','mode':mode,'C_lower_minus_inner':str(gap)})
derivs=[]
for mode in ['first_tail','tail_tail']:
 d=read(f'lower_alt014/{mode}_outer_theta_results.json');gap=clo-Q(d['ratio']);assert gap>0;rows.append({'family':'parity_theta_outer','mode':mode,'C_lower_minus_outer':str(gap)})
 d=read(f'lower_alt014/{mode}_derivative.json');L,H=bounds(d['right_derivative_coefficients']);assert H<0;derivs.append({'mode':mode,'right_derivative_interval':[str(L),str(H)],'scope':'saved coefficient arithmetic only'})
assert len(rows)==13
out={'status':'13_lower_comparisons_and_2_derivative_signs_passed','rows':rows,'derivatives':derivs,'scope':'no geometry union replay; finite outer, rational inner and local derivative scopes remain distinct','source_hashes':dict(sources),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start}
(p/'round006/check_lower_certificates.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['elapsed_seconds'])
