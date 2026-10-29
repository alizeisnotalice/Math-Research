from pathlib import Path
from fractions import Fraction as F
import json,hashlib
H=Path(__file__).resolve().parent;src=H.parent/'round048_banks/COVARIANCE.json';co=json.loads(src.read_text());p=F(co['parameters']['p']);n=int(co['parameters']['n']);K=co['parameters']['K'];rows=[]
for r in co['scales']:
 if r['integer_scale']:continue
 L=F(r['L']);a=p/(K*L);kappa=a*(2+F(40,9));diff=F(r['exact_H_diagonal'])-F(r['exact_H_cross_all_six_pairs']);E=n*a*a*diff/(4*(1+kappa));m=E.numerator//E.denominator
 cap=F(12,2**m)
 rows.append(dict(L=str(L),kappa=str(kappa),E_exact=str(E),E_float=float(E),floor_E=m,affinity_upper_power_of_two='2^(-'+str(m)+')',fourbank_uniform_anchor_epsilon='1/16',overlap_probability_rational_upper=str(cap),overlap_probability_float_reference=float(cap)))
out=dict(source_SHA256=hashlib.sha256(src.read_bytes()).hexdigest(),rows=rows,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Exact Fraction evaluation of candidate true-source Hellinger/anchor-overlap bounds, to be proved separately; no Gaussian assumption or C_n gain.')
(H/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
