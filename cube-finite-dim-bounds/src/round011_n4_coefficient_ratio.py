"""One fixed exact coefficient-ratio pass; no source or ledger mutation."""
from pathlib import Path
from fractions import Fraction as F
import json, hashlib
from flint import fmpq_mpoly_ctx, fmpq, arb, ctx
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'upper/u11_minimum_global_results.json'
SHA='43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SHA
C=fmpq_mpoly_ctx.get(('a','b','c','d'),ordering='lex');X=C.gens()
R=[[(1,0,0,0),(1,1,0,0),(2,1,1,0),(3,1,1,1)],
   [(2,1,1,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)],
   [(2,1,1,0),(1,1,0,0),(1,1,1,1),(3,1,1,1)]]
def poly(v):return C.from_dict({tuple(map(int,k.split(','))):fmpq(a) for k,a in v.items()})
rows=[]
for case,rays in zip(json.loads(SOURCE.read_text())['physical_cases'],R):
    sub=[sum((rays[j][i]*X[j] for j in range(4)),C.constant(0)) for i in range(4)]
    A=(X[0]**3*poly(case['numerator_dict'])).compose(*sub).to_dict()
    D=poly(case['denominator_dict']).compose(*sub).to_dict()
    assert D and all(c>=0 for c in D.values()) and all(c>=0 for c in A.values())
    ratios=[(F(str(A.get(e,0)))/F(str(c)),e) for e,c in D.items() if c>0]
    best,e=min(ratios)
    assert all(F(str(A.get(k,0)))-best*F(str(c))>=0 for k,c in D.items())
    rows.append({'case':case['case'],'coefficient_ratio_lower':str(best),'attaining_power':list(map(int,e)),'D_terms':len(D)})
bound=min(F(r['coefficient_ratio_lower']) for r in rows)
ctx.prec=256;g=(arb(bound.numerator)/bound.denominator)/arb.pi()**3
result={'status':'exact coefficient-ratio pass complete; source/continuity independent audit remains required',
        'source_sha256':SHA,'cases':rows,'a_cubed_pi_cubed_G_lower':str(bound),'g_interval':str(g),
        'strictly_above_old_g_by_22_over_7':bound>F(21,25000)*F(22,7)**3,
        'matches_known_point_17_over_640':bound==F(17,640),
        'new_weak_type_upper_registered':False,
        'execution_note':'First pass arithmetic succeeded but JSON fmpz exponent serialization failed; int conversion repaired, same fixed pass rerun.'}
(ROOT/'results/round011_n4_coefficient_ratio.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
