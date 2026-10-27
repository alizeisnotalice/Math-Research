"""Exact full-cost comparison gates; no old/new parameter or n=7 search."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

root=Path(__file__).resolve().parent
main=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')
geo=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds')
private=Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds')
data_path=main/'upper/round034/compute/c34a_results.json'
data=json.loads(data_path.read_text())
rows=[]
for item in data['results']:
    assert item['status']=='certified'
    q=F(item['certificate']['q']) # exact fraction; never q_approx/raw decimal
    assert 0<q<F(1,8)
    rows.append({'Y':item['Y'],'q':str(q),'required_L_at_least':max(64,item['Y']**2)})
assert len(rows)==12
maximum=max(rows,key=lambda row:F(row['q']))
assert maximum['Y']==256
assert F(maximum['q'])==F(74275164113,10**12)
assert F(data['baseline']['q'])==F(1,20)<F(1,8)
assert 36<41<121 and 5<9
assert F(39,320)<F(1,8)
assert 2<4 # sqrt(2)<2, no displayed decimal used.
L0=2**240
assert all(L0>=row['required_L_at_least'] for row in rows)
assert L0>=2**108
assert 240*F(2,3)>128 and 240*F(2,3)>3
assert 240<16**2
assert F(32,2**40)==F(1,2**35)
assert F(4,2**160)==F(1,2**158)
assert F(1,2**35)+F(1,2**158)<F(1,2**32)

paths={
 'old_parameter_fine_asymptotic':main/'upper/round030/parameter_theory/proof.md',
 'unit_core_all_t':main/'upper/round030/geometry_theory/REPORT.md',
 'unit_core_fixed_t_closure':main/'upper/round031/sigma_theory/REPORT.md',
 'R033':main/'upper/round033/typical_shell/PROOF.md',
 'R034_constant_family':main/'upper/round034/constant_theory/PROOF.md',
 'R034_exact_parameter_rows':data_path,
 'R035':geo/'upper/round035/energy_theory/PROOF.md',
 'R036_explicit_old_guard':main/'upper/round036/threshold_theory/THRESHOLD_PROOF.md',
 'R040_energy':geo/'upper/round040/energy_refinement/PROOF.md',
 'frozen_R042_full':geo/'upper/round042/second_moment_cut/PROOF.md',
 'frozen_R042_moment':private/'upper/round042/second_moment_audit/PROOF.md',
 'frozen_R042_threshold':private/'upper/round042/explicit_threshold/THRESHOLD_PROOF.md',
}
sources={name:{'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
         for name,path in paths.items()}
assert sources['frozen_R042_full']['sha256']=='722764654cd15a8bf7fb18a799c298b54ced2f0188f5554ccae173240ddf90a2'
assert sources['frozen_R042_moment']['sha256']=='425a383952f341f42ff0dee9651924b00fff61e958c2957f355e8cc9fdc72d00'
assert sources['frozen_R042_threshold']['sha256']=='a73bae47dbef7b56b95dcc28bf5dc368828265e91c0de49d03a3dcf518e8c0c5'
receipt={
 'status':'PASS: exact full-cost comparison gates',
 'explicit_comparison_domain':'log(n)>=2^240 intersected with each old theorem valid domain',
 'exact_R034_rows':rows,
 'exact_R034_maximum':maximum,
 'theory_q_upper':'39/320 < 1/8',
 'base_normalized_error_upper':'2^-35 + 2^-158 < 2^-32',
 'tail_proof':'analytic monotonicity when log(L)>3; not finite-sample extrapolation',
 'unit_core_comparison':'eventual only; no explicit start inferred from the fine asymptotic O term',
 'sources':sources,
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                  for name in ('PROOF.md','verify.py')},
 'no_production_search':True,'no_decimal_parameter_reconstruction':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'rows':len(rows),'max_q':maximum,
                  'proof_sha256':receipt['source_sha256']['PROOF.md'],
                  'three_R042_frozen_hashes':'unchanged',
                  'base':'2^-35+2^-158 < 2^-32'},ensure_ascii=False,indent=2))
