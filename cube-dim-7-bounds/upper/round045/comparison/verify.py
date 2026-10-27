"""Exact scalar/full-cost comparison gates, conditional on operator acceptance."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

root=Path(__file__).resolve().parent
main=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')
geo=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds')
soft=Path('/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds')
private=Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds')
L0=2**240
sqrtL0=2**120
params={
 'inverse_display':{'d':F(1,2**30),'sigma':F(1,2**55),'Delta':F(1,2**46),'c':F(1,2**47)},
 'soft_display':{'d':F(1,2**64),'sigma':F(1,2**106),'Delta':F(1,2**96),'c':F(1,2**97)},
 'inverse_domain_endpoint':{'d':F(1,2**24),'sigma':F(1,2**43),'Delta':F(1,2**40),'c':F(1,2**41)},
}
for name,v in params.items():
 assert 0<v['sigma']<F(1,2)
 assert v['c']==v['Delta']/2
 err=F(2,L0**2)+F(2)/(v['sigma']**2*L0**2)+F(16)/(v['sigma']**2*L0*sqrtL0)
 margin=v['c']-2*v['sigma']
 assert 0<err<margin
 v['scalar_margin']=margin;v['E_at_L0']=err
assert params['inverse_display']['sigma']==32*params['inverse_display']['d']**2
assert params['inverse_display']['Delta']==params['inverse_display']['d']/65536
assert params['soft_display']['sigma']==2**22*params['soft_display']['d']**2
assert params['soft_display']['Delta']==params['soft_display']['d']/2**32
assert params['inverse_domain_endpoint']['sigma']==32*params['inverse_domain_endpoint']['d']**2
assert params['inverse_domain_endpoint']['Delta']==params['inverse_domain_endpoint']['d']/65536
E42=F(4,L0)+F(1,2**31*2**80)
Epre=F(4,L0)+F(32,2**120)
assert E42==F(1,2**238)+F(1,2**111)
assert Epre==F(1,2**238)+F(1,2**115)
assert E42<params['soft_display']['c']
assert Epre<params['soft_display']['c']
assert F(4,L0)+params['soft_display']['Delta']<params['inverse_display']['c']
assert F(4,L0)+params['inverse_display']['Delta']<params['inverse_domain_endpoint']['c']
assert L0>=2**108
assert 240<16**2 and 240*F(2,3)>128
# Guard audit at maximum d; all smaller positive fixed d satisfy these algebra gates.
dmax=F(1,2**24)
assert 32*dmax**2<F(1,2)
assert 12*dmax<F(1,64)
assert 64*dmax**2==(dmax/65536)/4
assert 1280*F(1,2**20)==F(5,4096)<F(1,8)
# Exact old certified q values, never decimal displays.
data_path=main/'upper/round034/compute/c34a_results.json'
oldrows=json.loads(data_path.read_text())['results']
qs=[F(row['certificate']['q']) for row in oldrows]
assert len(qs)==12 and max(qs)==F(74275164113,10**12)
assert all(0<q<F(1,8) for q in qs)

paths={
 'candidate_inverse':geo/'upper/round045/soft_inverse_column/PROOF.md',
 'candidate_soft':soft/'upper/round045/layered_source/PROOF.md',
 'old_full_comparison':private/'upper/round042/comparison/PROOF.md',
 'accepted_R042_full':geo/'upper/round042/second_moment_cut/PROOF.md',
 'accepted_R042_threshold':private/'upper/round042/explicit_threshold/THRESHOLD_PROOF.md',
 'old_q_rows':data_path,
}
sources={k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()}
assert sources['candidate_inverse']['sha256']=='add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9'
assert sources['candidate_soft']['sha256']=='22db93b275416b90e9f15a0118a88cd412fee0c2d0b20c4cce1ed32769f259ac'
assert sources['accepted_R042_full']['sha256']=='722764654cd15a8bf7fb18a799c298b54ced2f0188f5554ccae173240ddf90a2'
assert sources['accepted_R042_threshold']['sha256']=='a73bae47dbef7b56b95dcc28bf5dc368828265e91c0de49d03a3dcf518e8c0c5'
receipt={
 'status':'PASS: exact conditional full-cost comparison algebra',
 'numeric_formula_domain':'log(n)>=2^240',
 'operator_status':'conditional on each new source theorem acceptance and its own unknown valid N; no new operator threshold certified',
 'parameters':{k:{kk:str(vv) for kk,vv in v.items()} for k,v in params.items()},
 'old_comparisons':['R042 complete root cost','R033/R034/R035/R036/R040 listed complete old families on valid domains'],
 'unit_core_all_t':'eventual only; old asymptotic threshold not quantified',
 'exact_old_max_q':str(max(qs)),
 'sources':sources,
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('PROOF.md','verify.py')},
 'no_production_search':True,'no_author_source_changes':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'proof_sha256':receipt['source_sha256']['PROOF.md'],
                  'complete_cost_order':'inverse_domain_endpoint < inverse_display < soft_display < valid listed old complete costs',
                  'operator_activation':'conditional; own unknown N retained'},ensure_ascii=False,indent=2))
