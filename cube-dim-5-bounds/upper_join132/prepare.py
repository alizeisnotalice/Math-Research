#!/usr/bin/env python3
"""Read-only source binding and finite132 plan, zero H and zero budget eval."""
import sys
sys.dont_write_bytecode=True
import ast,gzip,hashlib,json,os
from pathlib import Path

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
SCRATCH=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_join132')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p);return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def write(p,x):Path(p).write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')

def main():
    assert not (HERE/'registration.json').exists(),'inspect existing132;no duplicate setup'
    names=['upper_shorttime076/results.json','upper_shorttime076/report.md','upper_shorttime076/manifest.json',
      'upper_singular057/results.json','upper_singular057/report.md','upper_singular057/manifest.json',
      'upper_path112/numeric_hierarchy.json.gz','upper_path112/nodes.json.gz','upper_path112/node_certificates.json.gz','upper_path112/report.md','upper_path112/manifest.json','upper_path112/production_v2.py',
      'audit_path121/whole_results.json','audit_path121/report.md','audit_path121/manifest.json',
      'upper_alpha1_029/report.md','upper_alpha1_029/manifest.json',
      'upper_adaptive103/mean_floor_certificate_scope.json','upper_adaptive103/report.md','audit_adaptive113/floor_results.json','audit_adaptive113/report.md']
    locks={str(MAIN/n):sha(MAIN/n) for n in names}
    old=HERE.parent/'audit_g123_sol';m=read(old/'manifest.json');assert sha(old/'manifest.json')=='a961e0fd00c580cb495e951835aab2fa7d15408e754ff39670dc1cfcb9704c17'
    protected={str(old/'manifest.json'):sha(old/'manifest.json'),**{str(old/k):v for k,v in m['files'].items()}}
    for p,h in protected.items():assert sha(p)==h
    s=read(MAIN/'upper_shorttime076/results.json');e=next(x for x in s['cases'] if x['epsilon']=='1/512')
    h=read(MAIN/'upper_path112/numeric_hierarchy.json.gz');a=read(MAIN/'audit_path121/whole_results.json')
    assert a['status']=='pass' and a['final']['final_budget']['UH']==h['budget']['UH']
    assert h['budget']['UH']=='32424400796227698610199181381141800612723929/730750818665451459101842416358141509827966271488'
    assert h['weighted_volume']=='1' and len(h['frontier'])==94580 and len(h['nodes'])==97972
    leaf_volume_types=sorted({(h['nodes'][k]['level'],h['nodes'][k]['weight']) for k in h['frontier']})
    floor=read(MAIN/'upper_adaptive103/mean_floor_certificate_scope.json');p3=read(MAIN/'upper_singular057/results.json')['constants']['constants']['P3']
    inp=dict(epsilon='1/512',q='5/4',UH=h['budget']['UH'],E_U=e['E_norm']['upper'],
        C3_U=s['derived_constants']['C3']['upper'],C1_U=s['derived_constants']['C1']['upper'],P3_U=p3['upper_dyadic'],
        pi_L=s['derived_constants']['pi_outward_384_dyadic']['lower'],pi_U=s['derived_constants']['pi_outward_384_dyadic']['upper'],
        existing_I_cap='297621/6250000000',existing_exact_I_upper=a['final']['final_budget']['I_upper'],g='21/200000',
        threshold_L=floor['threshold_T_g105_lower'],threshold_U=floor['threshold_T_g105_upper'],
        constants_role='saved upper endpoints define a new upper-certificate object;not actual H/R values or lower bounds')
    write(HERE/'inputs.json',inp)
    for p in HERE.glob('*.py'):ast.parse(p.read_text())
    code={str(p):sha(p) for p in HERE.glob('*.py')}
    reg=dict(task='upper_join132_prepare',status='awaiting_primary_math_and_finite_budget_review',pid=os.getpid(),input_locks=locks,protected123_locks=protected,code_locks=code,
        inputs_sha256=sha(HERE/'inputs.json'),finite_cases=['zero_cross_floor','T_1over32','T_5over128','T_3over64','T_1over16','analytic_stationary'],
        rational_Ts=['1/32','5/128','3/64','1/16'],maximum_cases=6,workers=0,single_process=True,wall_cap=None,
        formula='UH+(5/4)*B(T)^(1/4)*log(1+256eps^2)/(512*pi_L)+E_U^(5/4)',
        B_formula='C3_U*(eps^-2-T^-2)/2+C1_U*log(T/eps)+P5_U/(4*T^4)',P5_formula='P3_U/(64*pi_L)',
        stationary_formula='T^2=2P5_U/(C3_U+sqrt(C3_U^2+4C1_U*P5_U));clip to eps boundary if needed',
        root_bits=384,log_bits=512,log_terms=192,new_H=0,new_R_pointwise=0,new_gradient=0,new_numeric_log=0,new_budget_cases=0,
        old_baseline_recomputed=False,actual_budget_authorized=False,global_C5_claim=False,
        acceptance='exact32 gate and new-I-upper vs accepted old exact I upper; if stationary certificate floor exceeds threshold,exclude only this fixed-envelope/fixed-E/g105 certificate family',
        reuse_scope='consume112/121 full IH and076 complete R norm;no old computation rerun;029 same generalL1/all-positive-scale centered-common-radius model',
        future_cell_metadata=dict(frontier_count=94580,level_weight_types=leaf_volume_types,volume_formula='weight/2^(4level)',new_mass_caps_evaluated=0,new_quarter_roots_evaluated=0))
    write(HERE/'registration.json',reg)
    print(json.dumps(dict(registration_sha256=sha(HERE/'registration.json'),input_files=len(locks),protected123_files=len(protected),finite_cases=6,new_H=0,new_log=0,new_budget_cases=0)))

if __name__=='__main__':main()
