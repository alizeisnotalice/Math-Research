#!/usr/bin/env python3
"""Lock proof/code and source-time metadata only;no real H/q fields read."""
from pathlib import Path
from fractions import Fraction as F
import ast,hashlib,json,os
H=Path(__file__).resolve().parent
B=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert not (H/'registration.json').exists()
    names=['primary/band_curvature172.md','upper_alpha1_029/report.md','audit_frame075/report.md','audit_curvature129/report.md',
      'upper_curvature153/report.md','upper_curvature153/penalty.py','upper_retarget155/core155.py','round018/primary_curvature155_proof.md',
      'upper_path112/scalar_registration.json','primary/finite_source_transfer.md']
    locks={str(B/n):sha(B/n) for n in names}
    times=json.loads((B/'upper_path112/scalar_registration.json').read_bytes())['time_endpoints'];assert len(times)==73
    selected=[t for t in times if F(t)>=F(1,64)];assert len(selected)==49 and selected[0]=='1/64' and selected[-1]=='1'
    for p in H.glob('*.py'):ast.parse(p.read_text())
    protected={str(H.parent/'audit_time160/manifest.json'):'6dc7bb171b762e6914a62f3a8b4da497344ebc0b14eb903c638b83f7da7545d1',
      str(H.parent/'upper_integral167/manifest.json'):'72c70c8a44edc3d4dbffca9342336101ab46ce70ae2ba0beb710a3729fbf7790',
      str(H.parent/'upper_integral167/execution_manifest.json'):'c36fb6329b5318c1b7c84871fc39feb840b8feed8edd49c2826d69298bb24f97'}
    for p,h in protected.items():assert sha(p)==h
    r={'task':'audit_bandcurvature174_theory_interface_only','pid':os.getpid(),'source_locks':locks,'code_locks':{str(p):sha(p) for p in H.glob('*.py')},
      'protected_old_manifests':protected,'binding_template_sha256':sha(H/'binding_template.json'),'all73_timekeys':times,'selected49_timekeys':selected,
      'new_symbolic_case':'positive_Kupper_over_Hlower_residual_only','symbolic_cases':1,'old153_weight_or_symbolic_baseline_replay':False,
      'real_H_fields_read':0,'real_q_bandproducts_corners_evaluations':0,'actual_registered_sources':0,'production173_status':'UNBOUND','wall_cap':None}
    with (H/'registration.json').open('x') as f:json.dump(r,f,sort_keys=True,indent=2);f.write('\n')
    print(json.dumps({'pid':os.getpid(),'registration_sha256':sha(H/'registration.json'),'real_H_q_corner_calls':0,'source_files':len(locks)}))
if __name__=='__main__':main()
