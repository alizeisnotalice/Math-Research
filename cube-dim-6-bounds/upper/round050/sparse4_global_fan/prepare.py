#!/usr/bin/env python3
"""Hash/pin/limit preparation only; no active system, kernel or fan construction."""
import sys
sys.dont_write_bytecode=True
import ast,datetime,hashlib,json,math
from pathlib import Path
import flint
HERE=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
a=ROOT/'Math-Research-cube6-theory/cube-dim-6-bounds/upper/round046/sparse4_first_cone_refinement'
f=ROOT/'Math-Research-cube6-literature/cube-dim-6-bounds/audit/round046/sparse4_local_fan_review'
b=HERE.parents[1]/'round044/sparse4_frame_preflight'
old=HERE.parents[1]/'round003/direct_kernel/fan_integer_input.txt'
pins={a/'manifest.json':'e0b943a42323ddb0f9cf2c937141bb252ee17c93e55927fd3de343cba6f74ca7',
      f/'report.md':'07dc3074ba73c184a8eea7f73f4e3a55bbfec55f583554a9f38b8570b448cd7d',
      f/'manifest.json':'6b347c5b9aa900f11893afbd0a7e706fda066c7891496d119eed7226bdd59edd',
      f/'results.json':'6bec898cb72118c8526b4363acd94cc7513aa94a95fad24e0ec0b61ed32eb1ad',
      b/'frame.json':'67a12cd67961baf5ee554c138c82d6b4b627701c66dc691cb6f8b2a823c9be1c',
      old:'7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'}
for p,h in pins.items():assert sha(p)==h,(p,sha(p))
paths=list(pins)
for folder in (a,b):
    m=json.loads((folder/'manifest.json').read_text())
    wh=m.get('files_sha256',m.get('whitelist_sha256'))
    for n,h in wh.items():
        p=folder/n;assert sha(p)==h;pins[p]=h
for folder,names in [(f,['proof_review.md','review.py','input_manifest.json','execution_receipt.json']),
                     (HERE.parents[2]/'round002/fan',[])]:
    for n in names:pins[folder/n]=sha(folder/n)
project=HERE.parents[2]
for rel in ['upper/round002/fan/exact_cells_results.json','upper/round002/fan/verify_exact_fan_results.json',
            'audit/upper002/independent_fan_audit_results.json','upper/round004/primary/integer_fan_bridge.json']:
    p=project/rel;pins[p]=sha(p)
runtime_root=Path(flint.__file__).parent
for p in [Path(sys.executable),Path(flint.__file__)]+sorted((runtime_root/'.dylibs').glob('*.dylib'))+sorted(runtime_root.rglob('*.so')):
    pins[p]=sha(p)
cross=json.loads((b/'results.json').read_text())['per_cone']
planned=[{'parent':x['cone_1_based'],'crossing_forms':x['split_form_count'],
          'active_systems':math.comb(x['split_form_count']+6,5)}for x in cross if x['cone_1_based']!=1]
assert len(planned)==20
total=sum(x['active_systems']for x in planned);assert total<=5000000
ast.parse((HERE/'build_global.py').read_text())
reg={'status':'REGISTERED_READY_NOT_STARTED','registered_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
     'release_status':'EXPLICIT_MAIN_RELEASE_AFTER_F046_FINAL_PASS',
     'release_text':'F046 final PASS gate fulfilled; parent explicitly releases B050 registered remaining20-parent finite batch, only after own registration/method/budget. No first-parent reenumeration or heat kernel.',
     'limits':{'new_parents':20,'active_systems':5000000,'full_cells':200000,'full_simplices':1000000,'actual_seconds':1200},
     'planned_per_parent_active_systems':planned,'planned_total_new_active_systems':total,
     'inherited_parent1_fan':str(a/'local_fan.json'),'inherited_parent1_sha256':sha(a/'local_fan.json'),
     'frame_path':str(b/'frame.json'),'frame_sha256':sha(b/'frame.json'),'old_fan_path':str(old),
     'source_sha256':{str(p):h for p,h in pins.items()},'builder_sha256':sha(HERE/'build_global.py'),
     'method_sha256':sha(HERE/'method.md'),'preparation_source_sha256':sha(__file__),
     'checker_sha256':sha(HERE/'check_saved.py'),
     'saved_verification_policy':'Read full classified ledger, exact determinants and saved rational witnesses, no new solve/vertex search; up to planned_total additional matrix checks, both stages inside actual1200s and combined5000000 matrix budget.',
     'dependencies':{'python':sys.version,'python_flint':flint.__version__,'flint':getattr(flint,'__FLINT_VERSION__',None)},
     'clock_policy':'Actual1200-second clock starts inside sole build batch after gate/hash reads; preparation/wait separate.',
     'output_policy':'Per-parent compressed full active classification ledger, complete vertices/face DAG/simplices, true152rates/determinants/volume; checkpoints preserve completed/uncovered. Final full rational fan only after all21 parents and shared faces pass.',
     'stop_policy':'Cap/time/failure stops once; retain partial ledger/checkpoint; no reset/retry/extension. No CPU or memory artificial shutdown setting.',
     'no_kernel_G_H_K_p_native_T_grid_optimizer_MATLAB':True,'new_C6_bound':False}
assert not(HERE/'registration.json').exists()
(HERE/'registration.json').write_text(json.dumps(reg,indent=2)+'\n')
print(json.dumps({'status':reg['status'],'registration_sha256':sha(HERE/'registration.json'),
                  'builder_sha256':reg['builder_sha256'],'planned_active_total':total,'planned_per_parent':planned,
                  'sources':len(pins),'actual_active_systems_so_far':0},indent=2))
