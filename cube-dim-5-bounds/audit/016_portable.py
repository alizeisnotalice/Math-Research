#!/usr/bin/env python3
"""Portable adapter for frozen audit 016; default is a cheap integrity check.

Usage:
  python3 audit/016_portable.py --project-root /path/to/cube-dim-5-bounds
  python3 audit/016_portable.py --project-root /path/to/cube-dim-5-bounds \
      --output-dir /tmp/audit016-replay --replay

--replay alone opts into the previously verified full 120-second computation.
Otherwise only imports, locked hashes, and the stored counterexample are
checked. All output goes to --output-dir or a newly created temporary folder.
The frozen source, result, and producer inputs are never written.

Required replay files: this adapter, 016_staggered_independent_audit.py,
016_staggered_independent_results.json, 010_independent_new_lower.py,
002_independent_lower_audit.py, and the three project-relative locked inputs
listed below. The separate 016 report is optional for execution.
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
import hashlib, json, tempfile

HERE = Path(__file__).resolve().parent
FROZEN = HERE/'016_staggered_independent_audit.py'
RESULT = HERE/'016_staggered_independent_results.json'
LOCKED = {
    '016_staggered_independent_results.json': '1ad913f408e9d7adf918e370b974f55d84806a51cecd9274873fca5858ba3932',
    '010_independent_new_lower.py': '264e5a639b92bb64ba7254dc0b7fb59f96843d0207cfd640512f0901631d7dd8',
    '002_independent_lower_audit.py': '8a5ece75606b1fdb478657955ca22d33a71c426ed0099d74257a46a5bc225cf3',
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=HERE.parent)
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--replay',action='store_true',help='run full exact geometry/union after the integrity checks')
    args=parser.parse_args()
    project=args.project_root.resolve()
    out=args.output_dir.resolve() if args.output_dir else Path(tempfile.mkdtemp(prefix='cube5-audit016-'))
    out.mkdir(parents=True,exist_ok=True)
    if out==HERE:
        raise ValueError('--output-dir must differ from the frozen audit directory')
    baseline=json.loads(RESULT.read_text())
    protected={FROZEN:baseline['source_sha256'],**{HERE/p:h for p,h in LOCKED.items()}}
    inputs={
        project/'lower_round010/two_row_staggered.py':baseline['producer_source_sha256'],
        project/'lower_round010/two_row_tau_1_8_results.json':baseline['producer_result_sha256'],
        project/'lower_round003/a2_equals_a3_results.json':baseline['complete_C_source_sha256'],
    }
    protected.update(inputs)
    for path,expected in protected.items():
        if sha(path)!=expected: raise ValueError(f'locked input changed: {path}')
    spec=spec_from_file_location('portable_frozen016',FROZEN)
    audit=module_from_spec(spec); spec.loader.exec_module(audit)
    audit.BASE=project; audit.PROD=project/'lower_round010'; audit.HERE=out
    assert audit.CORE_PATH==HERE/'010_independent_new_lower.py'
    Alg=audit.Alg
    def decode(data):
        return Alg(coeff=[Q(data.get(f'{i//5},{i%5}','0')) for i in range(25)])
    ex=baseline['counterexample']
    assert ex['retained_tag']==[0,0,[1,4,1,1,1]]
    point=tuple(map(decode,ex['point_coefficients']))
    box=tuple(map(decode,ex['box_endpoint_coefficients']))
    r=decode(ex['radius_coefficients']); s=audit.s; h=audit.h; a1=audit.a1
    assert r==audit.root(4)/2
    assert point==(r/2,s,s/4,s/4,s/4)
    assert box==(Alg(0),r,Alg(0),2*s,Alg(0),s/2,Alg(0),s/2,Alg(0),s/2)
    assert all(box[2*i]<point[i]<box[2*i+1] for i in range(5))
    assert 3*s-point[1]>r and 3*s>2*r
    assert s>r and 3*s/4>r
    assert -r/2<h<3*r/2 and h-a1 < -r/2 and h+a1>3*r/2
    assert ex['entire_periodic_measure_mass_at_this_radius']==1
    assert ex['actual_average_at_this_radius']=='1/4'
    if args.replay:
        audit.main()
    for path,expected in protected.items():
        assert sha(path)==expected, f'protected file changed: {path}'
    receipt={
        'status':'PASS_portable_016_import_hash_counterexample_check',
        'project_root':str(project),'output_directory':str(out),
        'full_volume_replayed':args.replay,
        'frozen_source_and_result_unchanged':True,
        'locked_project_inputs':{str(p.relative_to(project)):h for p,h in inputs.items()},
        'counterexample_strict_box_inclusion':True,
        'counterexample_entire_periodic_cube_mass':1,
        'counterexample_average':'1/4',
        'frozen_result_sha256':LOCKED[RESULT.name],
        'adapter_sha256':sha(Path(__file__)),
    }
    (out/'016_portable_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__': main()
