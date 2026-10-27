"""Freeze only accepted175 audit inputs; hashes are identities, not proof."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
HERE=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
AUD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-audit136-sol/cube-dim-5-bounds/audit_secant175')
def read(p):return json.loads(Path(p).read_bytes())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()
def put(name,d):
    with (HERE/name).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n')
def main():
    assert not (HERE/'registration.json').exists();files={}
    def bind(path,expected=None):
        p=Path(path);h=sha(p)
        if expected is not None:assert h==expected,str(p)
        files[str(p)]={'sha256':h,'bytes':p.stat().st_size}
    r=read(ROOT/'round020/upper175_registration.json');main=read(ROOT/'round020/primary_upper175_acceptance.json')
    assert main['source_registration_sha256']==sha(ROOT/'round020/upper175_registration.json')
    bind(ROOT/'round020/accept_upper175.py',r['code_sha256'])
    for row in r['sources'].values():bind(row['path'],row['sha256'])
    p,b,c=map(read,[AUD/'plan.json',AUD/'source_binding.json',AUD/'candidate_binding.json'])
    bind(AUD/'source_binding.json',p['source_binding_sha256']);bind(AUD/'candidate_binding.json',p['candidate_binding_sha256'])
    for f,h in p['code_locks'].items():bind(AUD/f,h)
    for locks in [b['source_locks'],c['candidate_locks']]:
        for row in locks.values():
            assert 'upper_band173' not in row['path'];bind(row['path'],row['sha256'])
    for row in b['shards']+[b['tree']]:bind(row['path'],row['sha256'])
    for f in AUD.rglob('*'):
        if f.is_file() and '__pycache__' not in f.parts:bind(f)
    for f in ['primary_upper175_acceptance.json','upper175_registration.json','primary_proof175.md',
       'matlab175_actual.json','matlab175_inputs.json','matlab175_registration.json','matlab175_result.json',
       'verify_upper175_matlab.m','secant_precision.md','secant_precision_registration.json','secant_precision_result.json']:
        bind(ROOT/'round020'/f)
    for f in ['primary/finite_source_transfer.md','primary/positive_face_singular_trace.md','primary/quotient_shorttime.md',
      'round019/primary_upper166_acceptance.json','round019/primary_proof166.md','round019/final_upper_review.json',
      'full_audit177/report.md','full_audit177/addendum166/report.md','second_audit178/report.md',
      'second_audit178/evidence178.tar.gz','second_audit178/evidence178_manifest.json',
      'strict_g180/portable/evidence/report.md',
      'strict_g180/portable/evidence/run001/result.json','strict_g180/portable/README.md']:
        bind(ROOT/f)
    # Main's result is archived, not a missing mathematical result. Lock both
    # its repository archive and original actual bytes, without regenerating it.
    bind(Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-secondaudit178-sol/cube-dim-5-bounds/second_audit178/run001/results.json'),
         'b22867e46c1d431b218e309227215852a3add59bad27a6ca3a77b783e702ce6b')
    bind(HERE/'check_final_upper.py');bind(HERE/'freeze_inputs.py')
    put('input_manifest.json',{'status':'frozen_real_accepted175_sources','files':files,
      'static175_source_locks':len(b['source_locks']),'actual_candidate_locks':len(c['candidate_locks']),
      'shards':len(b['shards']),'173_candidate_excluded':True,'hashes_are_identity_not_enclosure_truth':True})
    names={'accepted175':ROOT/'round020/primary_upper175_acceptance.json','accepted166':ROOT/'round019/primary_upper166_acceptance.json',
      'producer175':Path(r['sources']['producer175']['path']),'audit175':AUD/'run001/results.json',
      'matlab_result':ROOT/'round020/matlab175_result.json','matlab_inputs':ROOT/'round020/matlab175_inputs.json'}
    put('registration.json',{'task':'final_bounds020_upper175_only_while_research_paused','utc':datetime.now(timezone.utc).isoformat(),
      'pid':os.getpid(),'maximum_math_runs':1,'scope':'audit accepted175 final scalar and inheritance;no new model/search/173',
      'method':'integer_cross_products_no_Fraction_SymPy_or_root_algorithm','code_sha256':sha(HERE/'check_final_upper.py'),
      'input_manifest':{'path':str(HERE/'input_manifest.json'),'sha256':sha(HERE/'input_manifest.json')},
      'named_inputs':{k:str(v) for k,v in names.items()},'new_node_root_gain_DFS_H_calls':0,
      'planned_checks':'cleared I/g fourth power identity,savedroot4 inequalities,exact new-old gains,portable strict bound,all32 comparisons,MATLAB sourcebinding'})
    print(json.dumps({'files':len(files),'registration_sha256':sha(HERE/'registration.json'),
      'manifest_sha256':sha(HERE/'input_manifest.json')}))
if __name__=='__main__':main()
