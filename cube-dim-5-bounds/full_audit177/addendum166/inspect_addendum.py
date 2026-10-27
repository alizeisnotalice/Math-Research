"""Source/hash and new saved-field checks only. No gain, DFS or root replay."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,os,time
HERE=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T')
A=BASE/'Math-Research-cube5-audit136-sol/cube-dim-5-bounds/audit_combine169'
P=BASE/'Math-Research-cube-dim-5/cube-dim-5-bounds/round019'
def read(p):return json.loads(p.read_bytes())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()
def put(n,x):
    with (HERE/n).open('x') as f:json.dump(x,f,ensure_ascii=False,indent=2);f.write('\n')
def main():
    start=time.monotonic();assert not (HERE/'source_manifest.json').exists()
    original=read(HERE.parent/'frozen_manifest.json')
    for n,row in original['files'].items():assert sha(HERE.parent/n)==row['sha256']
    assert sha(HERE.parent/'frozen_manifest.json')=='7138ac40a84c53b4e4a739d89784705e8c987086600b46491eb586c7e29e2912'
    locks={}
    def lock(path,expected=None):
        p=Path(path);h=sha(p)
        if expected is not None:assert h==expected,str(p)
        if str(p) in locks:assert locks[str(p)]['sha256']==h
        locks[str(p)]={'sha256':h,'bytes':p.stat().st_size}
    plan,b,c=map(read,[A/'plan.json',A/'source_binding.json',A/'candidate_binding.json'])
    lock(A/'source_binding.json',plan['source_binding_sha256'])
    lock(A/'candidate_binding.json',plan['candidate_binding_sha256'])
    for n,h in plan['code_locks'].items():lock(A/n,h)
    for obj in [b['source_locks'],c['candidate_locks']]:
        for row in obj.values():lock(row['path'],row['sha256'])
    for row in b['shards']+[b['tree']]:lock(row['path'],row['sha256'])
    for p in A.rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts:lock(p)
    reg=read(P/'upper166_registration.json')
    lock(P/'accept_upper166.py',reg['code_sha256'])
    lock(P/'upper166_registration.json')
    for row in reg['sources'].values():lock(row['path'],row['sha256'])
    for n in ['primary_upper166_acceptance.json','primary_proof166.md']:lock(P/n)
    r=read(A/'run001/results.json');e=read(A/'execution_receipt.json')
    assert sha(A/'run001/results.json')=='75ac91b7c743dad1e705aa4b60dfce39b2ea8d033eb6becd011afb74e0dabe59'
    assert r['pid']==21868 and r['worker_exitcodes']==[0]*15
    assert (r['new_gain_cases'],r['nodes'],r['roots'])==(97971,97972,52360)
    assert e['observed_master_exitcode']==0 and e['actual_worker_exitcodes']==r['worker_exitcodes']
    assert r['source_hashes_before_after_equal'] and e['source_before_after']['all_hashes_equal']
    result=read(P/'primary_upper166_acceptance.json')
    old=read(Path(reg['sources']['accepted156']['path']))
    prod=read(Path(reg['sources']['producer166']['path']))
    assert result['pid']==22326 and result['status']=='primary_accepts_complete_I166_and_C5_upper'
    assert result['new_UH']==prod['new_UH']==r['new_UH']
    assert result['qualified_g']==old['qualified_g']
    assert result['same_complete_J139_upper']==old['same_complete_J139_upper']
    assert Q(result['complete_I_upper'])==Q(result['new_UH'])+Q(result['same_complete_J139_upper'])
    assert Q(old['complete_I_upper'])-Q(result['complete_I_upper'])==Q(result['exact_I_gain_vs156'])>0
    assert Q(result['best_upper']['fourth_power'])<Q(old['best_upper']['fourth_power'])
    newlo,newhi=map(Q,result['best_upper']['interval']);oldlo,oldhi=map(Q,old['best_upper']['interval'])
    assert 0<newlo<newhi<oldlo<oldhi
    assert Q(result['simple_rational_upper'])==Q(28936079,1000000)>=newhi
    assert result['registration_sha256']==sha(P/'upper166_registration.json')
    assert result['sources']==reg['sources'] and result['new_bound_root_isolations']==1
    assert result['general_L1_all_positive_common_radii'] and result['H128_floor_domain']=='complete_F_only'
    assert result['G_floor_domain']=='complete_Q_via_G_homogeneity' and not result['lower_changed']
    lock(HERE/'inspect_addendum.py')
    put('source_manifest.json',{'status':'actual_sources_and_dependency_locks_match',
      'original177_frozen_manifest_unchanged':True,'files':locks,
      'hashes_establish_identity_only':True,'169_actual_new_math_not_replayed':True})
    put('inspection_result.json',{'pid':os.getpid(),'seconds':time.monotonic()-start,
      'source_files':len(locks),'source_bytes':sum(v['bytes'] for v in locks.values()),
      'source_manifest_sha256':sha(HERE/'source_manifest.json'),
      'all_original177_files_unchanged':True,'169_result_matches_requested_sha':True,
      '169_source_locks_checked':len(b['source_locks']),'169_candidate_locks_checked':len(c['candidate_locks']),
      '169_independent_shards_checked':len(b['shards']),'169_tree_sha_checked':True,
      'new_saved_field_comparisons_pass':True,'new_H_R_mean_gain_slope_DFS_root_calls':0,
      'primary_new_U166_display_only':result['best_upper']['display_only'],
      'simple_rational_upper':result['simple_rational_upper'],
      'status':'first_review_addendum_166_169_and_primary_qualification_pass_not_second_audit'})
    print(json.dumps({'files':len(locks),'manifest_sha256':sha(HERE/'source_manifest.json'),
      'result_sha256':sha(HERE/'inspection_result.json'),'seconds':time.monotonic()-start}))
if __name__=='__main__':main()
