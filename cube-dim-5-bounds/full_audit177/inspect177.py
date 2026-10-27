"""Read/hash frozen sources and compare new166 rational fields; no old replay."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,os,subprocess,time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
BASE=Path('/Users/zhengzhihao/Desktop/T')
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_bytes())
def save(n,d):
    with (HERE/n).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n')

def main():
    start=time.monotonic()
    assert not (HERE/'source_manifest.json').exists()
    local_dirs=['primary','round018','upper_alpha1_029','upper_singular057','audit_frame075',
      'upper_shorttime076','upper_shift091','upper_range092','upper_path112','audit_path121',
      'upper_join132','upper_rearrangement133','upper_facepeak134','upper_increment139',
      'audit_curvature129','upper_curvature153','upper_retarget155','audit_retarget157',
      'upper_integral150','upper_integral156','audit_integral161','upper_time159','audit_time160',
      'upper_integral167','upper_global151','audit_final147']
    files={p for name in local_dirs for p in (REPO/'cube-dim-5-bounds'/name).rglob('*')
           if p.is_file() and '__pycache__' not in p.parts}
    files.update([REPO/'AGENTS.md',REPO/'cube-dim-5-bounds/rounds/018.md',HERE/'inspect177.py'])
    ext_dirs=[BASE/'Math-Research-cube5-curvature128-sol/cube-dim-5-bounds'/s
              for s in ['upper_combine166','upper_slope164','upper_moment162']]
    ext_dirs += [BASE/'Math-Research-cube5-audit136-sol/cube-dim-5-bounds/audit_moment163',
                 BASE/'Math-Research-cube5-auditcurvature129-sol/cube-dim-5-bounds/audit_slope165']
    # Code, theory, registration and final qualification; progress/chunks are bound
    # through original manifests, not copied or replayed here.
    for root in ext_dirs:
        files.update(p for p in root.iterdir() if p.is_file() and p.suffix in ('.py','.json','.md','.diff'))
    done=BASE/'cube-finite-dim-scratch/upper_combine166/run001/completion.json'
    files.add(done)
    records=[{'path':str(p),'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(files)]
    old=read(REPO/'cube-dim-5-bounds/round018/primary_integral156_acceptance.json')
    new=read(done)
    u156=Q(old['new_UH']);u166=Q(new['new_UH'])
    u162=Q(new['baseline162_UH']);u164=Q(new['baseline164_UH'])
    # This is a new166-vs-existing-field comparison, not an old H/mean/DFS/root replay.
    assert 0<u166<u164<u162<u156
    assert u162-u166==Q(new['exact_gain_vs162'])
    assert u164-u166==Q(new['exact_gain_vs164'])
    g=Q(old['qualified_g']);assert g>0
    save('source_manifest.json',{'task':'full_audit177_first_frozen_review',
      'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
      'files':records,'hashes_establish_identity_only':True,
      'old_H_interval_truth':'inherited accepted original rigorous arithmetic; not regenerated',
      '169_result_included':False,'173_175_included':False})
    save('inspection_result.json',{'pid':os.getpid(),'seconds':time.monotonic()-start,
      'source_files':len(records),'source_bytes':sum(x['bytes'] for x in records),
      'new166_exact_comparison':'0<UH166<UH164<UH162<UH156',
      'gain166_vs156':str(u156-u166),'gain166_vs162':str(u162-u166),
      'gain166_vs164':str(u164-u166),
      'new_H_R_J_gradient_log_mean_chord_corner_root_or_DFS_calls':0,
      'new166_comparison_only':True,'old_successful_baseline_replays':0,
      'formal_completed_target':'U156','166_numerical_qualification':'pending169_and_primary',
      'qualified_g':str(g),'source_manifest_sha256':digest(HERE/'source_manifest.json')})
    print(json.dumps({'files':len(records),'source_manifest_sha256':digest(HERE/'source_manifest.json'),
                      'result_sha256':digest(HERE/'inspection_result.json'),'seconds':time.monotonic()-start}))
if __name__=='__main__':main()
