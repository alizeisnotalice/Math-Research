"""Recheck all archived bounds/CI arithmetic without rerunning source counts."""
from pathlib import Path
import argparse,json,math,os,subprocess,sys,tempfile,time
import numpy as np
from replay_support import prepare,sha,HERE
SKIP={'python','python_version','numpy','seconds','runtime_seconds','elapsed_seconds'}
def compare(a,b):
    if isinstance(a,dict):
        assert set(a)-SKIP==set(b)-SKIP,(set(a)-set(b),set(b)-set(a))
        for k in set(a)-SKIP:
            if not k.endswith('_seconds'):compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-9),(a,b)
    else:assert a==b,(a,b)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output');args=parser.parse_args()
    started=time.time();prov=json.loads((HERE/'provenance.json').read_text())
    for name,digest in prov['payload_SHA256'].items():assert sha(HERE/name)==digest,name
    runs=[]
    with tempfile.TemporaryDirectory(prefix='cube-round041-')as td:
        root=prepare(td)
        def run(folder,script,output,expected_folder=None):
            original=HERE/(expected_folder or folder)/output
            expected=json.loads(original.read_text())
            target=root/folder/output;target.unlink()
            p=subprocess.run([sys.executable,str(root/folder/script)],cwd=root,
                env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
            if p.returncode:raise RuntimeError(p.stdout+p.stderr)
            got=json.loads(target.read_text());compare(expected,got)
            runs.append({'script':folder+'/'+script,'matches':True})
            return got
        run('round040_source','variance_probe.py','CANDIDATE_VARIANCE_RESULTS.json','variance')
        run('round041_mcd','probe.py','RESULTS.json')
        run('round041_mcd','uniform_probe.py','uniform_results.json')
        for script,output in [('postprocess_checks.py','postprocess_checks.json'),
                              ('mcd_checks.py','mcd_checks.json'),
                              ('uniform_lip_checks.py','uniform_lip_checks.json'),
                              ('method_uniform_ci_checks.py','method_uniform_ci_checks.json'),
                              ('replay_checks.py','replay_checks.json')]:
            run('round041_audit',script,output)
        # Re-run the final arithmetic using the recorded per-scale inputs.
        (root/'round041_post/FINAL_PLAN.json').unlink()
        run('round041_post','aggregate.py','FINAL_RESULTS.json')
        assert sha(root/'round041_post/FINAL_PLAN.json')==sha(HERE/'round041_post/FINAL_PLAN.json')
        # Numeric replay was compared above. Downstream frozen audits bind the original
        # result bytes, including its runtime metadata; restore those bytes in this
        # temporary view rather than compare a new runtime-dependent SHA as identical.
        (root/'round041_post/FINAL_RESULTS.json').write_bytes(
            (HERE/'round041_post/FINAL_RESULTS.json').read_bytes())
        ordinary=run('round041_mcd','offline_upper_check.py','OFFLINE_UPPER_AUDIT.json')
        uniform=run('round041_mcd','uniform_offline_check.py','uniform_offline_check.json')
        ci=run('round041_audit','final_ci_checks.py','final_ci_checks.json')
    out={'scope':'Stored arithmetic and one-dimensional fixtures, not a model replay',
         'model_replay_reexecuted':False,'new_independent_samples':0,'source_cache_builds':0,
         'shared_replay_NPZs_checked_against_actual_replay_metadata_SHA':True,
         'identity_checks_use_shared_archived_bytes':True,
         'payload_hash_count':len(prov['payload_SHA256']),'python':sys.version,'numpy':np.__version__,
         'runs':runs,'CI_count':ci['CI_count'],'max_CI_difference':ci['max_Decimal_CI_difference'],
         'full_FP_certificate':False,'Cn_log_gain_proved':False,'seconds':time.time()-started}
    if args.output:
        p=Path(args.output)
        if p.exists():raise RuntimeError('Refusing replacement')
        p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
