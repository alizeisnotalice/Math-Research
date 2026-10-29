"""Frozen stream-payload audit plus deterministic one-dimensional replays."""
from pathlib import Path
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile,time
import numpy as np
HERE=Path(__file__).resolve().parent
SKIP={'python','python_version','numpy','seconds','runtime_seconds','runtime_seconds_receipt',
      'superseded_own_review_SHA256'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compare(a,b):
    if isinstance(a,dict):
        assert set(a)-SKIP==set(b)-SKIP,(set(a)-set(b),set(b)-set(a))
        for k in set(a)-SKIP:
            if not k.endswith('_seconds'):compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-10),(a,b)
    else:assert a==b,(a,b)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output');args=parser.parse_args()
    started=time.time();prov=json.loads((HERE/'provenance.json').read_text())
    for name,digest in prov['payload_SHA256'].items():assert sha(HERE/name)==digest,name
    plan=json.loads((HERE/'round040_stream/PLAN.json').read_text())
    assert plan['CI_family_count']==504 and plan['N_per_R']==256
    old=HERE.parent/'038/round038_highdim/source_tables.py'
    runs=[]
    with tempfile.TemporaryDirectory(prefix='cube-round040-')as td:
        root=Path(td)
        for folder in ['round040_source','round040_audit','round040_stream']:shutil.copytree(HERE/folder,root/folder)
        (root/'round038_highdim').mkdir();shutil.copy2(old,root/'round038_highdim/source_tables.py')
        def run(folder,script,output,fields=None):
            expected=json.loads((HERE/folder/output).read_text())
            (root/folder/output).unlink()
            proc=subprocess.run([sys.executable,str(root/folder/script)],cwd=root,
                                env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
            if proc.returncode:raise RuntimeError(proc.stdout+proc.stderr)
            got=json.loads((root/folder/output).read_text())
            if fields:
                for key in fields:compare(expected[key],got[key])
            else:compare(expected,got)
            runs.append({'script':folder+'/'+script,'matches':True})
            return got
        stats=run('round040_source','read_stream_stats.py','STREAM_STATS_AUDIT.json')
        result=run('round040_audit','result_checks.py','result_checks.json')
        run('round040_audit','pilot_checks.py','pilot_checks.json')
        run('round040_audit','source_checks.py','source_checks.json')
        run('round040_source','validate_source.py','SOURCE_VALIDATION.json')
        run('round040_audit','correlation_probe.py','correlation_results.json')
        run('round040_audit','mixed_tilt_probe.py','mixed_tilt_results.json')
    out={'scope':'Frozen hashes, saved sample arithmetic and one-dimensional diagnostics only',
         'new_tail_samples':0,'full_cache_builds':0,'payload_hash_count':len(prov['payload_SHA256']),
         'python':sys.version,'numpy':np.__version__,'runs':runs,
         'CIs_recomputed':stats['one_sided_CIs_independently_recomputed'],
         'max_CI_difference':stats['max_CI_absolute_error'],
         'same_h_comparisons':result['same_h_comparisons'],
         'positive_same_h_MC_count':sum(c['same_h_MC_positive_count']for c in result['cases']),
         'positive_same_h_CI_endpoint_count':sum(c['same_h_endpoint_positive_count']for c in result['cases']),
         'original_peak_RSS_bytes':result['max_RSS_bytes'],'full_FP_certificate':False,
         'Cn_log_gain_proved':False,'seconds':time.time()-started}
    if args.output:
        p=Path(args.output)
        if p.exists():raise RuntimeError('Refusing replacement')
        p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
