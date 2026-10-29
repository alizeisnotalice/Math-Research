"""Replay deterministic round39 evidence in temporary copies; no tail MC."""
from pathlib import Path
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile,time
import numpy as np

HERE=Path(__file__).resolve().parent
SKIP={'seconds','python','python_version','numpy','source_input_directory'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),(a.keys(),b.keys())
        for k in a:
            if k not in SKIP:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-11),(a,b)
    else:assert a==b,(a,b)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true');parser.add_argument('--output');args=parser.parse_args()
    started=time.time();prov=json.loads((HERE/'provenance.json').read_text())
    for name,digest in prov['payload_SHA256'].items():assert sha(HERE/name)==digest,name
    channel=HERE/'round039_channel';posterior=HERE/'round039_posterior';coupled=HERE/'round039_coupled'
    reg=json.loads((channel/'registration.json').read_text())
    assert reg['script_SHA256']==sha(channel/'probe.py')
    source=HERE.parent/'038/round038_highdim/source_tables.py'
    assert sha(source)==reg['source_SHA256']
    data=json.loads((posterior/'RESULTS.json').read_text())
    assert sha(posterior/'probe.py')==data['source_SHA256']
    note=json.loads((posterior/'EXECUTION_NOTE.json').read_text())
    assert note['initial_source_SHA256']==sha(posterior/'probe_initial_linear.py')
    assert note['completed_source_SHA256']==sha(posterior/'probe.py')
    creg=json.loads((coupled/'registration.json').read_text())
    assert creg['script_SHA256']==sha(coupled/'probe.py')
    assert json.loads((coupled/'results.json').read_text())['plan_SHA256']==sha(coupled/'registration.json')
    old_order=json.loads((channel/'gaussian_order_registration.json').read_text())
    assert old_order['prior_results_SHA256']==sha(channel/'results.json')
    out={'scope':'Exact and deterministic numeric replays only; no tail MC or large cache',
         'payload_hash_count':len(prov['payload_SHA256']),'python':sys.version,'numpy':np.__version__,'full':args.full,'runs':[]}
    with tempfile.TemporaryDirectory(prefix='cube-round039-')as td:
        root=Path(td)
        for name in ['round039_channel','round039_posterior','round039_coupled']:shutil.copytree(HERE/name,root/name)
        (root/'round038_highdim').mkdir();shutil.copy2(source,root/'round038_highdim/source_tables.py')
        def run(folder,script,output,fields=None):
            original=json.loads((HERE/folder/output).read_text())
            (root/folder/output).unlink()
            result=subprocess.run([sys.executable,str(root/folder/script)],cwd=root,
                                  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
            if result.returncode:raise RuntimeError(result.stdout+result.stderr)
            actual=json.loads((root/folder/output).read_text())
            if fields:
                for key in fields:compare(original[key],actual[key])
            else:compare(original,actual)
            out['runs'].append({'script':folder+'/'+script,'replay_matches':True})
            return actual
        run('round039_posterior','probe.py','RESULTS.json')
        run('round039_channel','trace_checks.py','trace_results.json')
        run('round039_channel','multitable_checks.py','multitable_results.json')
        run('round039_channel','gap_checks.py','gap_results.json')
        run('round039_coupled','posterior_audit.py','posterior_audit.json')
        run('round039_coupled','probe.py','results.json')
        if args.full:
            import scipy
            out['scipy']=scipy.__version__
            run('round039_channel','probe.py','results.json')
            # The new source result includes fresh elapsed time; the order
            # replay correctly hashes that new result instead of the archive.
            run('round039_channel','gaussian_order_probe.py','gaussian_order_results.json',
                ['rows','three_candidate_vs_full_bin_max_abs_error','refinement_half_step','qualification'])
    # Independent finite consequence checks from archived Gaussian data.
    order=json.loads((channel/'gaussian_order_results.json').read_text())
    lower_checks=0
    for row in order['rows']:
        tau,R,a,EH=row['tau'],row['R'],row['a_R'],row['E_H'];s=math.sqrt(1-tau*tau)
        upper=3+2*tau*a/(s*math.sqrt(2*math.pi))
        assert EH<=upper+1e-9
        for eps in [.2,.4,.6,.8]:
            c,cp=1-eps,1-eps/2
            phi=math.exp(-.5*(cp*a)**2)/math.sqrt(2*math.pi)
            if R*(cp-c)*a*phi>1:
                d=1/(R*phi)
                lower=2*c*tau*a/(s*math.sqrt(2*math.pi))*math.exp(-tau*tau*d*d/(2*s*s))
                assert EH+1e-6>=lower
                lower_checks+=1
    out.update({'Gaussian_upper_checks':len(order['rows']),'Gaussian_applicable_finite_lower_checks':lower_checks,
                'full_FP_certificate':False,'Cn_log_gain_proved':False,'seconds':time.time()-started})
    if args.output:
        p=Path(args.output)
        if p.exists():raise RuntimeError('Refusing replacement')
        p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
