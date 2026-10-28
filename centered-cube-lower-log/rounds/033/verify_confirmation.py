"""Audit saved confirmation with stdlib; optional replay stays in a new directory.

Default does not sample. The two redundant raw natural-case JSON files are
reconstructed byte-for-byte from RESULTS.json for their registered hashes.
"""
from pathlib import Path
import argparse,hashlib,json,math,statistics,shutil,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent
def digest(data):return hashlib.sha256(data).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--resample',action='store_true',
        help='Replay the frozen plan in a new directory; never replace the saved confirmation.')
    args=parser.parse_args()
    reg=json.loads((HERE/'confirm_registration.json').read_text())
    results=json.loads((HERE/'RESULTS.json').read_text())
    payload={f"natural_v{c['metadata']['v']}.json":
             json.dumps(c,separators=(',',':')).encode() for c in results['cases']}
    checks={}
    for name,h in reg['preserved_prior_files_SHA256'].items():
        data=payload[name] if name in payload else (HERE/name).read_bytes()
        checks[name]=(digest(data)==h)
    checks['confirm.py']=digest((HERE/'confirm.py').read_bytes())==reg['confirm_script_SHA256']
    assert all(checks.values()),checks
    x=json.loads((HERE/'confirm_results.json').read_text())
    assert digest((HERE/'confirm_registration.json').read_bytes())==x['registration_SHA256']
    assert x['script_SHA256']==reg['confirm_script_SHA256']
    d=x['primary_D'];b,rho,alpha=x['b'],x['rho'],x['alpha']
    assert len(d)==reg['N']==1024 and min(d)>=-b and max(d)<=b
    def loge(theta):return math.fsum(math.log1p(rho*(z-theta)/(b+theta)) for z in d)
    lo,hi=-b+1e-10,b-1e-10;target=math.log(1/alpha)
    assert loge(lo)>target>loge(hi)
    for _ in range(100):
        mid=(lo+hi)/2
        if loge(mid)>=target:lo=mid
        else:hi=mid
    assert abs(lo-x['e_inverted_one_sided_lower'])<1e-12
    assert abs(loge(0)-x['log_e_at_theta_0'])<1e-12
    out=dict(N=len(d),mean=statistics.mean(d),SE_descriptive=statistics.stdev(d)/math.sqrt(len(d)),
        log_e_zero=loge(0),inverted_lower=lo,saved_lower_difference=lo-x['e_inverted_one_sided_lower'],
        registered_hash_checks=checks,
        scope='Read saved D only. No new sampling, no interval certification, no log lower bound.')
    def shifted_loge(epsilon,theta):
        return math.fsum(math.log1p(rho*(max(-b,z-epsilon)-theta)/(b+theta)) for z in d)
    sensitivity=[]
    for epsilon in [0.,.01,.05,.1,.2,.3,.4,.5]:
        left,right=-b+1e-8,b-1e-8
        for _ in range(100):
            mid=(left+right)/2
            if shifted_loge(epsilon,mid)>=target:left=mid
            else:right=mid
        sensitivity.append(dict(uniform_downward_error_budget=epsilon,
                               log_e_zero=shifted_loge(epsilon,0.),lower=left))
    left,right=0.,b
    for _ in range(100):
        mid=(left+right)/2
        if shifted_loge(mid,0.)>=target:left=mid
        else:right=mid
    out['rounding_sensitivity']=sensitivity
    out['critical_uniform_downward_budget']=left
    out['sensitivity_caution']='Sensitivity only; no displayed rounding-error budget has been certified.'
    print(json.dumps(out,indent=2))
    if args.resample:
        scratch=Path(tempfile.mkdtemp(prefix='cube033-confirm-replay-'))
        for name in ['natural_probe.py','RESULTS.json','confirm.py','confirm_registration.json']:
            shutil.copy2(HERE/name,scratch/name)
        for name,data in payload.items():(scratch/name).write_bytes(data)
        subprocess.run([sys.executable,str(scratch/'confirm.py'),'--output','replay_results.json'],
                       cwd=scratch,check=True)
        print(json.dumps({'replay_directory':str(scratch),
            'note':'Replay is separate. Retain every replay; it does not replace the unique registered primary result.'}))

if __name__=='__main__':main()
