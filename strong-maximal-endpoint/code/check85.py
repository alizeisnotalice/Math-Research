"""Final r085 exact replay and independent Wolfram/Python acceptance."""
import argparse,json,hashlib,re,multiprocessing as mp
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal as D,localcontext
from verify85 import constants85,job

def num(x):return D(re.sub(r'`[^*]*','',str(x)).replace('*^','E'))
def qnum(x):
    q=Q(x);return D(q.numerator)/D(q.denominator)
def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('output');ap.add_argument('--project-root',required=True);a=ap.parse_args();d=Path(a.directory);pr=Path(a.project_root)
    p=json.loads((d/'python/certificates85.json').read_text());w=json.loads((d/'wolfram/verification85.json').read_text());assert p['pass_all'] and w['pass'] and w['kernels_closed'] and w['exact_root_bracket'] and all(w['symbolic'].values());assert p['constants']==constants85()
    assert len(p['rows'])==len(w['rows'])==12 and [r['N'] for r in p['rows']]==[2**k for k in range(10,22)]
    with mp.get_context('spawn').Pool(12) as pool:rr=pool.map(job,[(r['N'],p['constants']['mu_interval']) for r in p['rows']])
    dif=[];cross=[];width=[];fields=['c','r','p','entropy','log_binomial','log_sphere_upper','entropy_limit','canonical_lost_mean_upper']
    for x,y,z in zip(p['rows'],w['rows'],rr):
        for f in ['N','k','c_interval','r_interval','p_interval','k_threshold_interval','log_binomial_interval','log_sphere_upper_interval','canonical_lost_mean_upper','binomial_integer']:assert x[f]==z[f]
        assert x['N']==y['N'] and x['k']==y['k'] and y['exact_integer_and_interval_checks']
        assert Q(x['k_threshold_interval'][1])<=x['k'] and 0<x['k']<x['N']//2
        assert Q(x['canonical_lost_mean_upper'])<Q(1,x['N']**2)
        py=x['precision_rows'];wl=y['rows'];assert [r['working_precision'] for r in py]==[80,160]==[r['working_precision'] for r in wl]
        for f in fields:
            dif.extend([abs(num(py[0][f])-num(py[1][f])),abs(num(wl[0][f])-num(wl[1][f]))]);cross.append(abs(num(py[1][f])-num(wl[1][f])))
            assert all(num(r[f]).is_finite() for r in py+wl)
        for f,g in [('c','c_interval'),('r','r_interval'),('p','p_interval'),('log_binomial','log_binomial_interval'),('log_sphere_upper','log_sphere_upper_interval')]:
            lo,hi=map(qnum,x[g]);assert lo<hi
            if f=='log_binomial':width.append(hi-lo)
            for row in py+wl:assert lo<num(row[f])<hi
        assert num(py[1]['entropy'])>num(py[1]['entropy_limit'])
    assert max(dif)<D('1e-50') and max(cross)<D('1e-50') and max(width)<D('1e-20')
    code=Path(__file__).resolve().parent
    out=dict(pass_all=True,scope='Exact exponential modular cost of retaining canonical witnesses; alternative-rectangle success is outside this round. No new C_N bound.',
      exact_root_and_constant_replay=True,constant_gates=p['constants']['gates'],Wolfram_symbolic_gates=w['symbolic'],exact_binomial_integer_replay=True,
      dimensions=[r['N'] for r in p['rows']],largest_binomial_bits=max(r['binomial_integer']['bits'] for r in p['rows']),
      Python_configured_workers=12,Python_actual_processes=p['actual_processes'],Python_replay_actual_processes=len({r['pid'] for r in rr}),Wolfram_configured_workers=12,Wolfram_actual_kernel_ids=sorted({r['kernel_id'] for r in w['rows']}),Wolfram_kernels_closed=w['kernels_closed'],
      working_precision=[80,160],stored_significant_digits=65,max_stored_precision_difference=str(max(dif)),max_stored_cross_software_difference=str(max(cross)),max_directed_log_binomial_interval_width=str(max(width)),
      exact_critical_entropy_interval=p['constants']['critical_entropy_interval'],exact_exponent_gap_interval=p['constants']['log_mu_interval'],goal_status='active',pause_request_revoked=True)
    out['source_sha256']={f:sha((pr/'code'/f) if f in ['verify84.py','sliver79.py'] else code/f) for f in ['verify85.py','verify85.wl','check85.py','verify84.py','sliver79.py']}
    out['input_sha256']={f:sha(d/f) for f in ['python/certificates85.json','wolfram/verification85.json']}
    out['prior_dependency_sha256']={f:sha(pr/f) for f in ['rounds/r010/report.md','rounds/r011/report.md','rounds/r028/report.md','rounds/r069/report.md','rounds/r083/exact83.md','rounds/r084/obstruction84.md']}
    Path(a.output).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if not k.endswith('sha256')},indent=2))
if __name__=='__main__':
    with localcontext() as c:c.prec=180;main()
