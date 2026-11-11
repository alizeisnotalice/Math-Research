"""r084 acceptance: exact rational replay, factorial hashes, two independent CAS traces."""
import argparse,json,hashlib,re,multiprocessing as mp
from pathlib import Path
from decimal import Decimal as D,localcontext
from fractions import Fraction as Q
from verify84 import constants,job
from slabs84 import job as slabjob

def number(x):return D(re.sub(r'`[^*]*','',str(x)).replace('*^','E'))
def qdec(x):
    x=Q(x);return D(x.numerator)/D(x.denominator)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('output');ap.add_argument('--project-root',required=True);args=ap.parse_args();p=Path(args.directory);project=Path(args.project_root)
    a=json.loads((p/'python/certificates84.json').read_text());s=json.loads((p/'python/slabs84.json').read_text());w=json.loads((p/'wolfram/verification84.json').read_text())
    assert a['pass_all'] and s['pass_all'] and w['pass'] and w['kernels_closed'] and all(w['symbolic'].values()) and w['exact_exponential_intervals']
    assert a['constants']==constants() and len(a['rows'])==len(w['rows'])==12 and len(s['rows'])==len(w['slab_checks'])==20
    assert [z['N'] for z in a['rows']]==[2**k for k in range(10,22)]
    with mp.get_context('spawn').Pool(12) as pool:
        replay=pool.map(job,[z['N'] for z in a['rows']]);sr=pool.map(slabjob,[(z['N'],z['p']) for z in s['rows']])
    for x,y in zip(a['rows'],replay):
        assert x['factorial']==y['factorial']
        for r,t in zip(x['rows'],y['rows']):
            for f in ['p','eta','parameter_layer_volume_lower','log_elementary_lower_interval','log_factorial_lower_interval']:assert r[f]==t[f]
    for x,y,z in zip(s['rows'],sr,w['slab_checks']):
        for f in ['N','p','probability_interval','proved_probability_lower','exact_probability_numerator_bits','exact_probability_denominator_bits','exact_probability_SHA256']:assert x[f]==y[f]
        assert z['N']==x['N'] and z['p']==x['p'] and z['exact_pass']
        lo,hi=map(Q,x['probability_interval']);assert lo>=Q(x['proved_probability_lower']) and hi<=1
    diff=[];cross=[];width=[];fields=['log_elementary_budget_lower','log_factorial_budget_lower','log_sharpened_budget_lower_numerical','log_actual_modular_upper','log_physical_residual_ratio_upper']
    for x,y in zip(a['rows'],w['rows']):
        n=x['N'];assert n==y['N']
        for row in x['rows']:
            p0=row['p'];wr=[r for r in y['rows'] if r['p']==p0];pr=row['precision_rows'];assert [z['working_precision'] for z in wr]==[80,160]==[z['working_precision'] for z in pr]
            for f in fields:
                diff.extend([abs(number(wr[0][f])-number(wr[1][f])),abs(number(pr[0][f])-number(pr[1][f]))]);cross.append(abs(number(wr[-1][f])-number(pr[-1][f])))
                assert all(number(r[f]).is_finite() for r in wr+pr)
            for f,g in [('log_elementary_budget_lower','log_elementary_lower_interval'),('log_factorial_budget_lower','log_factorial_lower_interval')]:
                lo,hi=map(qdec,row[g]);assert lo<hi;width.append(hi-lo)
                for trace in wr+pr:assert lo<number(trace[f])<hi
            assert number(wr[-1]['log_sharpened_budget_lower_numerical'])>number(wr[-1]['log_factorial_budget_lower'])>number(wr[-1]['log_elementary_budget_lower'])>0
            assert number(wr[-1]['log_actual_modular_upper'])<D(3).ln()+n*(D(21)/20).ln()
            assert number(wr[-1]['log_physical_residual_ratio_upper'])<-1000
    assert max(diff)<D('1e-50') and max(cross)<D('1e-50') and max(width)<D('1e-20')
    code=Path(__file__).resolve().parent
    result=dict(pass_all=True,scope='A full pointwise true-corner soft-cover budget obstruction; no new general C_N bound; fixed-N finiteness not settled.',
      theorem_N_min=1024,eta='N^-2',cap=2,pointwise_log_budget_growth='(3/2) N log N - O(N)',saturated_region_log_budget_growth='N log N - O(N)',
      elementary_p2='(5N/924)^N/(9N^2)',elementary_p3='(5N^(3/2)/1848)^N/(18N^3)',
      constant_gates=a['constants']['gates'],Wolfram_symbolic_gates=w['symbolic'],exact_exponential_interval_replay=True,exact_factorial_interval_replay=True,exact_slab_probability_replay=True,
      factorial_certificate_dimensions=[z['N'] for z in a['rows']],exact_slab_cases=20,largest_factorial_bits=a['rows'][-1]['factorial']['bits'],
      Python_configured_workers=12,Python_certificate_actual_processes=a['actual_processes'],Python_slab_actual_processes=s['actual_processes'],Python_replay_actual_processes=len({z['pid'] for z in replay}),
      Wolfram_configured_workers=w['workers'],Wolfram_actual_kernel_ids=sorted({z['kernel_id'] for z in w['rows']}|{z['kernel_id'] for z in w['slab_checks']}),Wolfram_kernels_closed=w['kernels_closed'],
      working_precision=[80,160],stored_significant_digits=65,max_stored_precision_difference=str(max(diff)),max_stored_cross_software_difference=str(max(cross)),max_directed_log_interval_width=str(max(width)),
      MATLAB_scope='Stable floating sanity check of original coarse formula only; no certificate depends on MATLAB. Owned audit MATLAB PID10076 exited.',goal_status='active',pause_request_revoked=True)
    result['source_sha256']={n:sha(code/n) for n in ['verify84.py','verify84.wl','slabs84.py','check84.py','sliver79.py']}
    result['input_sha256']={n:sha(p/n) for n in ['python/certificates84.json','python/slabs84.json','wolfram/verification84.json']}
    result['prior_dependency_sha256']={n:sha(project/n) for n in ['rounds/r028/report.md','rounds/r080/no_budget80.md','rounds/r081/smoothing81.md','rounds/r082/allocation82.md','rounds/r083/exact83.md'] if (project/n).exists()}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')},indent=2))
if __name__=='__main__':
    with localcontext() as ctx:ctx.prec=180;main()
