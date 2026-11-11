"""Final r083 acceptance: strict directed replay and independent completed numeric traces."""
import json,sys,re,hashlib,math,multiprocessing as mp,argparse
from decimal import Decimal as D,localcontext
from fractions import Fraction as Q
from pathlib import Path
from certify83 import run

def number(x):return D(re.sub(r'`[^*]*','',str(x)).replace('*^','E'))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('output');ap.add_argument('--project-root',required=True);args=ap.parse_args();p=Path(args.directory);project=Path(args.project_root)
    a=json.loads((p/'python/exact83.json').read_text());c=json.loads((p/'python/certificates83.json').read_text());w=json.loads((p/'wolfram/verification83.json').read_text())
    assert w['pass'] and w['kernels_closed'] and all(w['symbolic'].values())
    assert len(a['rows'])==len(w['rows'])==30 and len(c['rows'])==len(w['certificate_checks'])==32 and len(c['brackets'])==16
    with mp.get_context('spawn').Pool(12) as pool:rep=pool.map(run,c['rows'])
    for x,y in zip(c['rows'],rep):
        for f in ['mean_interval','scaled_modular_interval','log_norm_at_candidate_interval','partition_cells']:assert x[f]==y[f]
    idx={(z['N'],z['A'],z['side']):z for z in c['rows']};bi={(z['N'],z['A']):z for z in c['brackets']};pr=[];cr=[];tailm=[];tailc=[];res=[]
    for x,y in zip(a['rows'],w['rows']):
        assert (x['N'],x['A'])==(y['N'],y['A']) and [r['working_precision'] for r in y['rows']]==[60,120]
        for f in ['activation','log_true_corner_norm','scaled_source_modular']:
            pr.append(abs(number(y['rows'][0][f])-number(y['rows'][1][f])));cr.append(abs(D(str(x[f]))-number(y['rows'][1][f])))
        for r in y['rows']:
            assert number(r['scaled_source_modular'])>1 and number(r['radial_Fenchel_gap'])>0
            tailm.append(number(r['mean_tail_upper']));tailc.append(number(r['scaled_cost_tail_upper']));res.append(abs(number(r['mean_residual'])))
            if x['A']==2:assert number(r['relaxed_increment'])>0
            else:assert r['relaxed_increment'] is None
        if x['N']<=256:
            b=bi[x['N'],x['A']]
            assert D(b['activation_interval'][0])<number(y['rows'][-1]['activation'])<D(b['activation_interval'][1])
            assert D(b['log_true_norm_interval'][0])<number(y['rows'][-1]['log_true_corner_norm'])<D(b['log_true_norm_interval'][1])
    for r in w['certificate_checks']:
        assert r['all_inside'] and r['correct_mean_sign'];x=idx[r['N'],r['A'],r['side']]
        for f,g in [('mean','mean_interval'),('scaled_modular','scaled_modular_interval'),('log_norm_at_candidate','log_norm_at_candidate_interval')]:assert D(x[g][0])<number(r[f])<D(x[g][1])
    assert max(pr)<D('1e-35') and max(cr)<D('1e-8') and max(res)<D('1e-35') and max(tailm)<D('1e-300') and max(tailc)<D('1e-300')
    cap3=bi[64,3];assert D(cap3['log_true_norm_interval'][1])<D('54.57')
    expo=Q('71.79')-Q('54.57')-Q(1,64);assert expo==Q(27527,1600)
    esum=sum((expo**k/Q(math.factorial(k)) for k in range(81)),Q(0));assert esum>29_600_000
    result=dict(pass_all=True,symbolic_gates=w['symbolic'],completed_parameter_cases=30,base_mean_roots_at_two_precisions=60,relaxed_threshold_roots_at_two_precisions=30,
      verified_dimension_max=1048576,caps=[2,3],directed_root_brackets=16,directed_endpoints=32,complete_directed_replay=True,
      directed_precision=70,directed_final_cells=sum(x['partition_cells'] for x in c['rows']),directed_max_cells=max(x['partition_cells'] for x in c['rows']),
      Wolfram_working_precision=[60,120],Wolfram_endpoint_precision=90,max_precision_difference=str(max(pr)),max_float64_cross_difference=str(max(cr)),
      max_mean_tail_upper=str(max(tailm)),max_scaled_cost_tail_upper=str(max(tailc)),max_numeric_mean_residual=str(max(res)),
      acceptance_tolerances=dict(precision='1e-35',float64_cross='1e-8',mean_and_scaled_cost_tail='1e-300'),
      certified_log_corner64_cap3_upper='54.57',finite_soft_budget_log_lower=str(expo),finite_soft_budget_lower_integer=esum.numerator//esum.denominator,
      positive_exp_Taylor_terms=81,Python_directed_processes=c['actual_processes'],Python_directed_seconds=c['seconds'],Wolfram_workers=w['workers'],
      Wolfram_cached_completed_parameter_cases=w['cached_completed_cases'],Wolfram_final_certificate_phase_seconds=w['seconds'],Wolfram_kernels_closed=True,
      excluded_extension_dimensions=[2**30,2**40,2**60],excluded_extension_status='Each exceeded a 300-second per-case limit; integer-shape Gamma errors observed; none is counted as verified.',
      scope='Actual fixed physical-corner norm via existing r028 r=0 reduction. Constant-order asymptotic and radial-method loss3/2; not a new general strong-maximal upper bound. User requested pause after this round.')
    code=Path(__file__).resolve().parent;result['source_sha256']={n:hashlib.sha256((code/n).read_bytes()).hexdigest() for n in ['exact83.py','certify83.py','verify83.wl','check83.py','certify82.py']}
    result['input_sha256']={n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in ['python/exact83.json','python/certificates83.json','wolfram/verification83.json']}
    result['prior_dependency_sha256']={n:hashlib.sha256((project/n).read_bytes()).hexdigest() for n in ['rounds/r028/report.md','rounds/r082/asymptotic82.md','rounds/r082/wolfram/verification82.json','rounds/r054/scalar_proof54.md']}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')},indent=2))
if __name__=='__main__':
    with localcontext() as ctx:ctx.prec=180;main()
