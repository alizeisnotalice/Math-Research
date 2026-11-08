"""Final r082 acceptance: replay directed bounds, cross-check engines, prove finite budget lower."""
import json,sys,re,hashlib,math,multiprocessing as mp,argparse
from decimal import Decimal as D,localcontext
from pathlib import Path
from fractions import Fraction as Q
from certify82 import run

def number(x):return D(re.sub(r'`[^*]*','',str(x)).replace('*^','E'))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('output');ap.add_argument('--project-root',required=True);args=ap.parse_args()
    p=Path(args.directory);project=Path(args.project_root);a=json.loads((p/'python/corner82.json').read_text());c=json.loads((p/'python/certificates82.json').read_text());w=json.loads((p/'wolfram/verification82.json').read_text())
    assert w['pass'] and w['kernels_closed'] and all(w['symbolic'].values()) and len(a['rows'])==len(w['rows'])==30 and len(c['rows'])==len(w['certificate_checks'])==32
    with mp.get_context('spawn').Pool(12) as pool:replay=pool.map(run,c['rows'])
    for x,y in zip(c['rows'],replay):assert x['mean_interval']==y['mean_interval'] and x['partition_cells']==y['partition_cells']
    idx={(r['N'],r['A'],r['side']):r for r in c['rows']};prec=[];cross=[];tail=[];res=[]
    for x,y in zip(a['rows'],w['rows']):
        assert (x['N'],x['A'])==(y['N'],y['A']) and [r['precision'] for r in y['rows']]==[60,120]
        for field in ['q','log_charge_saving','gaussian_reference']:
            prec.append(abs(number(y['rows'][0][field])-number(y['rows'][1][field])))
            cross.append(abs(D(str(x[field]))-number(y['rows'][1][field])))
        for r in y['rows']:
            assert number(r['log_charge_saving'])>0
            tail.append(number(r['integration_tail_upper']));res.append(abs(number(r['mean_residual'])))
            if x['A']==2:assert number(r['relaxed_increment'])>0
            else:assert r['relaxed_increment'] is None and r['relaxed_threshold'] is None
        if x['N']<=256:
            q=number(y['rows'][-1]['q']);assert D(idx[x['N'],x['A'],'lower']['candidate'])<q<D(idx[x['N'],x['A'],'upper']['candidate'])
    for r in w['certificate_checks']:
        assert r['inside_directed_interval'] and r['correct_sign'];v=number(r['mean']);x=idx[r['N'],r['A'],r['side']]
        assert D(x['mean_interval'][0])<v<D(x['mean_interval'][1])
    assert max(prec)<D('1e-35') and max(cross)<D('1e-8') and max(tail)<D('3e-35') and max(res)<D('1e-35')
    # The older explicit cap-two finite source gives log(score)>71.79 (r054).
    assert D(idx[64,3,'upper']['candidate'])==D('56.91')
    expo=Q('71.79')-Q('56.91')-Q(1,64);assert expo==Q(118915,8000)
    terms=[expo**k/Q(math.factorial(k)) for k in range(81)];elo=sum(terms,Q(0));assert elo>2850000
    report=dict(pass_all=True,symbolic_gates=w['symbolic'],dimensions=15,caps=[2,3],Python_base_root_evaluations=30,Mathematica_base_root_evaluations=60,Mathematica_relaxed_root_evaluations=30,
      strict_directed_root_brackets=16,directed_root_endpoints=32,final_integral_cells=sum(x['partition_cells'] for x in c['rows']),max_partition_cells=max(x['partition_cells'] for x in c['rows']),
      directed_arithmetic_precision=70,complete_directed_replay=True,Mathematica_working_precision=[60,120],Mathematica_certificate_endpoint_precision=90,
      maximum_precision_difference=str(max(prec)),maximum_float64_high_precision_difference=str(max(cross)),
      maximum_omitted_mean_integral_tail_upper=str(max(tail)),maximum_calibration_numerical_residual=str(max(res)),
      acceptance_tolerance=dict(precision_difference='1e-35',float64_cross='1e-8',omitted_mean_tail='3e-35'),
      finite_soft_budget_log_lower=str(expo),finite_soft_budget_lower_integer=elo.numerator//elo.denominator,finite_exp_lower_Taylor_terms=81,
      Python_directed_processes=c['actual_processes'],Python_directed_seconds=c['seconds'],Mathematica_workers=w['workers'],Mathematica_seconds=w['seconds'],Mathematica_kernels_closed=True,
      scope='Improved capped fixed-physical-corner upper charge, exact within radial Fenchel/Markov class. Disproves r081 proposed soft-budget scale. Not a new general C_N upper/lower.')
    code=Path(__file__).resolve().parent
    report['source_sha256']={n:hashlib.sha256((code/n).read_bytes()).hexdigest() for n in ['corner82.py','certify82.py','verify82.wl','check82.py']}
    report['input_sha256']={n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in ['python/corner82.json','python/certificates82.json','wolfram/verification82.json']}
    report['prior_proof_sha256']={'rounds/r054/scalar_proof54.md':hashlib.sha256((project/'rounds/r054/scalar_proof54.md').read_bytes()).hexdigest()}
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if not k.endswith('sha256')},indent=2))
if __name__=='__main__':
    with localcontext() as ctx:ctx.prec=160;main()
