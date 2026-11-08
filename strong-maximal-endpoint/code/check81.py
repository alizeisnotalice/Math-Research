"""Acceptance gates for final r081 Python and Wolfram outputs."""
import json,sys,re,hashlib
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal,localcontext

def dec(s):return Decimal(re.sub(r'`[^*]*','',str(s)).replace('*^','E'))
def compare(a,b,field,errors):errors.append(abs(dec(a[field])-dec(b[field])))
def main():
    p=Path(sys.argv[1]);dest=Path(sys.argv[2]);a=json.loads((p/'python/soft81.json').read_text());e=json.loads((p/'python/exterior81.json').read_text());w=json.loads((p/'wolfram/verification81.json').read_text())
    assert w['pass'] and w['kernels_closed'] and all(w['symbolic'].values())
    assert len(a['cover'])==len(w['cover'])==23 and len(a['dilation'])==len(w['dilation'])==len(e['rows'])==len(w['exterior'])==15
    prec=[];cross=[];cells=0;max_width=Q(0)
    for x,y in zip(a['cover'],w['cover']):
        assert x['k']==y['k'] and y['exact_pass'] and x['prefix_families']==y['prefix_families']==x['k']+4
        assert Q(x['eta'])==Q(1,2**x['k']);cells+=x['prefix_families']
        assert Q(x['tail_weighted_node_sum_upper'])==Q(16*x['tail_nodes'],2**256)
        lo,hi=map(Q,x['budget_interval']);assert 0<lo<hi;max_width=max(max_width,hi-lo)
        for r,s in zip(x['precision_rows'],y['rows']):
            assert r['precision']==s['precision']
            for f in ['prefix_budget','normalized_by_eta_sqrt']:compare(r,s,f,cross)
        for f in ['prefix_budget','normalized_by_eta_sqrt']:
            compare(x['precision_rows'][0],x['precision_rows'][1],f,prec);compare(y['rows'][0],y['rows'][1],f,prec)
    for x,y in zip(a['dilation'],w['dilation']):
        assert x['N']==y['N'] and Q(x['eta'])==Q(1,x['N']**2)
        for r,s in zip(x['rows'],y['rows']):
            assert r['precision']==s['precision'] and r['cap_branch'] and s['cap_branch']
            assert dec(r['log_dilation_plus'])<=Decimal(1)/x['N'] and dec(r['log_dilation_minus'])<Decimal(4)/(3*x['N'])
            for f in ['log_dilation_plus','log_dilation_minus','cap_branch_increment','tau_critical']:compare(r,s,f,cross)
        for f in ['log_dilation_plus','log_dilation_minus','cap_branch_increment','tau_critical']:
            compare(x['rows'][0],x['rows'][1],f,prec);compare(y['rows'][0],y['rows'][1],f,prec)
    for x,y in zip(e['rows'],w['exterior']):
        assert x['N']==y['N'] and y['exact_pass'];lo,hi=map(Q,x['normalized_volume_interval']);assert 1<lo<hi
        for r,s in zip(x['rows'],y['rows']):
            assert r['precision']==s['precision'] and 0<dec(r['smooth_to_hard_volume'])<1
            for f in ['normalized_volume_series','log_exterior_volume','log_exterior_ratio_upper','smooth_to_hard_volume','scaled_smoothing_deficit','log_modular_relative_error_upper']:compare(r,s,f,cross)
        for f in ['normalized_volume_series','log_exterior_volume','log_exterior_ratio_upper','smooth_to_hard_volume','scaled_smoothing_deficit','log_modular_relative_error_upper']:
            compare(x['rows'][0],x['rows'][1],f,prec);compare(y['rows'][0],y['rows'][1],f,prec)
    assert max(prec)<Decimal('1e-30'),max(prec)
    assert max(cross)<Decimal('1e-30'),max(cross)
    assert max_width<Q(1,10**68),float(max_width)
    r=dict(pass_all=True,symbolic_gates=w['symbolic'],soft_cover_scales=23,certified_prefix_families=cells,
      exact_weight_and_budget_intervals=True,maximum_total_budget_interval_width=str(float(max_width)),
      dimensions=15,maximum_dimension=a['dilation'][-1]['N'],precision=[60,120],acceptance_tolerance='1e-30 absolute; Wolfram NIntegrate goals floor(precision/2)',
      cover_precision_evaluations_per_engine=46,dilation_precision_evaluations_per_engine=30,exterior_quadratures_per_engine=30,
      precision_difference_max=str(max(prec)),independent_engine_difference_max=str(max(cross)),
      Python_exterior_workers=len(e['worker_pids']),Python_exterior_seconds=e['seconds'],
      Wolfram_workers=w['workers'],Wolfram_seconds=w['seconds'],Wolfram_kernels_closed=w['kernels_closed'],
      last_scaled_exterior_smoothing_deficit=e['rows'][-1]['rows'][-1]['scaled_smoothing_deficit'],
      scope='Exact smoothing-norm comparison and fixed-dimensional true-column mechanism results. Exterior-sector volume only. No new general C_N upper/lower; sharp coefficient gap unchanged.')
    r['source_sha256']={n:hashlib.sha256((Path(__file__).resolve().parent/n).read_bytes()).hexdigest() for n in ['soft81.py','exterior81.py','verify81.wl','check81.py']}
    r['input_sha256']={n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in ['python/soft81.json','python/exterior81.json','wolfram/verification81.json']}
    dest.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if not k.endswith('sha256')},indent=2))
if __name__=='__main__':
    with localcontext() as c:c.prec=160;main()
