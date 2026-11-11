"""Final r087 acceptance; conditional numeric bounds do not instantiate unknown r022 constants."""
import argparse,json,hashlib,re,multiprocessing as mp
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal as D,localcontext
from verify87 import constants87,job

def num(x):return D(re.sub(r'`[^*]*','',str(x)).replace('*^','E'))
def qnum(x):
    z=Q(x);return D(z.numerator)/D(z.denominator)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('output');ap.add_argument('--project-root',required=True);a=ap.parse_args();d=Path(a.directory);pr=Path(a.project_root);p=json.loads((d/'python/conditional87.json').read_text());w=json.loads((d/'wolfram/verification87.json').read_text())
    assert p['pass_all'] and w['pass'] and w['kernels_closed'] and w['geometry_intervals'] and all(w['symbolic'].values());assert p['constants']==constants87();assert len(p['rows'])==48
    valid=[r for r in p['rows'] if r['certified_formula']];assert len(valid)==len(w['rows'])==w['conditional_validated_cases']==35 and w['excluded_inputs']==13
    with mp.get_context('spawn').Pool(12) as pool:replay=pool.map(job,[(r['N'],r['gamma'],p['constants']) for r in p['rows']])
    for x,y in zip(p['rows'],replay):assert {k:v for k,v in x.items() if k!='pid'}=={k:v for k,v in y.items() if k!='pid'}
    diff=[];cross=[];width=[]
    for x,y in zip(valid,w['rows']):
        assert (x['N'],x['gamma'])==(y['N'],y['gamma']) and y['pass'];py=x['precision_rows'];wl=y['rows'];assert [r['working_precision'] for r in py]==[80,160]==[r['working_precision'] for r in wl]
        f='log_conditional_efficiency_upper';diff.extend([abs(num(py[0][f])-num(py[1][f])),abs(num(wl[0][f])-num(wl[1][f]))]);cross.append(abs(num(py[1][f])-num(wl[1][f])));lo,hi=map(qnum,x['log_efficiency_upper_interval']);width.append(hi-lo)
        for z in py+wl:assert num(z[f]).is_finite() and lo<num(z[f])<hi
        assert all(z['inside_directed_interval'] for z in wl)
    assert max(diff)<D('1e-45') and max(cross)<D('1e-45') and max(width)<D('1e-20')
    code=Path(__file__).resolve().parent;deps=['verify85.py','verify86.py','verify84.py','sliver79.py']
    o=dict(pass_all=True,scope='Uniform receiver-restricted fixed-physical-corner efficiency theorem on a near-Q shell; no new general C_N bound. Finite-N numeric bounds are conditional on the explicitly assumed scalar deficit.',
      theorem_log_efficiency_upper='-(N/2) log N + (N/2) log log N + O(N)',constant_gates=p['constants']['gates'],Wolfram_symbolic_gates=w['symbolic'],geometry_interval_replay=True,conditional_exact_log_interval_replay=True,
      conditional_input_count=48,conditional_certified_cases=35,conditional_excluded_inputs=13,unknown_r022_constant_instantiated=False,certified_finite_N_onset=False,
      Python_configured_workers=12,Python_actual_processes=p['actual_processes'],Python_replay_actual_processes=len({r['pid'] for r in replay}),Wolfram_configured_workers=12,Wolfram_actual_kernel_ids=sorted({r['kernel_id'] for r in w['rows']}),Wolfram_kernels_closed=w['kernels_closed'],
      working_precision=[80,160],stored_significant_digits=65,max_stored_precision_difference=str(max(diff)),max_stored_cross_software_difference=str(max(cross)),max_directed_log_interval_width=str(max(width)),new_general_C_N_bound=False,goal_status='active',pause_request_revoked=True)
    o['source_sha256']={f:sha(pr/'code'/f if f in deps else code/f) for f in ['verify87.py','verify87.wl','check87.py']+deps};o['input_sha256']={f:sha(d/f) for f in ['python/conditional87.json','wolfram/verification87.json']}
    prior=['rounds/r022/report.md','archive/work/strong-maximal-r022-audit.md','rounds/r085/core85.md','rounds/r086/free_rectangles86.md','rounds/r083/exact83.md'];o['prior_dependency_sha256']={f:sha(pr/f) for f in prior};Path(a.output).write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({k:v for k,v in o.items() if not k.endswith('sha256')},indent=2))
if __name__=='__main__':
    with localcontext() as c:c.prec=180;main()
