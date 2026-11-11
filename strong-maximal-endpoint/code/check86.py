"""r086 acceptance: exact action-margin replay and independent clipped-log integrals."""
import argparse,json,hashlib,re,multiprocessing as mp
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal as D,localcontext
from verify86 import constants86,job

def num(x):return D(re.sub(r'`[^*]*','',str(x)).replace('*^','E'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('output');ap.add_argument('--project-root',required=True);a=ap.parse_args();d=Path(a.directory);pr=Path(a.project_root);p=json.loads((d/'python/moments86.json').read_text());w=json.loads((d/'wolfram/verification86.json').read_text())
    assert p['pass_all'] and w['pass'] and w['kernels_closed'] and w['input_guard_digits']==40 and all(w['symbolic'].values());assert p['constants']==constants86();assert len(p['rows'])==len(w['rows'])==72
    with mp.get_context('spawn').Pool(12) as pool:replay=pool.map(job,[(r['endpoint'],r['b']) for r in p['rows']])
    diff=[];cross=[];qerr=[];serr=[];slack=[]
    for x,y,z in zip(p['rows'],w['rows'],replay):
        assert (x['endpoint'],x['b'])==(y['endpoint'],y['b'])==(z['endpoint'],z['b']);assert x['precision_rows']==z['precision_rows']
        py=x['precision_rows'];wl=y['rows'];assert [r['working_precision'] for r in py]==[80,160]==[r['working_precision'] for r in wl]
        serr.extend(map(D,map(str,x['scipy_absolute_errors'])))
        for f in ['mean','log_mgf','moment_slack']:
            diff.extend([abs(num(py[0][f])-num(py[1][f])),abs(num(wl[0][f])-num(wl[1][f]))]);cross.append(abs(num(py[1][f])-num(wl[1][f])))
            assert all(num(r[f]).is_finite() for r in py+wl)
        for r in wl:assert r['pass'];qerr.append(abs(num(r['quadrature_formula_error'])));slack.append(num(r['moment_slack']))
    assert max(diff)<D('1e-40') and max(cross)<D('1e-40') and max(qerr)<D('1e-45') and max(serr)<D('1e-10') and min(slack)>0
    for r in p['constants']['rows']:assert 0<Q(r['mean_loss_gap_interval'][0])<Q(r['mean_loss_gap_interval'][1]);assert 0<Q(r['log_volume_inflation_per_dimension_interval'][0])<Q(r['log_volume_inflation_per_dimension_interval'][1])
    code=Path(__file__).resolve().parent
    o=dict(pass_all=True,scope='All alternative rectangles for one true source piece per relative direction, on the specified near-Q shell; neither general C_N nor fixed-physical-corner cover is settled.',
      exact_action_margin_replay=True,verified_rational_tolerances=[r['a'] for r in p['constants']['rows']],clipped_log_cases=72,endpoint_classes=18,moment_parameters=['1/100','1/25','1/10','1/4'],Wolfram_symbolic_gates=w['symbolic'],
      Python_configured_workers=12,Python_actual_processes=p['actual_processes'],Python_replay_actual_processes=len({r['pid'] for r in replay}),Wolfram_configured_workers=12,Wolfram_actual_kernel_ids=sorted({r['kernel_id'] for r in w['rows']}),Wolfram_kernels_closed=w['kernels_closed'],
      working_precision=[80,160],Wolfram_input_guard_digits=40,Wolfram_stored_significant_digits=45,Python_stored_significant_digits=65,max_stored_precision_difference=str(max(diff)),max_stored_cross_software_difference=str(max(cross)),
      max_Wolfram_quadrature_formula_error=str(max(qerr)),max_SciPy_absolute_error=str(max(serr)),minimum_numeric_moment_slack=str(min(slack)),
      certified_finite_N_onset=False,new_general_C_N_bound=False,MATLAB_scope='Formula sanity only; recorded MATLAB PID11609 exited; session version was not queried.',goal_status='active',pause_request_revoked=True)
    dep=['verify85.py','verify84.py','sliver79.py'];own=['verify86.py','verify86.wl','check86.py','audit86.m'];o['source_sha256']={f:sha(pr/'code'/f if f in dep else code/f) for f in own+dep}
    o['input_sha256']={f:sha(d/f) for f in ['python/moments86.json','wolfram/verification86.json']};o['prior_dependency_sha256']={f:sha(pr/f) for f in ['rounds/r085/core85.md','rounds/r085/validation85.json','rounds/r028/report.md','rounds/r069/report.md','rounds/r083/exact83.md']}
    Path(a.output).write_text(json.dumps(o,indent=2)+'\n');print(json.dumps({k:v for k,v in o.items() if not k.endswith('sha256')},indent=2))
if __name__=='__main__':
    with localcontext() as c:c.prec=180;main()
