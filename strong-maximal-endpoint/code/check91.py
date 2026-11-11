"""Replay r091 exact geometry and compare independent local numerical implementations."""
import argparse, hashlib, json, multiprocessing as mp, re, sys
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
from pathlib import Path
TOL=D('1e-45')
FIELDS=('mu','sigma','L','z','t0','t1','log_p0_over_h','log_p1_over_error',
        't0_over_z','t1_over_z','logPhi1_minus_d','logPhi2_minus_2d','fee_minus_d','fee_minus_d_over_L')
def dec(x):
    z=D(re.sub(r'`[^*]*','',str(x)).replace('*^','E'));assert z.is_finite();return z
def qdec(x):
    q=Q(x);return D(q.numerator)/D(q.denominator)
def rel(x,y):
    x,y=dec(x),dec(y)
    return D(0) if x==y==0 else abs(x-y)/max(abs(x),abs(y))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def resource(z):
    assert z['exit_code']==0 and z['stop_reason'] is None
    assert D(str(z['reserve_GiB']))>=1
    assert D(str(z['min_system_available_MiB']))>=D(str(z['reserve_GiB']))*1024
    assert D(str(z['peak_owned_RSS_MiB']))<=D(str(z['owned_RSS_budget_GiB']))*1024
    assert 'estimate' in z['scope'].lower() and z['samples']>0
    return z
def sans_pid(rows):return [{k:v for k,v in x.items() if k!='pid'} for x in rows]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('output');ap.add_argument('--project-root',required=True);a=ap.parse_args()
    d,root=Path(a.directory),Path(a.project_root);code=Path(__file__).resolve().parent
    sys.path.insert(0,str(root/'code'));sys.path.insert(0,str(code))
    from experiment91 import STRIPE_P,STRIPE_M,UNIT_M,SCALAR_K,stripe_job,unit_audit_job,scalar_job
    from verify85 import constants85
    names=('python/experiments91.json','python/resources91.json','wolfram/verification91.json','wolfram/resources91.json')
    py,pyres,wf,wfres=[json.loads((d/f).read_text()) for f in names]
    assert py['pass_all'] is True and wf['pass'] is True
    assert py['working_precision']==[80,160] and py['configured_workers']==12
    assert 1<=py['actual_workers']<=12 and 1<=len(set(py['worker_pids']))<=12
    assert len(wf['symbolic'])==12 and all(wf['symbolic'].values())
    assert wf['kernels_closed'] is True and 1<=wf['configured_workers']<=12
    assert 0<wf['per_case_memory_limit_MiB']<=1024
    pyresource,wfresource=resource(pyres),resource(wfres)
    mu=constants85()['mu_interval'];assert py['constants']['mu_interval']==mu
    assert py['counts']=={'stripe_cases':16,'stripe_rectangle_tests':64,'unit_audit_cases':4,'scalar_dimensions':12,'scalar_precision_rows':24}
    with mp.get_context('spawn').Pool(12) as pool:
        stripes=pool.map(stripe_job,[(str(p),n) for p in STRIPE_P for n in STRIPE_M])
        units=pool.map(unit_audit_job,UNIT_M)
        scalars=pool.map(scalar_job,[(k,mu) for k in SCALAR_K])
    for name,rows in [('stripe_rows',stripes),('unit_audit_rows',units),('scalar_rows',scalars)]:
        assert sans_pid(py[name])==sans_pid(rows),name
    ids={r['kernel_id'] for name in ('stripe_rows','unit_audit_rows','scalar_rows') for r in wf[name]}
    assert all(isinstance(i,int) and i>0 for i in ids) and len(ids)<=wf['configured_workers']
    exact_map={'base_density':'base_density','low_density':'low_density','density_error':'density_error',
               'proof_bound_2K0eps_over_bx':'proof_bound','actual_mean_1p25base_minus_0p25low':'actual_mean','strict_margin':'strict_margin'}
    count=nonzero=success=0;maxprice=D(0)
    def check_rect(p,w):
        nonlocal count,nonzero,success
        assert p['rect']==w['rectangle'] and p['actual_witness'] is w['actual_witness'] and w['pass'] is True
        for x,y in exact_map.items():assert Q(p[x])==Q(w[y]),(x,p[x],w[y])
        low,base=Q(p['low_density']),Q(p['base_density']);assert 0<=low<=base<=1
        err=Q(p['density_error']);assert err<=Q(p['proof_bound_2K0eps_over_bx'])
        assert (Q(p['strict_margin'])>0)==p['actual_witness']
        lo,hi=map(qdec,p['rectangle_price_interval']);assert lo<=hi
        v=qdec(low)+(qdec(base)-qdec(low))*D('1.25')*(1+D('1.25').ln())
        assert lo<=v<=hi
        assert [x['working_precision'] for x in p['precision_rows']]==[80,160]
        for row in p['precision_rows']:assert abs(dec(row['local_rectangle_price'])-v)<D('1e-60')
        count+=1;nonzero+=err>0;success+=p['actual_witness']
    wm={(Q(x['p']),x['m']):x for x in wf['stripe_rows']}
    assert len(wm)==16
    for case in py['stripe_rows']:
        w=wm[Q(case['p']),case['m']];assert w['pass'] is True
        assert [Q(x['low_share']) for x in case['global_low_shares']]==[Q(x) for x in w['global_low_shares']]==[Q(case['p'])]*2
        assert len(case['rectangles'])==len(w['rectangles'])==4
        for p,r in zip(case['rectangles'],w['rectangles']):assert p['id']==r['id'];check_rect(p,r)
        lo,hi=map(qdec,case['global_price_interval']);assert lo<=hi
        price=dec(w['price']);assert lo-D('1e-60')<=price<=hi+D('1e-60')
        for pr in case['global_price_precision_rows']:
            dif=rel(pr['global_price'],price);assert dif<TOL;maxprice=max(maxprice,dif)
    wm={x['m']:x for x in wf['unit_audit_rows']};assert len(wm)==4
    for case in py['unit_audit_rows']:
        w=wm[case['m']];assert w['pass'] is True and Q(w['global_low_share'])==Q(1,3)
        r=case['rectangles'][0];check_rect(r,w['rectangle_result'])
        assert Q(r['density_error'])==Q(r['audit_expected_density_error'])
        assert Q(r['actual_mean_1p25base_minus_0p25low'])==Q(r['audit_expected_mean'])
    assert count==68 and nonzero>0 and 0<success<count
    assert [x['k'] for x in py['scalar_rows']]==list(SCALAR_K)==[x['k'] for x in wf['scalar_rows']]
    maxrel=D(0);maxby={f:D(0) for f in FIELDS};within=D(0);maxres=D(0);tail_checks=[]
    import mpmath as m
    for p,w in zip(py['scalar_rows'],wf['scalar_rows']):
        assert p['N']==str(1<<p['k'])
        assert [x['working_precision'] for x in p['precision_rows']]==[80,160]==[x['working_precision'] for x in w['precision_rows']]
        for rows in (p['precision_rows'],w['precision_rows']):
            for field in FIELDS:
                dif=rel(rows[0][field],rows[1][field]);assert dif<TOL;within=max(within,dif)
        for pr,wr in zip(p['precision_rows'],w['precision_rows']):
            assert wr['pass'] is True and dec(wr['log_tail_residual'])<D(10)**(-wr['working_precision'])
            maxres=max(maxres,dec(wr['log_tail_residual']))
            for field in FIELDS:
                dif=rel(pr[field],wr[field]);assert dif<TOL,(p['k'],field,dif)
                maxby[field]=max(maxby[field],dif);maxrel=max(maxrel,dif)
            assert D(0)<dec(pr['t0'])<dec(pr['t1']) and D(0)<dec(pr['tau'])<D('.125')
            assert D(0)<dec(pr['p0'])<1 and D(0)<dec(pr['p1'])<1
            with m.workdps(110):
                z=m.mpf(pr['z']);s=m.mpf(pr['sigma']);L=m.mpf(pr['L'])
                residual=abs(m.log(m.erfc(z/m.sqrt(2))/2)-m.log(s*z)+L/2)
                assert residual<m.mpf('1e-55')
        pr=p['precision_rows'][-1];tail_checks.append({k:pr[k] for k in ('L','t0_over_z','t1_over_z','log_p0_over_h','log_p1_over_error','fee_minus_d_over_L')})
    for field in ('p0','p1'):
        vals=[dec(p['precision_rows'][-1][field]) for p in py['scalar_rows']]
        assert all(a>b for a,b in zip(vals,vals[1:])),field
    own=('experiment91.py','verify91.wl','check91.py');dep=('math_guard.py','verify85.py','verify84.py','sliver79.py')
    prior=['rounds/r011/report.md','rounds/r021/report.md','rounds/r022/report.md','rounds/r044/report.md','rounds/r089/core89.md','rounds/r090/core90.md','rounds/r090/audit90.md']
    output={'pass_all':True,'scope':'Exact two-dimensional finite geometry/price intervals and scalar tail/fee diagnostics; no finite-N physical theorem onset or new general norm coefficient.',
      'counts':py['counts'],'rational_rectangle_tests':count,'nonzero_boundary_errors':nonzero,'strict_success_tests':success,'failed_strict_success_tests':count-success,
      'Python_configured_workers':12,'Python_actual_stored_workers':py['actual_workers'],'Python_replay_scalar_workers':len({r['pid'] for r in scalars}),
      'Wolfram_configured_workers':wf['configured_workers'],'Wolfram_actual_kernel_ids':sorted(ids),'Wolfram_kernels_closed':True,'symbolic_gates':wf['symbolic'],
      'working_precision':[80,160],'relative_tolerance':str(TOL),'max_scalar_cross_relative_difference':str(maxrel),
      'max_within_software_precision_difference':str(within),'max_scalar_relative_by_field':{f:str(v) for f,v in maxby.items()},
      'max_global_price_cross_relative_difference':str(maxprice),'max_wolfram_log_tail_residual':str(maxres),'scalar_scale_samples':tail_checks,
      'Python_resource_monitor':pyresource,'Wolfram_resource_monitor':wfresource,'resources_are_sampled_estimates_not_OS_reservations':True,
      'source_sha256':{f:sha((code if f in own else root/'code')/f) for f in own+dep},'input_sha256':{f:sha(d/f) for f in names},
      'document_sha256':{f:sha(d/f) for f in ('core91.md','audit91.md')},'prior_dependency_sha256':{f:sha(root/f) for f in prior},
      'new_general_C_N_coefficient':False,'goal_status':'active','pause_request_revoked':True}
    Path(a.output).write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if 'sha256' not in k and k!='scalar_scale_samples'},indent=2))
if __name__=='__main__':
    with localcontext() as ctx:ctx.prec=180;main()
