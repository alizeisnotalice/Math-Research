#!/usr/bin/env python3
"""Finite plan first, then ONLY missing algebraic ratio points; zero H/log."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json,os,time
import pathmean as P

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pack(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return list(map(pack,x))
    return x
def write(p,x):p.write_text(json.dumps(pack(x),sort_keys=True,indent=2)+'\n')

def main():
    started=time.monotonic();print('pathmean111_pid',os.getpid(),flush=True)
    reg=json.loads((HERE/'registration.json').read_text())
    for p,h in reg['sources'].items():assert sha(p)==h,p
    oldcatpath=next(p for p in reg['sources'] if p.endswith('upper_adaptive103/ratio_catalog.json'))
    old91path=next(p for p in reg['sources'] if p.endswith('upper_shift091/results.json'))
    catalog=json.loads(Path(oldcatpath).read_text());origins={k:'saved103' for k in catalog}
    for row in json.loads(Path(old91path).read_text())['ratios']:
        if row['x'] not in catalog:
            catalog[row['x']]={k:row['R_new'][k] for k in ('lower','upper')}
            origins[row['x']]='saved091'
    catalog['0']={'lower':'1','upper':'1'};origins['0']='exact_identity'
    plans=[P.darboux_plan(case['offsets'],case['a'],case['panels_per_piece']) for case in reg['finite_cases']]
    allx=sorted({x for p in plans for x in p['ratio_nodes']})
    missing=[str(x) for x in allx if str(x) not in catalog]
    scalar={'task':'upper_pathmean111.finite_scalar_batch','registered_before_new_ratio_calls':True,
            'registration_sha256':sha(HERE/'registration.json'),'pathmean_sha256':sha(HERE/'pathmean.py'),
            'finite_x':list(map(str,allx)),'new_x':missing,
            'reused_x':{str(x):origins[str(x)] for x in allx if str(x) in catalog},
            'new_formula_count':len(missing),'bits':reg['new_ratio_precision_bits'],
            'new_H_jet_log_exp_calls':0,'wallcap':None}
    target=HERE/'scalar_registration.json'
    if target.exists():assert json.loads(target.read_text())==pack(scalar)
    else:write(target,scalar)
    write(HERE/'geometry_plans.json',plans)
    print('registered_new_ratio_count',len(missing),'reused',len(allx)-len(missing),flush=True)
    api_path=next(p for p in reg['sources'] if p.endswith('upper_shift091/ratio.py'))
    spec=importlib.util.spec_from_file_location('readonly_ratio091',api_path)
    api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)
    new=[]
    for x in missing:
        lo,hi=api.ratio_interval(F(x),reg['new_ratio_precision_bits'])
        catalog[x]={'lower':str(lo),'upper':str(hi)};origins[x]='new111'
        new.append({'x':x,'lower':str(lo),'upper':str(hi)})
    cases=[]
    for case,plan in zip(reg['finite_cases'],plans):
        result=P.evaluate_plan(plan,catalog)
        assert result['factor_upper']<=result['old085_factor_upper']
        if plan['profile']['h']>0:assert result['strict_curvature_improvement']
        cases.append({'name':case['name'],'result':result,
                      'strict_mean_factor_improvement':result['factor_upper']<result['old085_factor_upper']})
    used={str(x):{**catalog[str(x)],'origin':origins[str(x)]} for x in allx}
    write(HERE/'ratio_catalog.json',used)
    result={'status':'passed','pid':os.getpid(),'wall_seconds':time.monotonic()-started,
        'registration_sha256':sha(HERE/'registration.json'),'scalar_registration_sha256':sha(target),
        'pathmean_sha256':sha(HERE/'pathmean.py'),'certify_sha256':sha(__file__),
        'geometry_plans_sha256':sha(HERE/'geometry_plans.json'),'catalog_sha256':sha(HERE/'ratio_catalog.json'),
        'new_ratio_points':new,'new_formula_count':len(new),'reused_count':len(allx)-len(new),
        'cases':cases,'new_H_jet_log_exp_asinh_calls':0,'old_ratio_recomputations':0,
        'full_tree_repack':False,'global_I_or_C5_claim':False}
    write(HERE/'results.json',result)
    print(json.dumps({'status':'passed','pid':result['pid'],'wall_seconds':result['wall_seconds'],
        'new_ratio_points':len(new),'reused_points':len(allx)-len(new),
        'strict_mean_factor_improvements':sum(c['strict_mean_factor_improvement'] for c in cases)}))

if __name__=='__main__':main()
