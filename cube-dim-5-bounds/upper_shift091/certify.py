#!/usr/bin/env python3
"""Twelve lifts, six new ratio formulas and six new log qualifications."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import sys
import time
from flint import arb,ctx,fmpq
from ratio import chain_lift,ratio_interval,new_ratio_upper,combine_ratio_upper

HERE=Path(__file__).resolve().parent
CP=Path('/tmp/upper_shift091/checkpoint.json')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def serial(v):
    if isinstance(v,F):return str(v)
    if isinstance(v,dict):return {k:serial(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [serial(x) for x in v]
    return v


def save(out):
    CP.parent.mkdir(parents=True,exist_ok=True)
    tmp=CP.with_suffix('.tmp');tmp.write_text(json.dumps(serial(out),indent=2)+'\n');tmp.replace(CP)


def ar(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))


def ball(x):
    assert x.is_finite() and x>0
    lo,hi=F(str(x.lower().fmpq())),F(str(x.upper().fmpq()))
    assert 0<lo<hi
    return {'ball':str(x),'lower':str(lo),'upper':str(hi),
            'lower_display':float(lo),'upper_display':float(hi)}


def matched_source(row):
    data=json.loads(Path(row['source_path']).read_text())
    cell=next(c for c in data['cells'] if c['h']==row['source_cell_h'])
    band=next(b for b in cell['bands'] if b['start']==row['source_time_start'])
    assert cell['u_lower']==row['source_u_lower'] and cell['u_upper']==row['source_u_upper']
    h=max(abs(F(v)) for key in ['u_lower','u_upper'] for v in cell[key])
    assert h==F(cell['h'])
    assert h/F(band['start'])==F(row['x'])
    pre='R_old_' if 'R_old_lower' in band else 'R_'
    assert band[pre+'lower']==row['old_asinh_lower']
    assert band[pre+'upper']==row['old_asinh_upper']
    assert band.get('R_lip_lower')==row['old_exp_lower']
    assert band.get('R_lip_upper')==row['old_exp_upper']


def main():
    wall=time.monotonic();reg=json.loads((HERE/'registration.json').read_text())
    assert reg['status']=='registered_before_computation'
    for src in reg['sources']:assert sha(src['path'])==src['sha256']
    fp={'registration_sha256':sha(HERE/'registration.json'),
        'source_sha256':sha(__file__),'ratio_source_sha256':sha(HERE/'ratio.py')}
    if CP.exists():
        assert '--resume' in sys.argv,'Inspect existing checkpoint; do not replay'
        out=json.loads(CP.read_text());assert all(out[k]==v for k,v in fp.items())
        assert out['status']!='complete_chain_ratio_certificate'
        out.setdefault('prior_pids',[]).append(out['pid']);out['pid']=os.getpid()
    else:out={'id':reg['id'],'status':'armed','pid':os.getpid(),'python':sys.executable,**fp,'lifts':[],'ratios':[]}
    save(out);print(json.dumps({k:out[k] for k in ['status','pid','registration_sha256','source_sha256','ratio_source_sha256']}),flush=True)
    if '--arm' in sys.argv:input()
    begin=time.monotonic()
    try:
        for tail in reg['lift_cases'][len(out['lifts']):]:
            lift=chain_lift(tail)
            assert sum(lift['gaps'])==lift['radius']
            out['lifts'].append(serial(lift));save(out)
        assert ratio_interval(0)==(F(1),F(1))
        ctx.prec=384
        for i in range(len(out['ratios']),len(reg['fixed_ratios'])):
            x=F(reg['fixed_ratios'][i]);comparison=reg['saved_matched_comparisons'][i]
            assert x==F(comparison['x']);matched_source(comparison)
            X=ar(x)
            R=((1+X*X/16).sqrt()+X/4)**2*((1+X*X/256).sqrt()+X/16)**8
            logR=R.log()
            rb,lb=ball(R),ball(logR)
            lower,upper=ratio_interval(x,bits=reg['pure_integer_api_bits'])
            assert lower<=F(rb['lower'])<F(rb['upper'])<=upper
            assert F(lb['upper'])<x
            oldlo,oldup=F(comparison['old_asinh_lower']),F(comparison['old_asinh_upper'])
            assert upper<oldlo
            if comparison['old_exp_lower'] is not None:
                exp_lower=F(comparison['old_exp_lower']);exp_upper=F(comparison['old_exp_upper'])
                assert upper<exp_lower
                relative_exp=[lower/exp_upper,upper/exp_lower]
            else:relative_exp=None
            h,a=F(comparison['source_cell_h']),F(comparison['source_time_start'])
            assert new_ratio_upper(h,a)==upper
            assert combine_ratio_upper(h,a,[oldup])==upper
            row={'x':str(x),'R_new':rb,'log_R_new':lb,
                 'integer_api_interval':[str(lower),str(upper)],
                 'old_asinh_interval':[str(oldlo),str(oldup)],
                 'new_upper_below_saved_old_asinh_lower':True,
                 'log_upper_strictly_below_x':True,
                 'saved_exp_comparison':comparison['old_exp_lower'] is not None,
                 'new_to_saved_old_asinh_interval':[str(lower/oldup),str(upper/oldlo)],
                 'new_to_saved_exp_interval':serial(relative_exp),
                 'source_h':str(h),'source_time':str(a)}
            out['ratios'].append(row);out['status']='saved_scalar_prefix';save(out)
            print(json.dumps({'x':str(x),'R_new_display':float(upper),'log_R_new_display':lb['upper_display'],
                              'new_to_saved_exp_upper_display':float(relative_exp[1]) if relative_exp else None}),flush=True)
        out['status']='complete_chain_ratio_certificate';out['precision_bits']=384
        out['seconds_computation']=time.monotonic()-begin
        out['seconds_including_input_loading_and_arm_wait']=time.monotonic()-wall
        out['new_algebraic_R_formulas']=6;out['new_log_qualifications']=6
        out['old_R_exp_asinh_recomputations']=0;out['new_H_stencil_gradient_Green_calls']=0
        out['no_wall_time_cap']=True
        save(out);(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'status':out['status'],'pid':out['pid'],'seconds_computation':out['seconds_computation']}),flush=True)
    except Exception as error:
        out['status']='stopped_with_diagnostic';out['error']=repr(error);save(out);raise


if __name__=='__main__':main()
