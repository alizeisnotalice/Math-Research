#!/usr/bin/env python3
"""One exact affine-order check for the selected positive spacing step."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import os
import time
import shape_jets as S

HERE=Path(__file__).resolve().parent


def field_from_dict(values):
    c=[Q(0)]*36
    for ij,x in values.items():
        i,j=map(int,ij.split(','))
        c[6*i+j]=Q(x)
    return S.F(a=c)


def main():
    started=time.monotonic()
    paths=[HERE/'jet_results.json',HERE/'shape_jets.py',HERE/'step_registration.md',
           S.ARITHMETIC,Path(__file__),
           S.MAIN/'lower/round003/second_spacing_breakpoint_certificate_results.json']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    results=json.loads(paths[0].read_text())
    direction=results['positive_direction']
    assert direction['h']==[0,1,0,0,0,0] and direction['D']==direction['E']==0
    step=Q(1,1000000)
    raw,order,b,family=S.baseline_family(S.spacings())
    moving=order.index(1)
    aj=b[moving]
    endpoints={S.J(0),S.J(aj/2,aj/2)}
    for record in raw:
        q=record['starts'][moving]
        ell=record['ell'][moving]
        r=record['radius']
        endpoints.add(S.J((q+ell-1)*aj-r,(q+ell-1)*aj))
        endpoints.add(S.J(q*aj+r,q*aj))
    ordered=sorted(endpoints)
    separations=[]
    for left,right in zip(ordered,ordered[1:]):
        at_step=(right.c-left.c)+(right.d-left.d)*step
        lo,hi=S.A.bounds(at_step.a)
        assert lo>0,'selected finite step crosses or leaves an unresolved endpoint event'
        separations.append(lo)
    certificate=json.loads(paths[-1].read_text())
    base=field_from_dict(certificate['coefficients'])
    derivative=field_from_dict(direction['derivative_coefficients'])
    candidate=base+derivative*step
    lower,upper=S.A.bounds(candidate.a)
    old_upper=Q(results['current_complete_upper'])
    assert lower>old_upper,'one-step primary candidate fails the complete old upper'
    out={'status':'primary_affine_step_pass_pending_independent_inward_and_L1',
         'pid':os.getpid(),'step':str(step),'source':'unit rectangular lattice',
         'parameter_formula':['1+u^4/2','(1000001/1000000)*(u^3+u^4*v^3)/4',
                              '(1+u)/2','(1+u)/2','(1+u)/2','(1+u)/2'],
         'relations':['u^6=2, u>0','v^6=3, v>0'],
         'chosen_direction':{'h':direction['h'],'D':0,'E':0,'active_axis_for_representation':3},
         'full_raw_baseline_boxes':len(raw),
         'distinct_unclipped_affine_endpoints_moving_axis':len(ordered),
         'strict_adjacent_order_checks':len(separations),
         'minimum_separation_at_step_rational_lower':str(min(separations)),
         'affine_interval':'0<t<=1/1000000',
         'ratio_rule':'R(t)=R(0)+t*R1 for this retained full raw witness union',
         'ratio_coefficients':S.coefficients(candidate),
         'ratio_exact_interval':[str(lower),str(upper)],
         'ratio_float':float((lower+upper)/2),
         'current_complete_old_upper':str(old_upper),
         'strict_primary_margin_lower':str(lower-old_upper),
         'margin_float':float(lower-old_upper),
         'input_sha256':hashes,'union_calls':0,
         'elapsed_seconds':time.monotonic()-started,
         'actual_process_status':'finished; no MATLAB, child process, or step-size scan'}
    (HERE/'positive_candidate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','pid','step','full_raw_baseline_boxes',
          'strict_adjacent_order_checks','ratio_float','margin_float','union_calls','elapsed_seconds')}))


if __name__=='__main__':
    main()
