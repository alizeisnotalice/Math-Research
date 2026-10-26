#!/usr/bin/env python3
"""Finite symbolic/lift/moment certificate. No numerical heat-kernel calls."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import os
import sys
import time

import sympy as sp
from moment import box_statistics, band_mean_factor

HERE = Path(__file__).resolve().parent
CP = Path('/tmp/upper_mean085/checkpoint.json')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def serial(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {k:serial(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [serial(v) for v in value]
    return value


def save(out):
    CP.parent.mkdir(parents=True,exist_ok=True)
    temp=CP.with_suffix('.tmp')
    temp.write_text(json.dumps(serial(out),indent=2)+'\n')
    temp.replace(CP)


def symbolic():
    q,beta,z=sp.symbols('q beta z',real=True)
    scalar=((4+2*beta)*q-2*beta)/(1+q)**2
    major=(beta+2)**2/(4*(beta+1))
    identities={
        'score_square': 1-(2*z/(1+z*z))**2-(z*z-1)**2/(1+z*z)**2,
        'diagonal_upper_square': sp.Rational(9,8)-(6*q-2)/(1+q)**2-(3*q-5)**2/(8*(1+q)**2),
        'beta_majorant_square': major-scalar-((beta+2)*q-(3*beta+2))**2/(4*(beta+1)*(1+q)**2),
        'beta_monotonic_derivative': sp.diff(major,beta)-beta*(beta+2)/(4*(beta+1)**2),
        'face_equality_witness': scalar.subs({beta:sp.Rational(1,2),q:sp.Rational(7,5)})-sp.Rational(25,24),
        'general_equality_witness': scalar.subs({beta:1,q:sp.Rational(5,3)})-sp.Rational(9,8)
    }
    x=sp.Matrix(sp.symbols('x0:4'))
    Q=sp.eye(4)/4-sp.ones(4)/36
    y=sp.Matrix([x[j]-x[0]/6 for j in range(1,4)])
    identities['Q_schur_completion']=(x.T*Q*x)[0]-sp.Rational(5,24)*x[0]**2-(y.T*(sp.eye(3)/4-sp.ones(3)/36)*y)[0]
    for expression in identities.values():
        assert sp.cancel(expression)==0
    assert sp.Rational(25,48)*sp.Rational(2,3)==sp.Rational(25,72)
    assert sp.Rational(25,72)/sp.Rational(7,9)==sp.Rational(25,56)
    assert sp.Rational(25,48)*sp.Rational(24,5)==sp.Rational(20,21)*sp.Rational(21,8)
    return {'status':'seven_exact_symbolic_identities_zero','identities':list(identities),
            'all_zero':True,'symmetric_coefficient':'25/72','symmetric_new_old_ratio':'25/56',
            'universal_curvature_ratio_cap':'20/21'}


def cube_lift(raw):
    u=list(map(F,raw)); r=max(map(abs,u))
    V=[(1,)+tail for tail in product((-1,1),repeat=4)]
    if not r:
        coeff=[F(0)]*16
    else:
        weights={}
        for vertex in product((-1,1),repeat=5):
            weight=F(1)
            for sj,uj in zip(vertex,u):
                weight *= (1+sj*uj/r)/2
            assert weight>=0
            weights[vertex]=weight
        assert sum(weights.values())==1
        coeff=[r*(weights[v]-weights[tuple(-x for x in v)]) for v in V]
    assert [sum((a*v[j] for a,v in zip(coeff,V)),F(0)) for j in range(5)]==u
    assert sum(map(abs,coeff))==r
    assert sum(coeff)==u[0]
    if u[0]==0:
        assert max(map(abs,coeff))<=r/2
    return {'u':u,'radius':r,'coefficients':coeff,'l1_exact':r,
            'first_coordinate_zero':u[0]==0,'face_max_coefficient_cap':r/2 if u[0]==0 else None}


def main():
    begin=time.monotonic()
    reg=json.loads((HERE/'registration.json').read_text())
    assert reg['status']=='registered_before_computation'
    for row in reg['sources']:
        assert sha(row['path'])==row['sha256']
    fingerprint={'registration_sha256':sha(HERE/'registration.json'),
                 'source_sha256':sha(__file__),'moment_source_sha256':sha(HERE/'moment.py')}
    if CP.exists():
        assert '--resume' in sys.argv,'Existing checkpoint: inspect, do not rerun'
        out=json.loads(CP.read_text())
        assert all(out[k]==v for k,v in fingerprint.items())
        assert out['status']!='complete_finite_symbolic_lift_and_moment_certificate'
        out.setdefault('prior_pids',[]).append(out['pid']);out['pid']=os.getpid()
    else:
        out={'id':reg['id'],'status':'armed','pid':os.getpid(),**fingerprint,
             'python':sys.executable,'sympy_version':sp.__version__,'lifts':[],'cases':[]}
    save(out)
    print(json.dumps({k:out[k] for k in ['status','pid','registration_sha256','source_sha256','moment_source_sha256']}),flush=True)
    if '--arm' in sys.argv:input()
    try:
        if 'symbolic' not in out:
            out['symbolic']=symbolic();save(out)
        for raw in reg['lift_vectors'][len(out['lifts']):]:
            out['lifts'].append(serial(cube_lift(raw)));save(out)
        for case in reg['cases'][len(out['cases']):]:
            stats=box_statistics(case['bounds'],case['anchor_tail'])
            if 'expected_m2' in case:
                assert stats['E_max_abs_squared']==F(case['expected_m2'])
            if 'saved_geometry' in case:
                old=case['saved_geometry']
                assert stats['offsets']==list(zip(map(F,old['lo']),map(F,old['hi'])))
                assert stats['mu']==list(map(F,old['mu']))
                assert stats['eta_old_Q']==F(old['eta'])
                assert stats['h_infinity']==F(old['h_inf'])
                assert stats['volume']==F(old['volume'])
            # This is only an arithmetic API check with supplied R=1; no band is evaluated.
            assert band_mean_factor(F(1),F(1),stats)==1
            out['cases'].append({'id':case['id'],**serial(stats)})
            out['status']='saved_case_prefix';save(out)
        out['status']='complete_finite_symbolic_lift_and_moment_certificate'
        out['seconds_including_input_loading_and_arm_wait']=time.monotonic()-begin
        out['new_H_gradient_Green_calls']=0
        out['new_Arb_or_transcendental_calls']=0
        out['no_wall_time_cap']=True
        save(out)
        (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'status':out['status'],'pid':out['pid'],'cases':len(out['cases']),
                          'lifts':len(out['lifts']),'seconds':out['seconds_including_input_loading_and_arm_wait']}),flush=True)
    except Exception as error:
        out['status']='stopped_with_diagnostic';out['error']=repr(error);save(out);raise


if __name__=='__main__':main()
