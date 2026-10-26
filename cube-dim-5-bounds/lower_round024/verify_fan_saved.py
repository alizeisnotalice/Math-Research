#!/usr/bin/env python3
"""Saved raw-phase, endpoint-response and fan-certificate audit; no union jobs."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
import os,sys,time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from profile_variation import ENV,load,save,sha,rt,enc
def main():
    start=time.monotonic();F,u,v,bounds,_,_=field();zero=F(0)
    def decode(d):
        out=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));out[5*i+j]=Q(c)
        return F(a=out)
    @lru_cache(None)
    def radius(n):
        x=n;i=j=0
        while x%2==0:x//=2;i+=1
        while x%3==0:x//=3;j+=1
        assert x==1
        r=u**i*v**j/2;assert 32*r**5==F(n);return r
    a=(1+u**3/2,)+((1+u)/2,)*4
    stats=[];C=decode(load(ROOT/'lower_round003/a2_equals_a3_results.json')['coefficients'])
    specs=[('middle',load(P/'registration.json'))]+[(f'cell{i}',load(P/'fan_registration.json')) for i in range(6)]
    affine_pattern=((Q(1,2),Q(-1,2)),(Q(0),Q(1)),(Q(-1,2),Q(-1,2)))
    for name,reg in specs:
        for x in reg['inputs']:assert sha(Path(x['path']))==x['sha256']
        r=load(P/f'{name}_derivative.json');g=load(P/f'{name}_geometry.json');ck=load(P/f'{name}_checkpoint.json')
        pattern=tuple(map(Q,reg['profiles'][name]));assert sum(pattern)==0
        assert r['geometry_sha256']==sha(P/f'{name}_geometry.json')
        assert ck['status']=='COMPLETE_two_union_engines' and ck['completed_first_axis_slabs']==85
        assert len(g['mass_phase_catalog'])==3072
        for x in g['mass_phase_catalog']:
            ls=x['lengths'];j=x['phase'];dm=sum(pattern[(j+t)%3] for t in range(ls[0]))*prod(ls[1:])
            assert Q(x['mass_derivative'])==dm and x['N']==prod(ls)
        endpoint_checks=0;atomchecks=0
        for x in g['kept_jets']:
            ls=x['lengths'];starts=x['starts'];n=prod(ls);dm=Q(x['mass_derivative']);r0=radius(n);d=r0*dm/(5*n)
            atoms=set(product(*(range(j,j+ell) for j,ell in zip(starts,ls))))
            assert len(atoms)==n and sum((pattern[z[0]%3] for z in atoms),Q(0))==dm;atomchecks+=n
            j=starts[0];tail=prod(ls[1:]);dmA=sum(affine_pattern[(j+t)%3][0] for t in range(ls[0]))*tail;dmB=sum(affine_pattern[(j+t)%3][1] for t in range(ls[0]))*tail
            dA=r0*dmA/(5*n);dB=r0*dmB/(5*n)
            for k,(j,ell) in enumerate(zip(starts,ls)):
                domain=3*a[k] if k==0 else a[k]/2
                lo_const=(j+ell-1)*a[k]-r0;hi_const=j*a[k]+r0
                low_pair=max((zero,zero),(lo_const,-d));high_pair=min((domain,zero),(hi_const,d))
                assert low_pair<high_pair
                for pair,ep,sgn in ((low_pair,x['endpoints'][2*k],-1),(high_pair,x['endpoints'][2*k+1],1)):
                    assert ep['value']==enc(pair[0]) and ep['speed']==enc(pair[1])
                    expected=[enc(pair[1]) if q==k else {} for q in range(5)]
                    if name!='middle':
                        moving=(pair==(lo_const,-d)) if sgn==-1 else (pair==(hi_const,d))
                        expected += [enc(dA*sgn),enc(dB*sgn)] if moving else [{},{}]
                    assert ep['faces']==expected;endpoint_checks+=1
        derivative=decode(r['directional_derivative_coefficients'])
        assert decode(ck['partial_volume_jet']['value'])*Q(16,3)==C
        assert decode(ck['partial_volume_jet']['speed'])*Q(16,3)==derivative
        assert list(map(str,bounds(derivative.a)))==r['derivative_interval'] and Q(r['derivative_interval'][1])<0
        if name!='middle':
            A=decode(r['fan_cell_linear_coefficients']['constant']);B=decode(r['fan_cell_linear_coefficients']['slope']);t=Q(r['sample_t'])
            assert A+B*t==derivative
            assert decode(ck['partial_volume_jet']['faces'][5])*Q(16,3)==A and decode(ck['partial_volume_jet']['faces'][6])*Q(16,3)==B
        stats.append({'profile':name,'mass_phase_rows':3072,'endpoint_response_checks':endpoint_checks,'explicit_atom_checks':atomchecks})
    events=load(P/'events.json')['events']
    for event in events:
        n=decode(event['numerator']);d=decode(event['positive_denominator']);assert zero<d
        witness=event['comparison_witness'];A,B=map(decode,witness['left_velocity']);C,D=map(decode,witness['right_velocity'])
        assert (A-C)*d+(B-D)*n==zero
    fan=load(P/'fan_results.json');assert fan['normalized_fan_maximum']['strict_negative']
    assert Q(fan['normalized_fan_maximum']['interval'][1])<0
    out={'status':'PASS_saved_full_fan_profiles_and_responses','profiles':stats,'event_witness_equations_checked':len(events),'closed_continuity_and_maximum_source':'fan_results.json and finalize_fan.py','new_union_jobs':0,'no_finite_epsilon_candidate_started':True,'runtime':rt(start),'source_sha256':sha(Path(__file__)),'fan_result_sha256':sha(P/'fan_results.json')}
    save(P/'verification.json',out);print({k:v for k,v in out.items() if k!='profiles'})
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
