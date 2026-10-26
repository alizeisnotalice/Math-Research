#!/usr/bin/env python3
"""Audit saved period3 mass/endpoint/derivative certificates, without union reruns."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
import sys,json,os,time
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from period3_variation import ENV,sha,load,save,rt,enc

def main():
    start=time.monotonic();reg=load(P/'registration.json')
    assert sha(P/'period3_variation.py')==reg['source_sha256']
    for row in reg['inputs']:assert sha(Path(row['path']))==row['sha256']
    for oldround in ('lower_round021','lower_round022'):
        old=ROOT/oldround
        for row in load(old/'archive_whitelist.json')['files']:assert sha(old/row['path'])==row['sha256']
    F,u,v,bounds,_,_=field();a=(1+u**3/2,)+((1+u)/2,)*4
    def decode(d):
        out=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));out[5*i+j]=Q(c)
        return F(a=out)
    @lru_cache(None)
    def radius(n):
        t=n;e2=e3=0
        while t%2==0:t//=2;e2+=1
        while t%3==0:t//=3;e3+=1
        assert t==1
        z=u**e2*v**e3/2;assert 32*z**5==F(n);return z
    C=decode(load(ROOT/'lower_round003/a2_equals_a3_results.json')['coefficients'])
    audit=[];derivatives=[];total_atoms=0
    for direction,sign in (('plus',1),('minus',-1)):
        g=load(P/f'{direction}_geometry.json');r=load(P/f'{direction}_derivative.json');ck=load(P/f'{direction}_checkpoint.json')
        assert r['status']=='STRICT_NEGATIVE_DIRECTIONAL_DERIVATIVE'
        assert r['geometry_sha256']==sha(P/f'{direction}_geometry.json')
        assert ck['status']=='COMPLETE_two_union_engines' and ck['completed_first_axis_slabs']==ck['total_first_axis_slabs']==85
        catalog=g['mass_phase_catalog'];assert len(catalog)==3072
        assert len({(tuple(x['lengths']),x['phase']) for x in catalog})==3072
        for row in catalog:
            ls=row['lengths'];phase=row['phase'];n=prod(ls)
            direct=sign*sum((2,-1,-1)[j%3] for j in range(phase,phase+ls[0]))*prod(ls[1:])
            assert row['N']==n and row['mass_derivative']==direct
            assert row['first_axis_block_sum']==direct//prod(ls[1:]) and abs(direct)<=2*n
        endpoint_checks=0;pointchecks=0
        for row in g['kept_jets']:
            ls=row['lengths'];starts=row['starts'];n=prod(ls);dm=row['mass_derivative'];rad=radius(n);dr=rad*dm/(5*n)
            assert row['distinct_atom_count']==n and row['phase']==starts[0]%3
            atoms=set(product(*(range(j,j+ell) for j,ell in zip(starts,ls))))
            assert len(atoms)==n
            assert sum(sign*(2,-1,-1)[point[0]%3] for point in atoms)==dm
            pointchecks+=n
            for k,(j,ell) in enumerate(zip(starts,ls)):
                aa=a[k];domain=3*aa if k==0 else aa/2
                lo=max((F(0),F(0)),((j+ell-1)*aa-rad,-dr))
                hi=min((domain,F(0)),(j*aa+rad,dr))
                assert lo<hi
                for pair,endpoint in zip((lo,hi),row['endpoints'][2*k:2*k+2]):
                    assert endpoint['value']==enc(pair[0]) and endpoint['speed']==enc(pair[1])
                    assert endpoint['faces']==[enc(pair[1]) if t==k else {} for t in range(5)]
                    endpoint_checks+=1
        partial=ck['partial_volume_jet'];assert decode(partial['value'])*Q(16,3)==C
        derivative=decode(r['directional_derivative_coefficients']);derivatives.append(derivative)
        assert decode(partial['speed'])*Q(16,3)==derivative
        faces=[decode(x['coefficients']) for x in r['synchronous_coordinate_face_contributions']]
        assert sum(faces,F(0))==derivative
        for k,x in enumerate(faces):
            assert decode(partial['faces'][k])*Q(16,3)==x
            assert list(map(str,bounds(x.a)))==r['synchronous_coordinate_face_contributions'][k]['interval']
        assert list(map(str,bounds(derivative.a)))==r['derivative_interval'] and Q(r['derivative_interval'][1])<0
        assert r['counts']['active_raw']==179 and r['counts']['raw_candidates']==6144
        audit.append({'direction':direction,'mass_phase_rows_checked':3072,'kept_boxes_checked':len(g['kept_jets']),'endpoint_jets_checked':endpoint_checks,'explicit_distinct_atom_and_mass_sum_checks':pointchecks,'strict_negative_derivative_interval':r['derivative_interval'],'result_sha256':sha(P/f'{direction}_derivative.json')});total_atoms+=pointchecks
    assert derivatives[0]==derivatives[1]
    out={'status':'PASS_saved_period3_directional_certificates','checks':audit,'two_direction_derivative_coefficients_exactly_equal':True,'all021022_archive_files_unchanged':True,'old_source_and_failure_certificates_unchanged':True,'new_union_jobs':0,'actual_total_distinct_atom_checks':total_atoms,'no_finite_epsilon_candidate_started':True,'runtime':rt(start),'source_sha256':sha(Path(__file__))}
    save(P/'verification.json',out);print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
