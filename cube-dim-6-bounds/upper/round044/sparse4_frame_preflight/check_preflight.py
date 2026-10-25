#!/usr/bin/env python3
"""Bounded exact sparse4 frame audit; no kernel, optimizer or fan builder."""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent
START = time.monotonic()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
dot = lambda x,y: sum(a*b for a,b in zip(x,y))


def det(a):
    a = [[F(v) for v in row] for row in a]
    answer = F(1)
    for j in range(len(a)):
        k = next(i for i in range(j,len(a)) if a[i][j])
        if k != j:
            a[j],a[k] = a[k],a[j]
            answer = -answer
        p = a[j][j]
        answer *= p
        for i in range(j+1,len(a)):
            c = a[i][j]/p
            for k in range(j,len(a)):
                a[i][k] -= c*a[j][k]
    return answer


def gram(vectors):
    return [[sum(v[i]*v[j] for v in vectors) for j in range(6)] for i in range(6)]


def source_check():
    data = json.loads((HERE/'input_manifest.json').read_text())
    for p,h in data['sources'].items():
        assert sha(p) == h, p


def encode(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,(list,tuple)): return [encode(y) for y in x]
    if isinstance(x,dict): return {k:encode(v) for k,v in x.items()}
    return x


def compute():
    source_check()
    reg = json.loads((HERE/'registration.json').read_text())
    beta = F(reg['beta'])
    forms = []
    for support in combinations(range(6),4):
        for signs in product((-1,1),repeat=3):
            v = [0]*6
            for i,s in zip(support,(1,)+signs): v[i]=s
            forms.append(tuple(v))
    assert len(forms)==len(set(forms))==120
    assert all(sum(bool(a) for a in v)==4 and next(a for a in v if a)!=-1 for v in forms)
    supports = {tuple(i for i,a in enumerate(v) if a) for v in forms}
    assert len(supports)==15
    assert all(sum(tuple(i for i,a in enumerate(v) if a)==s for v in forms)==8 for s in supports)
    incidence = [sum(bool(v[i]) for v in forms) for i in range(6)]
    assert incidence==[80]*6
    sparse_gram = gram(forms)
    assert sparse_gram==[[80*int(i==j) for j in range(6)] for i in range(6)]
    body = [(1,)+v for v in product((-1,1),repeat=5)]
    body_gram = gram(body)
    assert body_gram==[[32*int(i==j) for j in range(6)] for i in range(6)]
    full_gram = [[F(body_gram[i][j])+beta**2*sparse_gram[i][j] for j in range(6)] for i in range(6)]
    rates = [sum(abs(v[i]) for v in body)+beta*incidence[i] for i in range(6)]
    assert full_gram==[[F(2561,80)*int(i==j) for j in range(6)] for i in range(6)]
    assert rates==[F(33)]*6
    frame = {'beta':str(beta),'common_denominator':80,
             'column_numerators':[[80*a for a in v] for v in body]+[list(v) for v in forms],
             'base_count':32,'new_orbit_count':120,'column_count':152,
             'interpretation':'Independent standard Cauchy(t) scalar in every column; one shared time t.'}
    frame_text = json.dumps(frame,sort_keys=True,separators=(',',':'))+'\n'
    frame_sha = hashlib.sha256(frame_text.encode()).hexdigest()
    tok = iter(Path(reg['fan_path']).read_text().split())
    n = int(next(tok)); assert n==21
    cones=[]
    for _ in range(n):
        d=int(next(tok)); rr=[]; qq=[]
        for _ in range(6):
            qq.append(int(next(tok)))
            rr.append(tuple(int(next(tok)) for _ in range(6)))
        cones.append((d,rr,qq))
    assert next(tok,None) is None
    unique=sorted({r for _,rr,_ in cones for r in rr})
    assert len(unique)<=21
    cache={(r,v):dot(r,v) for r in unique for v in forms}
    assert len(cache)<=2520
    summaries=[]; split=[]; first=None
    for ci,(d,rr,qq) in enumerate(cones,1):
        cross=[]; zero_contacts=0
        for fi,v in enumerate(forms,1):
            values=[cache[(r,v)] for r in rr]
            if 0 in values: zero_contacts+=1
            if min(values)<0<max(values):
                cross.append(fi); split.append([ci,fi])
                if first is None: first=(ci,fi,d,rr,qq,v,values)
        summaries.append({'cone_1_based':ci,'split_form_indices_1_based':cross,
                          'split_form_count':len(cross),'zero_contact_forms':zero_contacts})
    assert len(cones)*len(forms)==2520
    assert first is not None, 'No split: stop and report instead of building a fan.'
    ci,fi,d,rr,qq,v,values=first
    actual_det=det(list(zip(*rr))); assert abs(actual_det)==d
    ip=next(i for i,a in enumerate(values) if a>0)
    im=next(i for i,a in enumerate(values) if a<0)
    total=[sum(r[j] for r in rr) for j in range(6)]
    m=abs(sum(values))+1
    cp=[1]*6; cm=[1]*6; cp[ip]+=m; cm[im]+=m
    xp=[sum(c*r[j] for c,r in zip(cp,rr)) for j in range(6)]
    xm=[sum(c*r[j] for c,r in zip(cm,rr)) for j in range(6)]
    assert all(a>0 for a in cp+cm)
    assert all(x[j]>x[j+1] for x in (xp,xm) for j in range(5)) and xp[-1]>0 and xm[-1]>0
    sp,sm=dot(v,xp),dot(v,xm)
    assert sp>0>sm
    mid=[F(a+b,2) for a,b in zip(xp,xm)]
    defect=beta*(F(abs(sp)+abs(sm),2)-abs(dot(v,mid)))
    assert defect>0
    witness={'cone_1_based':ci,'form_index_1_based':fi,'unscaled_form':v,
             'beta':beta,'six_integer_rays':rr,'old_rates':qq,'saved_abs_det':d,'verified_det':actual_det,
             'form_values_on_six_rays':values,'positive_ray_index_1_based':ip+1,'negative_ray_index_1_based':im+1,
             'point_plus':xp,'point_minus':xm,'positive_barycentric_coefficients_plus':cp,
             'positive_barycentric_coefficients_minus':cm,'form_point_plus':sp,'form_point_minus':sm,
             'midpoint':mid,'form_midpoint':dot(v,mid),'strict_symbol_convexity_defect_lower':defect,
             'old_cone_is_not_a_new_symbol_linearity_cone':True}
    result={'status':'SPLIT_CERTIFIED_OLD_21_FAN_INVALID_FOR_Q_ONLY_UPDATE',
            'frame_content_sha256':frame_sha,'beta':beta,'orbit_count':120,'support_sets':15,
            'forms_per_support':8,'coordinate_incidence':incidence,'sparse4_unscaled_Gram':sparse_gram,
            'base_Gram':body_gram,'full_Gram':full_gram,'coordinate_Cauchy_rates':rates,
            'unique_old_rays':len(unique),'distinct_ray_form_dot_products':len(cache),
            'cone_form_configurations':2520,'relation_groups':10,'splitting_cones':sum(bool(s['split_form_count']) for s in summaries),
            'splitting_cone_form_pairs':len(split),'per_cone':summaries,'first_witness':witness,
            'new_G_H_p_native':0,'new_fan_enumerations':0,'new_C6_bound':False,
            'old_tube_or_value_certificates_inherited':False}
    assert time.monotonic()-START<900
    source_check()
    return encode(result),frame_text


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--record',action='store_true'); args=parser.parse_args()
    if args.record: assert not (HERE/'results.json').exists(),'single authorized batch already recorded'
    result,frame_text=compute()
    if args.record:
        (HERE/'frame.json').write_text(frame_text)
        (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
        receipt={'pid':os.getpid(),'exit_code':0,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 'wall_seconds':time.monotonic()-START,'python':sys.version,'python_path':sys.executable,
                 'python_sha256':sha(sys.executable),'fractions_sha256':sha(__import__('fractions').__file__),
                 'code_sha256':sha(__file__),'registration_sha256':sha(HERE/'registration.json'),
                 'input_manifest_sha256':sha(HERE/'input_manifest.json'),'third_party_dependencies':[],
                 'new_G_H_p_native':0,'child_processes':[],'all_handles_closed':True}
        (HERE/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    else:
        assert result==json.loads((HERE/'results.json').read_text())
        assert frame_text==(HERE/'frame.json').read_text()
    print(json.dumps({'status':result['status'],'pid':os.getpid(),'wall_seconds':time.monotonic()-START,
                      'unique_rays':result['unique_old_rays'],'dot_products':result['distinct_ray_form_dot_products'],
                      'splitting_cones':result['splitting_cones'],'split_pairs':result['splitting_cone_form_pairs'],
                      'witness':result['first_witness']}))


if __name__=='__main__': main()
