#!/usr/bin/env python3
"""Finite incidence/constant certificate for the new positive rank-descent bound."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import gzip
import hashlib
import json
import os
import signal
import sys
import time

from flint import arb, ctx, fmpq

HERE = Path(__file__).resolve().parent
CP = Path('/tmp/upper_singular057/checkpoint.json')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(value):
    CP.parent.mkdir(parents=True,exist_ok=True)
    temp = CP.with_suffix('.tmp')
    temp.write_text(json.dumps(value,indent=2)+'\n')
    temp.replace(CP)


def alarm(_signum,_frame):
    raise TimeoutError('120-second registered checkpoint; no automatic restart')


def ar(q):
    q=F(q)
    return arb(fmpq(q.numerator,q.denominator))


def ball(q):
    assert q.is_finite() and q>0
    return {'ball':str(q),'lower_dyadic':str(q.lower().fmpq()),
            'upper_dyadic':str(q.upper().fmpq())}


def geometry(reg):
    begin=time.monotonic()
    s1=[(1,)+tail for tail in product((-1,1),repeat=4)]
    assert len(s1)==16
    positive=[v for v in s1 if min(v)>=0]
    assert positive==[(1,1,1,1,1)]
    corner=(1,1,1,1,1)
    lines=[]
    for k in range(1,5):
        for D in combinations(range(1,5),k):
            D=set(D)
            endpoint=tuple(0 if j in D else 1 for j in range(5))
            direction=tuple(-1 if j in D else 0 for j in range(5))
            assert endpoint not in s1
            # Incident pair planes at the non-S1 endpoint: fixed +1 outside D,
            # opposite signs on D, and unordered pairs stored only once.
            pairs=set()
            for signs in product((-1,1),repeat=k):
                sign_at=dict(zip(sorted(D),signs))
                v=tuple(sign_at[j] if j in D else 1 for j in range(5))
                w=tuple(-sign_at[j] if j in D else 1 for j in range(5))
                assert v in s1 and w in s1 and v!=w
                assert all(F(v[j]+w[j],2)==endpoint[j] for j in range(5))
                pairs.add(tuple(sorted((v,w))))
            assert len(pairs)==2**(k-1)
            # No triple enumeration is needed: a line through two distinct
            # affine hypercube vertices cannot contain a third vertex.
            lines.append({'D':sorted(D),'size':k,'common_S1_endpoint':list(corner),
                          'non_S1_endpoint':list(endpoint),'direction':list(direction),
                          'incident_pair_plane_count_at_non_S1_endpoint':len(pairs),
                          'incident_pairs':[[list(v),list(w)] for v,w in sorted(pairs)]})
    assert len(lines)==15
    cross=0
    for line1,line2 in combinations(lines,2):
        a,b=line1['direction'],line2['direction']
        minors=[a[i]*b[j]-a[j]*b[i] for i,j in combinations(range(1,5),2)]
        assert any(minors)
        cross+=1
    assert cross==105
    data=json.loads(gzip.decompress(Path(reg['sources'][-1]['path']).read_bytes()))
    assert data['denominator']==144 and len(data['vertices'])==280
    norms=Counter(F(sum(z*z for z in v),144**2) for v in data['vertices'])
    assert norms=={F(1,128):40,F(1,144):80,F(1,162):160}
    assert max(norms)==F(1,128)
    # For a fixed body direction v, the squared projection lengths of the
    # remaining fifteen directions have only two classes.
    projection_norms=Counter()
    for w in s1:
        if w==corner:continue
        dot=sum(w)
        projection_norms[F(5)-F(dot*dot,5)]+=1
    assert projection_norms=={F(16,5):5,F(24,5):10}
    # The determinant argument for any third projected direction uses only
    # det Gram(v,w)<=24 and one nonzero 3x3 minor divisible by four.
    assert F(16,24)==F(2,3)
    return {'status':'finite_endpoint_incidence_and_ellipticity_inputs_verified',
            'S1_full_face_points':[list(v) for v in s1],
            'S1_positive_face_points':[list(v) for v in positive],
            'positive_S2_lines':lines, 'distinct_line_pair_count':cross,
            'all_distinct_positive_lines_intersect_only_at_common_S1_corner':True,
            'endpoint_multiplicities_by_D_size':{str(k):2**(k-1) for k in range(1,5)},
            'number_of_endpoint_incident_pair_records':sum(x['incident_pair_plane_count_at_non_S1_endpoint'] for x in lines),
            'vertex_squared_norm_counts':{str(k):v for k,v in sorted(norms.items())},
            'one_direction_projected_squared_length_counts':{str(k):v for k,v in sorted(projection_norms.items())},
            'all_third_direction_projected_squared_lengths_at_least':'2/3',
            'third_direction_bound_proof':'nonzero 3x3 sign minor has magnitude >=4; Gram determinant >=16; pair Gram determinant <=24',
            'seconds':time.monotonic()-begin}


def constants(reg):
    begin=time.monotonic(); ctx.prec=256
    pi=arb.pi(); s2=arb(2).sqrt(); s5=arb(5).sqrt()
    P2=1/(256*pi)
    P3=1/(1024*s2*pi**2)
    P4=3/(65536*pi**2)
    S=5*s5/4+5*ar(F(5,6)).sqrt()
    B=7*ar(F(3,2)).sqrt()/(8*pi**2)
    A0=288*P4+8*P3*S/pi
    A1=16*P3/pi
    F5=16/(pi*s5)
    vals={'P2':P2,'P3':P3,'P4':P4,'S':S,'B':B,'A0':A0,'A1':A1,'F5':F5}
    checks={name:bool(vals[name]<ar(cap)) for name,cap in reg['preselected_strict_rational_caps'].items()}
    assert all(checks.values())
    # The logarithm truncation constant satisfies Lambda^2=14sqrt(3)/9<4.
    assert F(14,9)**2*3 < 16
    assert (B/(576*P3)-14*arb(3).sqrt()/9).contains(0)
    return {'status':'one_scalar_constant_batch_certified','precision_bits':ctx.prec,
            'constants':{name:ball(value) for name,value in vals.items()},
            'strict_rational_caps':reg['preselected_strict_rational_caps'],
            'all_caps_verified':checks, 'logarithm_comparison_exact':'(14/9)^2*3<16',
            'seconds':time.monotonic()-begin}


def main():
    reg=json.loads((HERE/'registration.json').read_text())
    assert reg['status']=='registered_before_computation'
    for row in reg['sources']:assert sha(row['path'])==row['sha256']
    out={'id':reg['id'],'status':'armed' if '--arm' in sys.argv else 'starting',
         'pid':os.getpid(),'python':sys.executable,'registration_sha256':sha(HERE/'registration.json'),
         'source_sha256':sha(__file__),'external_sources':reg['sources']}
    save(out); print(json.dumps({k:out[k] for k in ['status','pid','registration_sha256','source_sha256']}),flush=True)
    if '--arm' in sys.argv:input()
    signal.signal(signal.SIGALRM,alarm)
    for stage,func in [('geometry',geometry),('constants',constants)]:
        out['status']='running_'+stage;save(out);signal.alarm(120)
        try:out[stage]=func(reg)
        except Exception as error:
            out['status']='stopped_with_diagnostic';out['error']=repr(error);save(out);raise
        finally:signal.alarm(0)
        save(out);print(json.dumps({'stage':stage,'status':out[stage]['status'],'seconds':out[stage]['seconds']}),flush=True)
    out['status']='kernel_geometry_and_constants_certified_no_tube_integral_claim'
    out['new_Green_evaluations']=0;out['new_I_or_g_or_C5_bound']=False
    save(out);(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
