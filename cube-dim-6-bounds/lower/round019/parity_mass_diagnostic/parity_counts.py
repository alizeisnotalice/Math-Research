#!/usr/bin/env python3
"""Exactly 12 rational parity-count cases; no geometry/union/search."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import os
import time


def parity(k):
    return 1 if k%2==0 else -1


def main():
    started=time.process_time()
    eps=Q(1,32)
    rest_mass=6
    results=[]
    assert (1+eps)+(1-eps)==2
    for q in (-1,0,1):
        for ell in (1,2,3,4):
            sigma=0 if ell%2==0 else parity(q)
            direct=rest_mass*sum((1+eps*parity(k) for k in range(q,q+ell)),Q(0))
            formula=rest_mass*(ell+eps*sigma)
            negative=rest_mass*sum((1-eps*parity(k) for k in range(q,q+ell)),Q(0))
            neighbor=rest_mass*sum((1+eps*parity(k) for k in range(q+1,q+ell+1)),Q(0))
            reflected_q=-q-ell+1
            reflected=rest_mass*sum((1+eps*parity(k) for k in range(reflected_q,reflected_q+ell)),Q(0))
            assert direct==formula and negative==neighbor and reflected==direct
            assert direct+neighbor==2*rest_mass*ell
            results.append({'start':q,'length':ell,'rest_mass':rest_mass,
                            'parity_sum':sigma,'mass_at_1_32':str(direct),
                            'neighbor_mass_at_1_32':str(neighbor),
                            'relative_radius_derivative_at_zero':str(Q(sigma,6*ell)),
                            'reflection_start':reflected_q,'status':'PASS'})
    assert len(results)==12
    here=Path(__file__).resolve()
    out={'status':'PASS','pid':os.getpid(),'parity_count_cases':12,
         'epsilon':str(eps),'period_mass':2,'results':results,
         'CPU_seconds':time.process_time()-started,
         'full_union_calls':0,'geometry_scans':0,'optimization_scans':0,
         'baseline_replays':0,'new_C6_bound':False,'no_MATLAB':True,
         'code_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),
         'registration_sha256':hashlib.sha256(here.with_name('registration.json').read_bytes()).hexdigest(),
         'actual_process_status':'finished; stop after this diagnostic'}
    here.with_name('parity_counts_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='results'},indent=2))


if __name__=='__main__':
    main()
