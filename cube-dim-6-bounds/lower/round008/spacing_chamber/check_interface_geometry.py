#!/usr/bin/env python3
"""Both boundary checks imply the affine family has common-radius witnesses."""
from fractions import Fraction as Q
from pathlib import Path
from math import prod
import hashlib
import json
import os
import time
import chamber as C

HERE=Path(__file__).resolve().parent

def main():
    started=time.monotonic()
    path=HERE/'chamber_interface.json'
    packet=json.loads(path.read_text())
    values={name:(C.field(x['c0']),C.field(x['c1'])) for name,x in packet['endpoints'].items()}
    s0,s1=C.field(packet['spacing2']['c0']),C.field(packet['spacing2']['c1'])
    checks=0
    for t in (Q(0),Q(1)):
        a=(C.ONE+C.u**4/2,s0+s1*t,*((C.ONE+C.u)/2 for _ in range(4)))
        for box,ell in zip(packet['boxes'],packet['box_lengths']):
            r=C.root_radius(prod(ell))
            for axis,(ids,l,aj) in enumerate(zip(box,ell,a)):
                q=-(l//2) if l%2 else 1-l//2
                lo,hi=(values[name][0]+values[name][1]*t for name in ids)
                assert C.ZERO<=lo<=hi<=aj/2
                for x in (lo,hi):
                    for index in (q,q+l-1):
                        assert -r<=x-index*aj<=r
                        checks+=1
    out={'status':'PASS','pid':os.getpid(),'boxes':len(packet['boxes']),
         'boundary_parameters':['0','1'],'source_endpoint_checks':checks,
         'claim':'all center/source differences and cell faces are affine, so both boundary inequalities imply validity for all t in [0,1]',
         'common_radius':'r=(prod lengths)^(1/6)/2, independent of t and the same in all axes',
         'union_calls':0,'input_sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest(),
                                        str(HERE/'chamber.py'):hashlib.sha256((HERE/'chamber.py').read_bytes()).hexdigest()},
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'elapsed_seconds':time.monotonic()-started,'actual_process_status':'finished'}
    (HERE/'geometry_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
