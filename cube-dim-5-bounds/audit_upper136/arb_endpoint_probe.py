#!/usr/bin/env python3
"""Two new artificial rationals test only the narrow Arb endpoint API path."""
from datetime import datetime, timezone
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import os
import time
from flint import arb, ctx, __version__ as flint_version

HERE=Path(__file__).resolve().parent

def write(name,data):
    with (HERE/name).open('x') as stream:
        json.dump(data,stream,indent=2,sort_keys=True)
        stream.write('\n')

def main():
    assert not (HERE/'arb_probe_results.json').exists()
    codehash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    write('arb_probe_registration.json',{
        'task':'audit_upper136_narrow_Arb_endpoint_probe',
        'utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),
        'code_sha256':codehash,'precision_bits':384,
        'new_artificial_rationals':['31/97','97/31'],
        'old133_134_135_toys_repeated':0,'physical_endpoint_cases':0,
        'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,
    })
    begin=time.monotonic();ctx.prec=384;ctx.threads=1
    out=[]
    for text in ('31/97','97/31'):
        value=Q(text)
        root=(arb(value.numerator)/arb(value.denominator)).sqrt().sqrt()
        assert root.is_finite()
        lo=Q(str(root.lower().fmpq()));hi=Q(str(root.upper().fmpq()))
        assert 0<lo<hi and lo**4<=value<=hi**4
        out.append({'rational':text,'root_interval':[str(lo),str(hi)],
            'exact_fourth_power_residuals':[str(value-lo**4),str(hi**4-value)]})
    write('arb_probe_results.json',{
        'task':'audit_upper136_narrow_Arb_endpoint_probe','status':'PASS_two_new_artificial_rationals',
        'pid':os.getpid(),'seconds':time.monotonic()-begin,'flint_version':flint_version,
        'precision_bits':ctx.prec,'results':out,
        'physical_endpoint_cases':0,'old133_134_135_toys_repeated':0,
        'API_scope':'exact rational construction; sqrt twice; finite test; lower/upper to exact fmpq',
        'general_interval_API_or_kernel_validation_claim':False,
    })
    print('PASS: two new artificial rationals; narrow Arb sqrt/endpoints API only')

if __name__=='__main__':main()
