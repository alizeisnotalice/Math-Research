#!/usr/bin/env python3
"""Read-only one-point candidate reevaluation from the locked existing cache.

This utility was not invoked during the frozen experiment. Each future
invocation is one extra full point computation, and must be budgeted separately.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import argparse,json,runpy
from pathlib import Path
from fractions import Fraction
import numpy as np

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('coordinates',nargs=7,help='exact decimal/rational strings or floathex strings'); a=ap.parse_args()
    x=np.array([float.fromhex(v) if '0x' in v.lower() else float(Fraction(v)) for v in a.coordinates])
    if x[0]!=1 or np.any(x[1:]<0) or np.any(x[1:]>1): raise ValueError('Require x=(1,y), y in [0,1]^6')
    mod=runpy.run_path(str(Path(__file__).with_name('probe.py')))
    g,gr,c,d=mod['evaluate'](x)
    print(json.dumps({'kind':'one read-only full-point numerical candidate evaluation, no shape draw','full_point_calls':1,'floathex':[v.hex() for v in x],'Ghat':g,'gradient_free':gr[1:].tolist(),'blocks':c.reshape(8,1024).mean(axis=1).tolist(),'empirical_SEM':float(c.std(ddof=1)/np.sqrt(8192)),'max_single_sample':float(c.max()),'max_single_share':float(c.max()/c.sum()),'diagnostics':d,'strict_global_g':False,'strict_global_minimum':False,'new_C7_bound':False},indent=2))
if __name__=='__main__': main()
