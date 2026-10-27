#!/usr/bin/env python3
"""One fixed eight-small-integer root-boundary check; no actual node inputs."""
from pathlib import Path
import ast,json,os,time
from math import isqrt
from common import sha,read,put
HERE=Path(__file__).resolve().parent

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();a.output.mkdir(exist_ok=False)
    cases=[0,1,15,16,17,80,81,82]
    put(a.output/'registration.json',{'pid':os.getpid(),'cases':cases,'selfcheck_sha256':sha(__file__),
        'consumer_sha256':sha(HERE/'consumer.py'),'actual_node_root_gain_DFS_cases':0})
    start=time.monotonic();rows=[]
    for name in ['consumer.py','prepare.py','common.py','selfcheck.py']:ast.parse((HERE/name).read_text())
    tree=ast.parse((HERE/'consumer.py').read_text())
    imports=[n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
    assert set(imports)=={'pathlib','fractions','math','common'}
    for v in cases:
        r=isqrt(isqrt(v));assert r**4<=v<(r+1)**4;rows.append({'input':v,'root':r,'left_residual':v-r**4,'next_residual':(r+1)**4-v})
    assert 'skip_initial=count if resume else 0' in (HERE/'consumer.py').read_text()
    put(a.output/'results.json',{'pid':os.getpid(),'status':'syntax_import_scope_and8new_small_root_boundaries_pass',
        'seconds':time.monotonic()-start,'rows':rows,'actual_node_root_gain_DFS_cases':0,'old_symbolic_or_baseline_replay':False})
