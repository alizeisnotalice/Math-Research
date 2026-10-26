#!/usr/bin/env python3
"""One parent-authorized point-freeze only, before the MC input is frozen."""
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
START=time.monotonic()

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def canonical(values):
    vals=tuple(sorted((abs(Q(str(t))) for t in values),reverse=True))
    return vals if len(vals)==6 and vals[0]==1 and all(t<=1 for t in vals) else None

def walk(node, found):
    if isinstance(node,dict):
        for key in ('numerators','point_numerators','quantized_numerators','final_numerators','starting_numerators'):
            den=node.get('denominator',node.get('starting_denominator'))
            nums=node.get(key)
            if den and isinstance(nums,list) and len(nums)==6:
                point=canonical([Q(int(n),int(den)) for n in nums])
                if point is not None: found.add(point)
        for value in node.values(): walk(value,found)
    elif isinstance(node,list):
        if len(node)==6 and all(isinstance(v,(str,int,float)) and not isinstance(v,bool) for v in node):
            try: point=canonical(node)
            except (ValueError,ZeroDivisionError): point=None
            if point is not None: found.add(point)
        for value in node: walk(value,found)

def main():
    from scipy.stats import qmc
    import scipy
    import numpy as np
    spec=json.loads((HERE/'point_freeze_inputs.json').read_text())
    assert spec['candidate_count']==32 and spec['take_first_eligible']==16 and spec['denominator']==4096
    old=set(); sources=[]
    for relative in spec['old_tables']:
        path=PROJECT/relative
        if path.suffix=='.txt':
            lines=path.read_text().splitlines(); den,count=map(int,lines[0].split())
            assert len(lines[1:])==count
            for line in lines[1:]:
                point=canonical([Q(int(v),den) for v in line.split()])
                if point is not None: old.add(point)
        else: walk(json.loads(path.read_text()),old)
        sources.append({'path':relative,'sha256':sha(path)})
    raw=qmc.Sobol(d=5,scramble=True,seed=20261027).random_base2(5)
    candidates=np.floor(raw*4096).astype(np.int64)
    points=[]; rejected=[]; seen=set(old)
    for i,row in enumerate(candidates):
        nums=[4096]+sorted(map(int,row),reverse=True)
        if not all(nums[j]>nums[j+1]>0 for j in range(5)):
            rejected.append({'candidate':i,'reason':'not_strict_positive_order'}); continue
        point=canonical([Q(v,4096) for v in nums])
        if point in seen:
            rejected.append({'candidate':i,'reason':'canonical_duplicate'}); continue
        seen.add(point)
        gap=min(nums[j]-nums[j+1] for j in range(5))
        points.append({'id':len(points),'candidate_index':i,'numerators':nums,'denominator':4096,
                       'min_adjacent_gap':str(Q(gap,4096)),
                       'distance_infinity_to_every_S2_pairplane_at_least':str(Q(gap,4096))})
        if len(points)==16: break
    assert len(points)==16
    payload={'status':'FROZEN_16_NEW_ORDERED_DYADIC_POINTS','points':points,'old_canonical_points':len(old),
             'limited_old_tables':sources,'rejected_before_selection':rejected,
             'source_sha256':sha(Path(__file__).resolve()),'input_sha256':sha(HERE/'point_freeze_inputs.json'),
             'scipy_version':scipy.__version__,'numpy_version':np.__version__,
             'no_kernel_or_MC_calls':True,'g_certificate':False}
    (HERE/'points.json').write_text(json.dumps(payload,indent=2)+'\n')
    receipt={'status':'PASS','exit_code':0,'wall_seconds':time.monotonic()-START,'children_created':0,
             'points_sha256':sha(HERE/'points.json'),'point_count':16,'MC_samples':0,
             'source_sha256':payload['source_sha256'],'input_sha256':payload['input_sha256']}
    (HERE/'point_freeze_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
    return 0

if __name__=='__main__':
    sys.exit(main())
