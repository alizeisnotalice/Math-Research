#!/usr/bin/env python3
from pathlib import Path
from functools import lru_cache
import json
import satellite_screen as S
HERE=Path(__file__).resolve().parent
EPS=1/1024
S.A=S.A*(1+EPS)**.2;S.HALF=S.A/2
_,_,B,_=S.boxes_for(0,0,EPS)
@lru_cache(None)
def vol(boxes):
 if not boxes:return 0.
 cuts=sorted({x for b in boxes for x in b[:2]})
 if len(boxes[0])==2:return sum(hi-lo for lo,hi in zip(cuts,cuts[1:]) if any(b[0]<=lo and hi<=b[1] for b in boxes))
 return sum((hi-lo)*vol(tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]}))) for lo,hi in zip(cuts,cuts[1:]))
r=float(32*vol(B)/(1+EPS));c=json.loads((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_text())['display_only'][0]
out={'status':'floating_scaling_calibration','collapsed_orbits':'coincident primary and satellite aggregated to mass1025/1024 at each site','ratio_float':r,'C_display':c,'absolute_difference':abs(r-c),'formula':'common dilation factor=(1025/1024)^(1/5) preserves weak ratio exactly by change of variables','passed_tolerance_1e_12':abs(r-c)<1e-12,'initial_execution_note':'Initial inline run completed arithmetic but JSON serialization rejected NumPy bool; standalone script converts ratio to builtin float.'}
assert out['passed_tolerance_1e_12']
(HERE/'scaling_calibration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
