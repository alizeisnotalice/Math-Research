#!/usr/bin/env python3
"""Exact comparison of three outer sufficient tests at the existing kernel."""
import argparse
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser()
ap.add_argument('--core',type=Path)
a=ap.parse_args()
p=a.core or ROOT.parent/'round008'/'certify.py'
if not p.exists():
    p=ROOT.parent.parent/'Math-Research-lowdim-agent'/'cube-lowdim-bounds'/'experiments'/'round008'/'certify.py'
spec=importlib.util.spec_from_file_location('outer_comparison_core',p)
c=importlib.util.module_from_spec(spec);sys.modules[spec.name]=c;spec.loader.exec_module(c)
R=Q(379,200);e2=Q(3,4);e4=Q(21,100);c.R=R
H2=8-8*c.log_point(2);H4=Q(136,3)-64*c.log_point(2)
old_base=Q(4)/(R*(R-1))
D=2*(2*c.tail_series(c.I.point(1))+c.I.point(1)/R)/(R-1)
# Independent closed form of the same diagonal tail, with exact log intervals.
closed=2*(-2*R*c.log_point(1-1/R)-2+1/R)/(R-1)
assert max(D.lo,closed.lo)<=min(D.hi,closed.hi)
slacks={'old_common_constant':old_base-(e2+e4)*H2,
        'separate_powers_old_base':old_base-e2*H2-e4*H4,
        'separate_powers_diagonal_base':D-e2*H2-e4*H4}
assert all(v.lo>0 for v in slacks.values())
assert (e4*(H2-H4)).lo>0
assert (D-old_base).lo>0
print(json.dumps({'status':'exact_outer_comparison_verified',
 'parameters':{'R':str(R),'e2':str(e2),'e4':str(e4)},
 'diagonal_base':D.pair(),'closed_form_overlap':True,
 'slacks':{k:{'interval':v.pair(),'lower_display':float(v.lo)} for k,v in slacks.items()},
 'unchanged_cost':str(R*R/(R-1)-e2/6-e4/15),
 'scope':'outer sufficient tests for existing parameters only; no new weak-type bound'},indent=2))
