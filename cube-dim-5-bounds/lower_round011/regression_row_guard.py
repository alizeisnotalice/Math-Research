#!/usr/bin/env python3
"""Exact N=4 interior counterexample to round010's row interval."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'lower'))
from certify_rectangular import field
F,u,v,bounds,_,_=field(); a=1+u**3/2; s=(1+u)/2; h=a/8; r=u**2/2
center=(r/2,s,s/4,s/4,s/4)
old_box=(F(0),r,F(0),2*s,F(0),s/2,F(0),s/2,F(0),s/2)
assert all(old_box[2*q]<x<old_box[2*q+1] for q,x in enumerate(center))
atoms=[(h if k%2 else F(0),k*s,F(0),F(0),F(0)) for k in range(4)]
contained=[all(p-r<=x<=p+r for p,x in zip(atom,center)) for atom in atoms]
assert contained==[False,True,False,False]
# Correct y interval is empty because the rows span more than the diameter.
assert 3*s-r>r
assert 2*s>r
out={'status':'PASS_exact_regression_counterexample','old_source_tag':[0,0,[1,4,1,1,1]],'N':4,'common_radius':'u^2/2','center':['u^2/4','s','s/4','s/4','s/4'],'old_box':'[0,r] x [0,2s] x [0,s/2]^3','center_strictly_inside_old_box':True,'selected_atoms_contained':contained,'missing_row3_distance':'2*s > r','correct_y_interval':'[3*s-r,r] is empty','u_relation':'u^5=2','old_source_sha256':hashlib.sha256((ROOT/'lower_round010/two_row_staggered.py').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'regression_row_guard.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
