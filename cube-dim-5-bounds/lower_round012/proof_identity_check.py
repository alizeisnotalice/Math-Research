#!/usr/bin/env python3
"""Actual symbolic validation of the two scalar residual identities."""
import hashlib,json,sys,time
from pathlib import Path
import sympy as s
p=Path(__file__);started=time.monotonic();R,L,d=s.symbols('R L d',positive=True)
one=L-L/R*(R-d);two=(R+L)/2-d-L/R*(R-d)
assert s.expand(one-L*d/R)==0
assert s.simplify(two-(R-L)*(R-2*d)/(2*R))==0
out={'status':'PASS_exact_symbolic_identities','runtime':sys.executable,'sympy_version':s.__version__,'residual_1':str(s.factor(one)),'residual_2':str(s.factor(two)),'sign_hypotheses':['R>0','0<L<=R','0<=d<=R/2'],'proof_scope':'symbolic identities only; hypotheses and volume argument must be checked mathematically','wall_seconds':time.monotonic()-started,'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
p.with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
