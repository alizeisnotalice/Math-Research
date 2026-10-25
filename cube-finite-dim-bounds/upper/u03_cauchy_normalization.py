"""U03 alpha=1 one-dimensional resolvent normalization, same two registered grids.

Jump density is |z|^-2, so A=-pi|D|. Standard Cauchy Green and lifetime
must both be divided by pi relative to A=-|D|.
"""
from pathlib import Path
import hashlib
import json
from math import pi,log,sqrt
import numpy as np
from scipy.linalg import inv
from u03_restricted_scan import generator

rows=[]
for ncells in (65,129):
    q,x,h,killing,defect=generator(1.,ncells)
    resolvent=inv(-q)
    source=ncells//2
    target=source+round(.5/h)
    a=float(x[target])
    lifetime=float(resolvent[source].sum())
    green=float(resolvent[source,target]/h)
    exact_green=log((1+sqrt(1-a*a))/a)/(pi*pi)
    rows.append({'N':ncells,'a':a,'lifetime':lifetime,'lifetime_exact':1/pi,
                 'lifetime_relative_error':lifetime*pi-1,
                 'green':green,'green_exact':exact_green,
                 'green_relative_error':green/exact_green-1})
print(json.dumps({'status':'normalization diagnostic, not error bounds',
                  'rows':rows,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
