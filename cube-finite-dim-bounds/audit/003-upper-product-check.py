#!/usr/bin/env python3
"""Independent exact finite identities for the direct CND audit proof."""
import hashlib,json
from pathlib import Path
import sympy as s
A,B,C=s.symbols('A B C',real=True)
rect=s.Matrix([[0,A,B,C],[A,0,C,B],[B,C,0,A],[C,B,A,0]])
vecs=[s.Matrix(v) for v in [(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1)]]
rect_forms=[s.expand((v.T*rect*v)[0]) for v in vecs]
assert rect_forms==[4*(A-B-C),4*(B-A-C),4*(C-A-B)]
tri=s.Matrix([[0,C,A],[C,0,B],[A,B,0]])
w=s.Matrix([1,1,-2]);tri_form=s.expand((w.T*tri*w)[0])
assert tri_form==2*C-4*A-4*B
qj,qk,gj,gk=s.symbols('qj qk gj gk',real=True)
fst=gj*qk+gk*qj
assert s.expand(fst-gj-(-gj*(1-qk)+gk*qj))==0
# Algebraic rearrangement from fst-gj >= -gk.
assert s.expand((fst-gj+gk)-(gk*(1+qj)-gj*(1-qk)))==0
result={'engine':'SymPy '+s.__version__,'scope':'exact finite identities only; cutoff/DCT are analytic proofs','rectangle_forms':list(map(str,rect_forms)),'three_point_form':str(tri_form),'product_identity':str(fst-gj),'status':'all exact identities PASS','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('003-upper-product-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
