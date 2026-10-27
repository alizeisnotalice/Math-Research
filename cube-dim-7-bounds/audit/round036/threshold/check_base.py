#!/usr/bin/env python3
"""Independent coarse rational base certificate only; no production run."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import platform
import sys

root=Path(__file__).resolve().parent
p=F(3,8)
assert p**8<F(1,2500)
T=F(1,2500**4)  # p^32<T; every p^m for m>=32 is also <T.
U=F(1,10**9)
D=F(1,10**6)
Q=F(99,100)
assert T<U
assert 12*p**64<D and 1-F(3,4)*D>Q
z=F(128)
Em=z*(3*z*U+F(3,4)*D)/(4*Q)
Er=(F(13,100)+48*D*U+4*D*U)/(2*Q)
Ec=F(17,200)+D+64*D*U
assert Em<F(1,1000) and Er<F(33,500) and Ec<F(43,500)
E=F(153,1000)
V=F(1,10**6)
assert Em+Er+Ec<E
assert (D+U/100+8*U)/2<V
Perr=E+V+V/(1-V)*(z/4+E)
Bad=(F(2000000,3)*128**2*12+40+51200)*T
Tau=86400*T
Rem=131072*(385+147457)*T
assert Perr<F(1,5)
assert Bad<F(1,250)<F(1,100)<F(1,36)
assert Tau<F(1,100)
assert Rem<F(1,1000)
domains={
 'd_lt_quarter':D<F(1,4),
 '3z_over_L_lt_quarter':3*z*U<F(1,4),
 'q0_gt_half':Q>F(1,2),
 'shift_over_d_lt_one':1/(100*z)+8*U<1,
 'source_rho_gt_four':Q/(2*U)>4,
 'target_corner_gap':1/(1+F(3,4)*D)-2*U>F(1,2),
 'target_cheb_lt_quarter':F(20000,3)*D**3<F(1,4),
 'Ld_gt_one':1/U>1,
 'L2d_over2_gt_two':1/(2*U)>2,
 'max_loss_lt_one':V<1,
}
assert all(domains.values())
assert 128**3>24*4096 and F(128**2,8)>=16*128
assert F(3,4)**5<F(1,4)
source=Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round036/threshold_theory')
frozen=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round035/energy_theory/PROOF.md')
sources={str(f):hashlib.sha256(f.read_bytes()).hexdigest()
         for f in [source/'THRESHOLD_PROOF.md',source/'verify_threshold.py',
                   source/'verification.json',frozen]}
assert sources[str(source/'THRESHOLD_PROOF.md')]=='064e3e6d697220817642250961180fcf2f747dc36446eb693344fb74f719bed1'
assert sources[str(frozen)]=='f7d98eb21b89d83264053fda2c594387387facfdcbdbec7cf83a41e5dd55f0bd'
receipt={
 'status':'PASS','python':sys.version,'platform':platform.platform(),
 'auditor_pid':os.getpid(),
 'arithmetic':'stdlib Fraction; coarse positive rational bounds only',
 'Perr_upper':str(Perr),'Bad_upper':str(Bad),'Tau_upper':str(Tau),'Rem_upper':str(Rem),
 'domain_checks':domains,'source_sha256':sources,
 'scope':'one exact base certificate plus analytical all-tail review; no Arb, boxes, search, n7, MATLAB or production replay',
}
(root/'base_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
