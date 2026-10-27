#!/usr/bin/env python3
"""Independent small rational threshold checks; no producer import or search."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib
import json
import os
import platform
import sys

S=2**36
cE=Q(1,2**10)
H=2**20
fac=factorial(8)
checks={
 'shift_at_s2':2**13>1600,
 'target_Cheb_at_S':Q(20000,3*S**3)<Q(1,4),
 'target_corner_at_s4':Q(16,19)-Q(1,32)>Q(1,2),
 'n6_corner_factor':Q(3,4)**5<Q(1,4),
 'row_mean_at_s8':Q(9,28)<Q(1,3),
 'row_tolerance_at_s8':Q(19,100)+Q(4,8**4)<Q(1,5),
 'row_log_at_s8':Q(4,8**4)+Q(4,8**3)<Q(1,10),
 'product_error_total':Q(1,3)+Q(1,5)+Q(1,10)<2,
 'source_corner_at_S':Q(819200,S**3)<Q(1,36),
 'raw_energy_gap':Q(1,9*64)>cE,
 'second_moment_cut':H>=4*144/cE,
 'degree_cut':H>36,
 'deficit_constant':cE/(2*H)==Q(1,2**31),
 'tau_Taylor_coefficient':100*4**8*fac<2**38,
 'tau_at_S':Q(2**38,S**4)<cE/4,
 'moment_explicit_gate':S**3>=64,
}
R=Q(2,S)+Q(2,S**5)+Q(fac,S**7)+Q(2**9*fac,S**5)
checks['scalar_at_S']=R<Q(1,2**32)
assert all(checks.values()),[k for k,v in checks.items() if not v]
base=Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042')
files=[base/'explicit_threshold/THRESHOLD_PROOF.md',base/'explicit_threshold/verify_threshold.py',base/'explicit_threshold/verification.json',base/'second_moment_audit/PROOF.md',Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round042/second_moment_cut/PROOF.md')]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
assert hashes[str(files[0])]=='a73bae47dbef7b56b95dcc28bf5dc368828265e91c0de49d03a3dcf518e8c0c5'
assert hashes[str(files[3])]=='425a383952f341f42ff0dee9651924b00fff61e958c2957f355e8cc9fdc72d00'
assert hashes[str(files[4])]=='722764654cd15a8bf7fb18a799c298b54ced2f0188f5554ccae173240ddf90a2'
receipt={'status':'PASS','threshold_s':'2^36','threshold_log_n':'2^108',
 'python':sys.version,'platform':platform.platform(),'auditor_pid':os.getpid(),
 'arithmetic':'stdlib Fraction and integer powers; no producer code imported',
 'checks':checks,'scalar_positive_remainder_upper_at_S':str(R),
 'source_sha256':hashes,
 'scope':'one exact base plus analytical uniform all-tail proof; no production run, Arb, n7 or search'}
out=Path(__file__).resolve().parent/'exact_receipt.json'
out.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
