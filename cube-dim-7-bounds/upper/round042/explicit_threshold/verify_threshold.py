"""Minimal exact threshold gates; no production numerical search."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math

root=Path(__file__).resolve().parent
S=2**36
fact8=math.factorial(8)
checks={
 'x_ge_L_over2_s2':1-F(1,2**2)>=F(1,2),
 'zeta_t_le_h_over2_s2':2**13>=1600,
 'target_source_geometry_s8':F(8**3,4)>4,
 'h_shift_lt_d_over2_s2':F(3,200*2)<F(1,2),
 'target_Cheb_at_S':F(20000,3*S**3)<F(1,4),
 'target_nd_gt2_s2':F(2**5,2)>2,
 'target_corner_gap_s4':F(16,19)-F(1,32)>F(1,2),
 'target_n_ge6_s2':F(2**6,2)>=6,
 'target_corner_n6':F(3,4)**5<F(1,4),
 'mean_error_s8':F(9,28)<F(1,3),
 'row_shift_error_s8':F(19,100)+F(4,8**4)<F(1,5),
 'row_log_error_s8':F(4,8**4)+F(4,8**3)<F(1,10),
 'row_total_error_lt2':F(1,3)+F(1,5)+F(1,10)<2,
 'source_corner_at_S':F(819200,S**3)<F(1,36),
 'raw_energy_cE':9*64<1024,
 'truncation_constants':2**20>=4*144*1024,
 'max_degree_constants':2**20>=36,
 'tau_Taylor8_integer':100*4**8*fact8<2**38,
 'tau_power_at_S':F(2**38,S**4)<F(1,2**12),
 'moment_L_gate':S**3>=64,
}
Rscalar=F(2,S)+F(2,S**5)+F(fact8,S**7)+F(2**9*fact8,S**5)
checks['scalar_remainder_at_S']=Rscalar<F(1,2**32)
assert all(checks.values()),[k for k,v in checks.items() if not v]
full=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round042/second_moment_cut/PROOF.md')
moment=root.parent/'second_moment_audit'/'PROOF.md'
fullhash=hashlib.sha256(full.read_bytes()).hexdigest()
momenthash=hashlib.sha256(moment.read_bytes()).hexdigest()
assert fullhash=='722764654cd15a8bf7fb18a799c298b54ced2f0188f5554ccae173240ddf90a2'
assert momenthash=='425a383952f341f42ff0dee9651924b00fff61e958c2957f355e8cc9fdc72d00'
receipt={
 'status':'PASS: all exact gates and scalar threshold',
 'threshold_s':'2^36', 'threshold_log_n':'2^108',
 'checks':checks,
 'tail_scope':'all tails bounded analytically by monotone rational functions and powers in THRESHOLD_PROOF.md',
 'frozen_full_proof_sha256':fullhash,
 'frozen_moment_sha256':momenthash,
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                  for name in ('THRESHOLD_PROOF.md','verify_threshold.py')},
 'no_Arb_search':True,'no_production_computation':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
