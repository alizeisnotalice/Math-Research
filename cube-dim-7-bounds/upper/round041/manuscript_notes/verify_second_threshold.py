"""Exact manuscript threshold gates, adapted from accepted R042 source; no production numerical search or n=7 work."""
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

checks['direct_144_moment_loss_s8']=F(1,16)+F(1,8**5)+F(77,200)<1
assert all(checks.values())
receipt={
 'status':'PASS',
 'threshold_s':'2^36','threshold_log_n':'2^108',
 'checks':checks,'check_count':len(checks),
 'tail_scope':'Analytic monotone rational/power and Taylor bounds in manuscript; no finite sample inference.',
 'scope':'Only exact algebra at the displayed bases, including the direct C2=144 moment compensation; no n=7, parameter search or production computation.',
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
(root/'second_threshold_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
