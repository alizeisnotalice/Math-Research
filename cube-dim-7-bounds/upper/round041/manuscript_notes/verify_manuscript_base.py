"""Exact Fraction base point for manuscript explicit-loglog-threshold subsection; adapted from frozen R036 verifier, without external absolute-path dependencies."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

root=Path(__file__).resolve().parent
z=F(128)
p=F(3,8)
powp=lambda m:p**m
D=12*powp(64)
Q=1-F(3,4)*D
Em=z*(3*z*powp(128)+F(3,4)*D)/(4*Q)
Er=(F(13,100)+48*D*powp(512)+4*D*powp(128))/(2*Q)
Ec=F(17,200)+D+64*D*powp(512)
V=(D+powp(64)/100+8*powp(640))/2
E=Em+Er+Ec
Perr=E+V+V/(1-V)*(z/4+E)
Bad=F(2000000,3)*z*z*12*powp(32)+F(400,121)*12*powp(32)+51200*powp(160)
Tau=7200*12*powp(32)
Rem=131072*((3*z+1)*powp(64)+(9*z*z+1)*powp(192))
margins={
 'd_lt_quarter':F(1,4)-D,
 'x_ge_3L_over4':F(1,4)-3*z*powp(128),
 'q_gt_half':Q-F(1,2),
 'h_zeta_lt_d_ratio':1-1/(100*z)-8*powp(576),
 'target_rho_gt_4':Q/(2*powp(128))-4,
 'corner_target_face_gap':1/(1+F(3,4)*D)-2*powp(128)-F(1,2),
 'target_cheb_below_quarter':F(1,4)-F(20000,3)*D**3,
 'source_Ld_gt_one':1/powp(64)-1,
 'target_nd_gt_two':1/(2*powp(192))-2,
 'max_loss_lt_one':1-V,
 'uniform_product_error_lt_two':2-Perr,
 'source_bad_below_one36':F(1,36)-Bad,
 'tau_below_half_row':1-Tau,
 'W_remainder_below_delta_half':1-Rem,
}
for name,margin in margins.items():
    assert margin>0,name
assert Perr<F(1,2)
assert Bad<F(1,100)
assert Tau<F(1,100)
assert Rem<F(1,100)
extra={
 'sqrt128_lt_12':128<12**2,
 'e_gt_8over3_Taylor':F(1)+1+F(1,2)+F(1,6)==F(8,3),
 'L_gt_64_at_base':1/powp(128)>64,
 'rq_mean_gap_coefficient':F(1,2)+F(9,64)<F(3,4),
 'band_h_over_d_lt_half':1/(100*z)<F(1,2),
 'post_corner_n6':F(3,4)**5<F(1,4),
 'source_corner_L_ge_4096z':128**3>24*4096,
 'source_corner_sqrtL_ge_16z':F(128**2,8)>=16*128,
 'source_corner_new_coefficient':F(25,8)<51200,
 'kappa_n6_x48':48*6**2>2,
 'energy_final_constant':81<128,
}
assert all(extra.values())

receipt={"status":"PASS","method":"Python Fraction exact arithmetic","checks":len(margins),"all_14_margins_positive":all(x>0 for x in margins.values()),"Perr_lt_half":Perr<F(1,2),"Bad_lt_1_over_100":Bad<F(1,100),"Tau_lt_1_over_100":Tau<F(1,100),"Rem_lt_1_over_100":Rem<F(1,100),"extra_checks":extra,"scope":"Only explicit threshold base-point algebra; no dimension-seven verification, optimization or production recomputation."}
(root/'exact_base_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False,indent=2))
