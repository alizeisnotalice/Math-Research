#!/usr/bin/env python3
"""Seal two exact saved-data obstacles; no special functions or optimizer."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import gzip,hashlib,json,os,sys,time
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
S=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes')
T=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round128/true_beta1_joint_moment')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dec(q):
    with localcontext()as c:c.prec=24;return str(Decimal(q.numerator)/Decimal(q.denominator))
def main():
    require_nice=os.getpriority(os.PRIO_PROCESS,0);assert require_nice>=10
    start=time.monotonic();r=read(S/'results.json');cones=read(S/'cone_table.json')['cones'];up=read(T/'upstream.json');cand=read(T/'candidate.json');cert=read(T/'candidate_certificate.json')
    lo,U=map(F,r['moment_interval']);fU=F(r['pi_factor_interval'][1]);beta=F(r['old_beta2_interval'][1]);gap=U-beta
    assert beta==F(up['beta_by_h'][5][1])
    assert U-lo==sum(F(r[k])for k in ['weighted_box_width','tail_upper','pi_factor_width'])
    widths=[]
    for cone,row in zip(cones,r['leaf_files']):
        path=S/row['path'];assert sha(path)==row['sha256']
        with gzip.open(path,'rt')as f:data=json.load(f)
        for leaf in data['leaves']:
            a,b=map(F,leaf['integral_interval']);widths.append(fU*F(cone['alpha'])*(b-a))
    assert len(widths)==5376 and sum(widths)==F(r['weighted_box_width'])
    widths.sort(reverse=True);s15=sum(widths[:15]);s2188=sum(widths[:2188]);s2189=s2188+widths[2188]
    assert s15<=s2188<=gap<s2189
    x=F(cand['numerators'][32],int(cand['denominator']));assert x==F(2514712686865,2**44)
    HU5=F(cert['M_power6_upper']);HL5=F(cert['prefix_only_H_power5_lower']);theta=F(cert['threshold_M_power6_exact'])
    assert 0<HL5<=HU5 and theta==F(8,9)**6*24576**5/65000**6
    Zlo=1-x+x*lo/beta
    assert Zlo*HL5>theta
    reqU=beta*(1-(1-theta/HU5)/x);reqL=beta*(1-(1-theta/HL5)/x)
    assert lo>reqL>=reqU>0
    out={'status':'PASS_TWO_EXACT_SAVED_DATA_OBSTACLES','precision':'stdlib Fraction only; Decimal display only','nice':require_nice,
      'patch_obstacle':{'scope':'replace at most2188 original130 leaf integral upper bounds, keeping all other130 leaves, tail and pi factor unchanged; even exact integral evaluation cannot improve old beta2 upper',
        'old_upper_minus_beta2_exact':str(gap),'top15_recoverable_width_exact':str(s15),'top2188_recoverable_width_exact':str(s2188),'top2189_recoverable_width_exact':str(s2189),
        'top15_cannot_improve_beta2':True,'necessary_minimum_patched_old_leaves':2189},
      'fixed128_coefficient_ray_obstacle':{'scope':'only the existing128 unnormalized coefficient ray with only h1/128 beta2 changed and all coefficients jointly renormalized; lower proved with saved prefix-H lower, not saved M upper',
        'weight_x_exact':str(x),'Z_formula':'1-x+x*b/old_beta2','true_beta2_known_lower_exact':str(lo),
        'g0_beta2_threshold_from_saved_upper_display_only':dec(reqU),
        'g0_beta2_threshold_from_prefix_lower_display_only':dec(reqL),
        'best_optimistic_fixed_ray_objective_power6_lower_exact':str(Zlo*HL5),
        'best_optimistic_fixed_ray_lower_over_g0_display_only':dec(Zlo*HL5/theta),
        'strict_full_fixed_ray_g0_obstacle':True,'all_reoptimized_weights_excluded':False},
      'new_functional_kernel_erfc_exp_log_root_optimizer_calls':0,'actual_M6_lower_bound':False,'global_g0_proved':False,'new_C6_bound':False,
      'source_hashes':[{'path':str(p),'sha256':sha(p)}for p in [S/'results.json',S/'cone_table.json',T/'upstream.json',T/'candidate.json',T/'candidate_certificate.json']],
      'script_sha256':sha(Path(__file__)),'wall_inside_checker_seconds':time.monotonic()-start}
    (P/'obstacle_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'minimum_patch_leaves':2189,'fixed_ray_prefix_lower_over_g0':dec(Zlo*HL5/theta),'threshold_upper_display':dec(reqU),'threshold_prefix_lower_display':dec(reqL),'new_functional_calls':0}))
if __name__=='__main__':main()
