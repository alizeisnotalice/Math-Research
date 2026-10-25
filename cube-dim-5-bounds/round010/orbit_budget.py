"""New orbit-corrected I expression via exact endpoints; no Gamma/kernel replay."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,time
P=Path(__file__).resolve().parents[1];began=time.monotonic();iterations=0
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def interval_root(x,k,digits=100):
    global iterations
    den=10**digits;N=x.numerator*den**k//x.denominator
    lo=0;hi=1<<((N.bit_length()+k-1)//k)
    while lo+1<hi:
        mid=(lo+hi)//2;iterations+=1
        if mid**k<=N:lo=mid
        else:hi=mid
    a,b=F(lo,den),F(lo+1,den);assert a**k<=x<b**k
    return a,b
def pi_interval():
    def atan(q):
        a=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(80)),F(0))
        return a,a+F(1,161*q**161)
    a,b=atan(5),atan(239)
    return 16*a[0]-4*b[1],16*a[1]-4*b[0]
def display(x,digits=60):
    den=10**digits;n=x.numerator*den//x.denominator
    return f'{n//den}.{n%den:0{digits}d}'
source=P/'upper_alpha1_029/check_results.json';a=json.loads(source.read_text())
j=a['independent_512bit_checks']['J']['independent_interval'];jl,ju=F(j['lower_fraction']),F(j['upper_fraction'])
old=a['independent_512bit_checks']['I']['independent_interval'];oldlo,oldhi=F(old['lower_fraction']),F(old['upper_fraction'])
alpha_l,alpha_u=interval_root(F(1,2),4);dl,du=2*(1-alpha_u)/21,2*(1-alpha_l)/21
pl,pu=pi_interval();assert F(3)<pl<pu<F(22,7)
jnew_l,jnew_u=jl-du,ju-dl;assert jnew_l>1
# Fourth power removes every fractional power from the full I expression.
lo4=jnew_l**5/(60**4*32**5*pu**7)
hi4=jnew_u**5/(60**4*32**5*pl**7)
ilo=interval_root(lo4,4)[0];ihi=interval_root(hi4,4)[1]
assert ihi<oldlo
cap_den=10**12;cap=F((ihi.numerator*cap_den+ihi.denominator-1)//ihi.denominator,cap_den);assert ihi<cap
point_path=P/'astra_green019/point_0.json';point=json.loads(point_path.read_text());gu=F(point['G_upper_dyadic'])
assert point['point']==['1','1','1','2/3','1/16'] and point['pair_rank_counts']=={'3':120}
assert all(F(v)==0 for v in point['principal_coefficients_inside_common_negative_factor'].values())
gate=interval_root(F(32**4,3125)*gu**5,4)
reqg=interval_root(F(3125,32**4)*cap**4,5)
cost_floor=interval_root(3125*ilo**4/gu**5,4)
assert cost_floor[0]>32 and gate[1]<ilo
files=[Path(__file__),P/'round010/registration.md',source,P/'astra_orbit038/report.md',
       point_path,P/'audit_zero041/report.md',P/'audit_zero041/result.json']
assert iterations<5000 and time.monotonic()-began<1
out={'status':'PASS_exact_new_auxiliary_I_expression_arithmetic_038_independent_theory_qualification_separate',
     'formula':'I_orbit=c1^(5/4)*(pi^2/60)*(J_flag-2*(1-2^(-1/4))/21)^(5/4), c1=1/(32*pi^3)',
     'fourth_power_formula':'I_orbit^4=(J_flag-delta)^5/(60^4*32^5*pi^7)',
     'J_flag_full_source_interval':[str(jl),str(ju)],'delta_interval':[str(dl),str(du)],
     'pi_machin80_interval':[str(pl),str(pu)],'new_I_budget_expression_interval':[str(ilo),str(ihi)],
     'new_I_strict_simple_rational_cap':str(cap),'new_I_upper_display':display(ihi),
     'old_complete_I_expression_interval':[str(oldlo),str(oldhi)],
     'strict_complete_old_lower_minus_new_upper':str(oldlo-ihi),
     'essential_g_upper_from_independently_checked_zero_ray_point':str(gu),
     'necessary_I_for_cost32_interval':[str(z) for z in gate],
     'necessary_I_threshold_upper_display':display(gate[1]),
     'g_required_using_new_simple_I_cap_interval':[str(z) for z in reqg],
     'g_required_using_new_simple_I_cap_upper_display':display(reqg[1]),
     'best_possible_cost_using_new_expression_lower_interval':[str(z) for z in cost_floor],
     'best_possible_cost_lower_display':display(cost_floor[0]),
     'new_C5_upper':False,'formal_C5_upper':'32','continuation_upper_difference_exact':'0',
     'limitations':['I interval encloses a proven upper-budget expression, not actual I',
                    'Point is an upper obstacle for legal g, not a global g lower bound',
                    '038 theorem qualification audited separately; no C5 breakthrough'],
     'source_sha256':{str(p.relative_to(P)):sha(p) for p in files},
     'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began,'integer_root_iterations':iterations,
                'arithmetic':'Python Fraction/int only','Gamma_or_Green_evaluations':0,'new_samples':0}}
(P/'round010/orbit_budget.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','new_I_upper_display','new_I_strict_simple_rational_cap',
    'necessary_I_threshold_upper_display','g_required_using_new_simple_I_cap_upper_display',
    'best_possible_cost_lower_display','runtime']},indent=2))
