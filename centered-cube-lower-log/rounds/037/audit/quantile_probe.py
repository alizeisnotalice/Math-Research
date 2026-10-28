"""Exact binomial/true-box query diagnostic, not a Gaussian or tail MC."""
from pathlib import Path
from fractions import Fraction as F
import json,math

HERE=Path(__file__).resolve().parent
def partition(n,R):
    nums=[math.comb(n,k)for k in range(n+1)];den=1<<n
    cuts=[];s=0;j=0
    for b in range(1,R):
        while s*R<b*den:s+=nums[j];j+=1
        cuts.append(j-1)
    endpoints=[-1]+cuts+[n]
    weights=[F(sum(nums[endpoints[b]+1:endpoints[b+1]+1]),den)for b in range(R)]
    return endpoints,weights
rows=[]
geometry_checks=0
central_checks=0
for m in range(65):
    central=F(math.comb(2*m,m),4**m)
    assert central*central*(3*m+1)<=1
    assert F(2*m+1,2*m+2)**2<=F(3*m+1,3*m+4)
    central_checks+=1
for n,R in [(63,2),(63,4),(63,8),(255,2),(255,4),(255,8),(255,16),
            (1023,2),(1023,4),(1023,8),(1023,16),(1023,32)]:
    K=4;cuts,pi=partition(n,R)
    assert sum(pi,F(0))==1 and min(pi)>0
    b=R-1;delta=F(1,64*K*R);L=1+2*delta
    theta=F(2*b+1,2*R);atom=theta/K;query=(atom-F(1,2))%1
    # The macro plateau is simultaneously the all-minus b0 source support:
    # half-period translation maps each selected plus atom to table0 minus.
    for j in range(K//2):
        shifted=((j+theta)/K+F(1,2))%1
        table0_minus=(j+K//2+(1-F(1,2*R)))/K
        assert shifted==table0_minus
    clip=lambda z:min(F(1),max(F(0),z))
    for offset in [-delta/4,F(0),delta/4]:
        x=query+offset
        for table in range(R):
            th=F(2*table+1,2*R);parts=[F(0),F(0)]
            for cell in range(K):
                a=(cell+(th if cell<K//2 else 1-th))/K
                mass=sum((clip((x-a-sh+L/2+delta/2)/delta)-clip((x-a-sh-L/2+delta/2)/delta)
                          for sh in [-1,0,1]),F(0))/(K*L)
                parts[int(cell>=K//2)]+=mass
            assert parts==[(F(1,2)+(F(1,K)if table==b else 0))/L,F(1,2)/L]
            geometry_checks+=1
    # The extra interval has length 2delta and covers exactly one plus atom
    # from table b. On its plateau other tables have zero extra mass.
    total=(1+F(1,K))/L;p=(F(1,2)+F(1,K))/(1+F(1,K))
    assert p==F(3,5)
    threshold=cuts[-2]+1
    pb=F(sum(math.comb(n,s)*3**s*2**(n-s)for s in range(threshold,n+1)),5**n)
    likelihood=(1+F(1,K))**n
    # Divide all table likelihoods by the common factor L^-n.
    joint=likelihood*pb+1-pi[-1]
    prior=likelihood*pi[-1]+1-pi[-1]
    ratio=joint/prior
    sorted_pi=sorted(pi,reverse=True)
    prior_density_weak=max(F(k+1)*w for k,w in enumerate(sorted_pi))
    # For an aligned extra window of one complete coarse-cell length,
    # all tables instead have identical A and PB laws, forcing joint/prior=1.
    aligned_ratio=sum(pi,F(0))
    plateau_one_coord=F(K,2)*delta
    log_plateau_weak_upper=math.log(2)+n*(math.log(1+1/K)-math.log(L)+math.log(float(plateau_one_coord)))
    rows.append({'n':n,'K':K,'R':R,'cuts':cuts,'pi_exact':list(map(str,pi)),
                 'delta':str(delta),'L':str(L),'query_coordinate':str(query),
                 'p_selected':str(p),'selected_PB_mass':float(pb),
                 'pi_selected':float(pi[-1]),'PB_over_prior':float(pb/pi[-1]),
                 'true_joint_over_prior_query_ratio':float(ratio),
                 'selected_likelihood_prior_weight':float(likelihood*pi[-1]/prior),
                 'density_joint_weak_exact':'1',
                 'density_prior_weak_exact':str(prior_density_weak),
                 'density_prior_weak_display':float(prior_density_weak),
                 'full_scale_ratio_on_plateau_lower':float(1/max(pi)),
                 'full_scale_ratio_on_plateau_upper':float(1/pi[0]),
                 'all_tables_identical_query_ratio':str(aligned_ratio),
                 'plateau_one_coordinate_measure':str(plateau_one_coord),
                 'plateau_interval_width_per_atom':str(delta),
                 'log_upper_fixed_L_tail_contribution_on_plateau':log_plateau_weak_upper,
                 'full_scale_M_on_plateau_log':n*math.log(64*R),
                 'full_scale_weak_contribution_on_plateau_exact':str(F(1,2)**n),
                 'full_scale_weak_contribution_on_plateau_log':-n*math.log(2)})
out={'scope':'Exact Fraction binomial masses and actual positive-box centered-cube query; finite K4 source, not original matched source; no Gaussian or MC.',
     'exact_coordinate_box_geometry_checks':geometry_checks,
     'exact_central_binomial_induction_checks':central_checks,
     'quantile_convention':'cut b=min k with CDF(k)>=b/R; ties yield actual bin masses, never replaced by 1/R.',
     'rows':rows}
(HERE/'quantile_results.json').write_text(json.dumps(out,indent=2)+'\n')
for x in rows:print(x['n'],x['R'],'pi',x['pi_selected'],'PB',x['selected_PB_mass'],'query ratio',x['true_joint_over_prior_query_ratio'],'log fixed-L plateau',x['log_upper_fixed_L_tail_contribution_on_plateau'],'log full-M plateau',x['full_scale_weak_contribution_on_plateau_log'])
