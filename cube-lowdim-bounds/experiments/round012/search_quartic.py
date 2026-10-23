#!/usr/bin/env python3
"""Round12 bounded true-2D LP study adding psi4. Numerical discovery only.

Exactly 3 LPs: coarse quartic, enriched quartic, enriched e4=0 control.
7500 distinct joint points; four generators per point = 30000 evaluations.
No p4 angular maximum or continuous feasibility assumption is made.
"""
from __future__ import annotations
from fractions import Fraction as Q
import json
import math
from pathlib import Path
import platform
import sys
import time
import numpy as np
import scipy
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parent
RQ=Q(379,200);R=float(RQ)
INCUMBENT=Q(6938554331,1790000000)
LP_BUDGET=3
GENERATOR_EVAL_BUDGET=30000
COUNTS={'joint_points':0,'AH_R':0,'Apsi2':0,'Apsi3':0,'Apsi4':0}


def exterior_slice_series(p,t):
    """For t>2, h_p=2 int_(−1)^1 (1−|y|)^p/(t−y)^2 dy.
    Even moments give term_0=4/((p+1)t²) and a positive convergent series.
    This is a numerical stability device, not an interval certificate.
    """
    term=4./((p+1)*t*t);total=term
    for m in range(32):
        term*=((2*m+3)*(2*m+2))/((2*m+p+3)*(2*m+p+2)*t*t)
        total+=term
    return total

def hp(p,t):
    if t<=0:raise ValueError('Positive x/a is required.')
    if t>2:return exterior_slice_series(p,t)
    if p==2:
        if t==1:return 4*(2-2*math.log(2))
        return 4*(2-(1-t)*math.log(abs(1-t))-(1+t)*math.log1p(t)+2*math.log(t))
    if t==1:
        minus=0.;plus=-math.log(2.)
    elif t<1:
        minus=math.log(t)-math.log1p(-t)
        plus=math.log(t)-math.log1p(t)
    else:
        minus=-math.log1p(-1/t)
        plus=-math.log1p(1/t)
    if p==3:return 6*(3+(1-t)**2*minus+(1+t)**2*plus)
    if p==4:return 8*(2*t*t+11/3+(1+t)**3*plus+(1-t)**3*minus)
    raise ValueError(p)

def bump_generator(p,u,v):
    def term(a,x):return 0. if a<=0 else a**(p-1)*hp(p,x/a)
    return term(1-v,u)+term(1-u,v)

def tail_f(c):
    y=c/R
    if abs(y)<.25:
        # Same convergent series as the round008 core, here in floating point.
        total=0.;power=1.
        for n in range(40):
            total+=power/((n+2)*R);power*=y
        return total
    return -(R/c**2)*math.log1p(-c/R)-1/c

def tail_generator(r,d):
    return 2*(2*tail_f(r)+tail_f(d)+tail_f(-d))/(R-1)

def evaluate(points):
    rows=[];rhs=[]
    for r,q in points:
        u=r*(1+q)/2;v=r*(1-q)/2
        rows.append([bump_generator(p,u,v) for p in (2,3,4)])
        rhs.append(tail_generator(r,r*q))
    count=len(points);COUNTS['joint_points']+=count
    for key in ('AH_R','Apsi2','Apsi3','Apsi4'):COUNTS[key]+=count
    if sum(COUNTS[k] for k in ('AH_R','Apsi2','Apsi3','Apsi4'))>GENERATOR_EVAL_BUDGET:
        raise AssertionError('Generator budget exceeded')
    return np.asarray(rows),np.asarray(rhs)

def coarse_points():
    radii=np.r_[np.linspace(.005,.995,100),np.linspace(1.002,R-.002,50)]
    angles=np.arange(30)/30.
    return [(float(r),float(q)) for r in radii for q in angles],radii

def enrichment_points(radii):
    # Fifteen angle midpoints spanning [0,1), plus three increasingly near axes.
    qextra=np.r_[(2*np.arange(15)+.5)/30.,.99,.9999,1-1e-8]
    points=[(float(r),float(q)) for r in radii for q in qextra]
    rextra=[.998,.999,.9999,.99999,1-1e-8,1.,1+1e-8,1.00001,1.0001,1.001]
    qboundary=np.r_[np.arange(27)/27.,.99,.9999,1-1e-8]
    points += [(float(r),float(q)) for r in rextra for q in qboundary]
    assert len(points)==3000
    return points

def globals_rows():
    b=2-2*math.log(2)
    return np.array([[R*(R-1)*b]*3,[(R-1)/4]*3]),np.array([1.,R])

def active_rows(A,b,points,x,dual):
    slack=b-A@x
    descriptions=([{'kind':'true_2d_generator','r':r,'d_over_r':q,
                     'u':r*(1+q)/2,'v':r*(1-q)/2} for r,q in points]
                  +[{'kind':'global_exterior_sufficient'},{'kind':'global_kernel_dominance_sufficient'}])
    order=np.argsort(slack)[:16]
    return [{**descriptions[int(i)],'slack':float(slack[i]),
             'max_benefit_dual':float(-dual[i]),'coefficients':A[i].tolist(),'rhs':float(b[i])}
            for i in order]

def solve(name,rows,rhs,points,control=False):
    extra,extra_rhs=globals_rows()
    A=np.vstack([rows,extra]);b=np.r_[rhs,extra_rhs]
    c=np.array([-1/6,-1/10,-1/15])
    bounds=[(0,None),(0,None),(0,0) if control else (0,None)]
    result=linprog(c,A_ub=A,b_ub=b,bounds=bounds,method='highs',
                   options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    out={'name':name,'e4_forced_zero':control,'status':'finite_grid_numerical_only',
         'joint_points':len(points),'total_constraints':len(b),'success':bool(result.success),
         'solver_status':int(result.status),'solver_message':result.message,'iterations':int(result.nit)}
    if not result.success:return out
    x=result.x;dual=np.asarray(result.ineqlin.marginals)
    cost=R*R/(R-1)+float(c@x)
    out.update({'e2':float(x[0]),'e3':float(x[1]),'e4':float(x[2]),
                'cost_candidate':cost,'strict_incumbent':str(INCUMBENT),
                'candidate_minus_incumbent':cost-float(INCUMBENT),
                'minimum_grid_generator_slack':float(np.min(rhs-rows@x)),
                'global_sufficient_slacks':(extra_rhs-extra@x).tolist(),
                'max_benefit_primal_objective':float(-result.fun),
                'max_benefit_dual_objective':float(-b@dual),
                'numeric_dual_gap':float(result.fun-b@dual),
                'active_constraints':active_rows(A,b,points,x,dual),
                'lower_bound_reduced_costs':np.asarray(result.lower.marginals).tolist(),
                'upper_bound_reduced_costs':np.asarray(result.upper.marginals).tolist()})
    return out

def violations(rows,rhs,points,x):
    slack=rhs-rows@x
    order=np.argsort(slack)
    selected=[int(i) for i in order if slack[i]<-1e-8]
    return {'violation_count_below_minus_1e_minus_8':len(selected),
            'minimum_slack':float(slack[order[0]]),
            'worst_points':[{'r':points[i][0],'d_over_r':points[i][1],
                             'slack':float(slack[i])} for i in selected[:12]]}

def slice_argument_summary(points):
    # Geometry-only bookkeeping; this makes no additional generator calls.
    arguments=[]
    for r,q in points:
        u=r*(1+q)/2;v=r*(1-q)/2
        for a,x in ((1-v,u),(1-u,v)):
            if a>0:arguments.append(x/a)
    exterior=[t for t in arguments if t>2]
    return {'all_positive_a_slice_arguments_min':min(arguments),
            'all_positive_a_slice_arguments_max':max(arguments),
            't_gt_2_slice_count_per_bump':len(exterior),
            't_gt_2_min':min(exterior),'t_gt_2_max':max(exterior),
            'exterior_moment_series_terms':33,
            'exterior_moment_series_indices':[0,32],
            'exterior_remainder_exact_arithmetic_upper_for_each_slice':
                {str(p):str(Q(4,3*(p+1))*Q(1,2**66)) for p in (2,3,4)},
            'remainder_scope':'Analytic omitted-series bound only, excluding floating-point roundoff; no interval certificate.',
            'tail_f_small_argument_series_terms':40,'tail_f_series_condition':'abs(c/R)<1/4'}

def save(data):
    data['evaluation_counters']=dict(COUNTS)
    data['generator_evaluations_total']=sum(COUNTS[k] for k in ('AH_R','Apsi2','Apsi3','Apsi4'))
    (ROOT/'search_result.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')

def main():
    if (ROOT/'search_result.json').exists():
        raise SystemExit('Recorded round12 search already exists; refusing duplicate budget use.')
    start=time.monotonic()
    out={'round':12,'status':'running_finite_numerical_search_only','R':str(RQ),
         'incumbent_strict_upper':str(INCUMBENT),'LP_budget':LP_BUDGET,
         'generator_evaluation_budget':GENERATOR_EVAL_BUDGET,
         'generator_counting':'each true-2D point costs one AH_R and one Apsi_p for p=2,3,4',
         'family':'nonnegative e2,e3,e4; psi_p=(1-r)_+^p; fixed R',
         'formula_source':'slice_derivation.txt; check_slice.py and slice_check.json are independent inputs, not overwritten',
         'grid':{'coarse':'100 inner radii .005..995, 50 outer 1.002..R-.002, 30 q angles j/30',
                 'enrichment':'150 original radii x 18 added angles + 10 near-r1 radii x 30 angles',
                 'r1_both_sides':[.998,.999,.9999,.99999,1-1e-8,1.,1+1e-8,1.00001,1.0001,1.001],
                 'near_axis_q':[.99,.9999,1-1e-8]},
         'environment':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__},
         'LPs':[]}
    points,radii=coarse_points();rows,rhs=evaluate(points)
    coarse=solve('coarse_quartic',rows,rhs,points);out['LPs'].append(coarse);save(out)
    print(json.dumps({k:v for k,v in coarse.items() if k!='active_constraints'},indent=2),flush=True)
    if not coarse['success']:raise SystemExit('Coarse LP failed; no expansion.')
    added=enrichment_points(radii);added_rows,added_rhs=evaluate(added)
    assert len(set(points+added))==7500
    out['coarse_candidate_on_new_points']=violations(added_rows,added_rhs,added,
                                                     np.array([coarse[k] for k in ('e2','e3','e4')]))
    dense_points=points+added;dense_rows=np.vstack([rows,added_rows]);dense_rhs=np.r_[rhs,added_rhs]
    dense=solve('enriched_quartic',dense_rows,dense_rhs,dense_points);out['LPs'].append(dense);save(out)
    print(json.dumps({k:v for k,v in dense.items() if k!='active_constraints'},indent=2),flush=True)
    control=solve('enriched_e4_zero_control',dense_rows,dense_rhs,dense_points,True)
    out['LPs'].append(control)
    if dense['success'] and control['success']:
        out['shape_comparison']={'enriched_e4':dense['e4'],
             'quartic_cost_minus_e4zero_control':dense['cost_candidate']-control['cost_candidate'],
             'quartic_new_shape_selected':dense['e4']>1e-9,
             'scope':'Only the finite enriched LP; no continuous feasibility/optimality claim.'}
    out['slice_argument_and_stability_summary']=slice_argument_summary(dense_points)
    out['status']='finite_search_complete_numerical_only'
    out['LPs_completed']=len(out['LPs']);out['elapsed_seconds']=time.monotonic()-start
    save(out)
    print(json.dumps({k:v for k,v in control.items() if k!='active_constraints'},indent=2),flush=True)
    print(json.dumps({'status':out['status'],'evaluation_counters':COUNTS,
                      'generator_evaluations_total':out['generator_evaluations_total'],
                      'shape_comparison':out.get('shape_comparison'),
                      'coarse_candidate_on_new_points':out['coarse_candidate_on_new_points'],
                      'elapsed_seconds':out['elapsed_seconds']},indent=2),flush=True)
if __name__=='__main__':main()
