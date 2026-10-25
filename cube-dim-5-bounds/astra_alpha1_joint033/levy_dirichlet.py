#!/usr/bin/env python3
"""Three fixed beta=1/2 rank-two quotient negative moments."""
import hashlib,itertools,json,math,os,sys,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.special import betaincinv,gammaln
from scipy.stats import qmc
from flint import arb,ctx

P=Path(__file__).resolve().parent
SOURCE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper-audit-cont/cube-dim-5-bounds/upper_rank2_interface')
START=time.monotonic()
S=F(247,2560)
RHO=F(247,256)

def digest(path):
    p=Path(path);return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}

def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator

def dump(name,data):
    tmp=P/(name+'.tmp');tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n');tmp.replace(P/name)

def geometry(row):
    m=row['independent_class_count'];lam=sp.Rational(row['normalization_lambda'])
    c=sp.Matrix(row['Q_offdiag_coefficient_matrix'])
    kappas=[]
    for cl in row['classes']:
        kappas.append(sp.simplify(sum(sp.sqrt(sp.Rational(a)) for a in cl['Xi_original_coefficients'].values())))
    assert all(k.is_Integer and k>0 for k in kappas)
    h=[sp.simplify(k*k/4) for k in kappas]
    terms=[]
    for i,j in itertools.combinations(range(m),2):
        comp=[k for k in range(m) if k not in (i,j)]
        coeff=sp.simplify(c[i,j]*h[i]*h[j]/lam)
        assert coeff.is_Rational and coeff>0
        terms.append({'complement':comp,'coefficient':str(coeff),'pair':[i,j]})
    # Verify rational quotient polynomial at all-one ORIGINAL S, not X=1.
    xs=[sum(sp.Rational(v) for v in cl['Xi_original_coefficients'].values()) for cl in row['classes']]
    assert sum(c[i,j]*xs[i]*xs[j] for i,j in itertools.combinations(range(m),2))==lam
    return {'triple':row['representative_J'],'Gram_det':row['Gram_det'],'theta':row['triple_theta'],
            'm':m,'Lambda':str(lam),'kappa_beta_half':[str(k) for k in kappas],
            'h_kappa_squared_over_four':[str(x) for x in h],'normalized_R_terms':terms}

def gamma_factor(m):
    s=float(S);a=.5+s
    return math.exp(m*(gammaln(a)-gammaln(.5))+gammaln(m/2+2*s)-gammaln(m*a))

def evaluate(row):
    model=geometry(row);m=model['m'];s=float(S);a=.5+s
    factor=gamma_factor(m)
    model['tilted_gamma_shape']=str(F(1,2)+S)
    model['radial_gamma_ratio_arguments']=[str(F(m,2)+2*S),str(F(m,2)+m*S)]
    model['total_gamma_prefactor_diagnostic']=factor
    if m==3:
        assert {t['coefficient'] for t in model['normalized_R_terms']}=={'1/3'}
        ctx.prec=256
        z=ar(S);aa=ar(F(1,2)+S)
        k=3**z*(aa.gamma()/(arb(1)/2).gamma())**3*(arb(3)/2+2*z).gamma()/(arb(3)/2+3*z).gamma()
        model['method']='exact_constant_simplex_polynomial_R_equals_1_over_3'
        model['K_certified_ball']=str(k);model['K_diagnostic']=float(k)
        model['Dirichlet_polynomial_negative_moment']=3**s
    else:
        unit=qmc.Sobol(d=m-1,scramble=True,seed=33033).random_base2(16)
        probs=np.empty((len(unit),m));remaining=np.ones(len(unit))
        for j in range(m-1):
            z=betaincinv(a,(m-j-1)*a,unit[:,j])
            probs[:,j]=remaining*z;remaining*=1-z
        probs[:,-1]=remaining
        assert np.all(probs>0) and np.all(np.isfinite(probs))
        polynomial=np.zeros(len(unit))
        for t in model['normalized_R_terms']:
            polynomial+=float(F(t['coefficient']))*np.prod(probs[:,t['complement']],axis=1)
        assert np.all(polynomial>0) and np.all(np.isfinite(polynomial))
        values=polynomial**(-s)
        model.update({'method':'one_scrambled_Sobol_Dirichlet_not_certified','seed':33033,'samples':len(unit),
                      'Dirichlet_dimension':m-1,'Dirichlet_polynomial_negative_moment':float(np.mean(values)),
                      'K_diagnostic':factor*float(np.mean(values)),
                      'R_range_in_finite_sample':[float(np.min(polynomial)),float(np.max(polynomial))],
                      'smallest_sample_simplex_coordinate':float(np.min(probs)),
                      'simplex_sum_max_error':float(np.max(np.abs(np.sum(probs,axis=1)-1))),
                      'eight_block_K_diagnostic':[factor*float(np.mean(v)) for v in np.array_split(values,8)],
                      'clipping':False})
    return model

def main():
    assert json.loads((P/'registration.json').read_text())['status']=='registered_before_computation'
    manifest=json.loads((SOURCE/'manifest.json').read_text())
    for f in manifest['white_list']:
        assert digest(SOURCE/f['path'])['sha256']==f['sha256'],f['path']
    data=json.loads((SOURCE/'coefficients.json').read_text())
    rows=[]
    if (P/'checkpoint.json').exists():
        cp=json.loads((P/'checkpoint.json').read_text())
        assert cp['script_sha256']==digest(__file__)['sha256'];rows=cp['rows']
    for row in data['models'][len(rows):]:
        if time.monotonic()-START>120:
            print('120 second checkpoint; rerun to resume',flush=True);return
        rows.append(evaluate(row));dump('checkpoint.json',{'script_sha256':digest(__file__)['sha256'],'pid':os.getpid(),'rows':rows})
        print(json.dumps(rows[-1],indent=2),flush=True)
    s=float(S);rho=float(RHO)
    m=math.exp(gammaln(1+2*s)-gammaln(1+s))
    J=4*math.gamma(1.5)*m**(3/rho)*math.exp(sum(float(F(row['theta']))/rho*math.log(row['K_diagnostic']) for row in rows))
    c1=1/(32*math.pi**3)
    I=c1**1.25*math.pi**2/60*J**1.25
    result={'status':'three_fixed_beta_half_joint_moments_diagnostic_complete',
            'parameters':{'beta':'1/2','b':'1/10','rho':str(RHO),'s':str(S)},
            'models':rows,'J_joint_budget_diagnostic':J,'I_joint_budget_diagnostic':I,
            'below_68e_minus6_diagnostic_only':bool(I<6.8e-5),
            'source_geometry_lock_verified':True,'source_files':[digest(SOURCE/x) for x in ('manifest.json','coefficients.json')],
            'script':digest(__file__),'registration':digest(P/'registration.json'),
            'runtime':{'pid':os.getpid(),'python':sys.executable,'seconds':time.monotonic()-START},
            'new_certified_I_upper_bound':False,'new_C5_upper_bound':False,
            'warning':'Only m=3 Gamma expression is certified; m=4,5 Sobol moments and resulting I are diagnostics, never caps.'}
    dump('results.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='models'},indent=2),flush=True)

if __name__=='__main__':main()
