"""Pure rational certificate against the posterior-only stopping relaxation, not G_c."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from hashlib import sha256
import json

HERE=Path(__file__).resolve().parent
REG=HERE/'posterior_only_weak_obstruction_registration_20261007.json'
OUT=HERE/'posterior_only_weak_obstruction_results_20261007.json'
assert REG.exists() and not OUT.exists(), "Preserve registration and previous outputs"

def qr(x): return str(x)
def receipt(x):
    def digest(k):
        return sha256(k.to_bytes(max(1,(k.bit_length()+7)//8),'big')).hexdigest()
    return {'rational':str(x),'numerator_binary_sha256':digest(x.numerator),
            'denominator_binary_sha256':digest(x.denominator)}
def basis(n,k,r):
    return r**k*(1-r)**(n-k)
def derivative(n,k,r):
    if k==0: return -n*(1-r)**(n-1)
    if k==n: return n*r**(n-1)
    return k*r**(k-1)*(1-r)**(n-k)-(n-k)*r**k*(1-r)**(n-k-1)

rounds=[]; total=0
for n in (8,32,128):
    tests=[]
    def check(name,condition):
        assert condition,(n,name)
        tests.append(name)
    lam=Q(3,2)
    p=[Q(k,n) for k in range(n+1)]
    b=[basis(n,k,p[k]) for k in range(n+1)]
    bins=[comb(n,k) for k in range(n+1)]
    v=[Q(bins[k])*b[k]/lam for k in range(n+1)]
    amp=[lam/b[k] for k in range(n+1)]
    D=sum(Q(bins[k])*b[k] for k in range(n+1))
    H=lam*max([Q(1)]+[1/(Q(bins[k])*b[k]) for k in range(1,n+1)])
    ve=v.copy(); ve[0]=1/H
    ae=amp.copy(); ae[0]=H
    rv=sorted(set([Q(0),Q(1),Q(1,n+1),Q(1,n),Q(1,4),Q(1,2),Q(3,4),1-Q(1,n)]))
    for k in range(n+1):
        check(f'layer_mass_{k}',v[k]*amp[k]==bins[k])
        check(f'enhanced_layer_mass_{k}',ve[k]*ae[k]==bins[k])
        check(f'positive_finite_box_{k}',0<v[k]<=Q(2,3) and 0<ve[k]<=Q(2,3))
        # E_k=[2k,2k+v_k] x [0,1]^(n-1): positive gaps guarantee disjointness.
        if k<n: check(f'box_gap_{k}',2*k+v[k]<2*(k+1) and 2*k+ve[k]<2*(k+1))
        check(f'peak_height_{k}',amp[k]*b[k]==lam)
        check(f'enhanced_peak_{k}',ae[k]*b[k]==(H if k==0 else lam))
        check(f'height_cap_{k}',ae[k]<=bins[k]*H)
        posterior=[Q(0)]*(n+1); posterior[k]=amp[k]*b[k]/lam
        mean=sum(j*posterior[j] for j in range(n+1))
        variance=sum((j-mean)**2*posterior[j] for j in range(n+1))
        check(f'posterior_normalization_{k}',sum(posterior)==1)
        check(f'posterior_mean_{k}',mean==n*p[k])
        check(f'posterior_variance_{k}',variance==0 and variance<=n*p[k]*(1-p[k]))
        if 0<k<n:
            check(f'stationary_{k}',derivative(n,k,p[k])==0)
            for side in (-1,1):
                r=p[k]+Q(side,2*n)
                check(f'derivative_sign_{k}_{side}',derivative(n,k,r)*side<0)
                check(f'nearby_peak_{k}_{side}',basis(n,k,r)<b[k])
            for r in (Q(1,4),Q(1,2),Q(3,4)):
                identity=r**(k-1)*(1-r)**(n-k-1)*(k-n*r)
                check(f'derivative_identity_{k}_{r}',derivative(n,k,r)==identity)
        else:
            check(f'endpoint_global_peak_{k}',b[k]==1)
    for r in rv:
        mass=sum(v[k]*amp[k]*basis(n,k,r) for k in range(n+1))
        masse=sum(ve[k]*ae[k]*basis(n,k,r) for k in range(n+1))
        check(f'total_mass_{r}',mass==1)
        check(f'enhanced_total_mass_{r}',masse==1)
        check(f'binomial_identity_{r}',sum(Q(bins[k])*basis(n,k,r) for k in range(n+1))==1)
    check('initial_mass_one',ve[0]*ae[0]==1)
    check('H_strictly_exceeds_lambda',H>lam)
    check('enhanced_weak_sup_above_one',D-1+lam/H>1)
    sup_enh=D-1+lam/H
    thresholds=[]; last=Q(0); laste=Q(0)
    for m in (1,2,4,8,16,32,64):
        t=lam*(1-Q(1,2**m))
        vol=sum(v[k] for k in range(n+1) if amp[k]*b[k]>t)
        vole=sum(ve[k] for k in range(n+1) if ae[k]*b[k]>t)
        ratio=t*vol; ratioe=t*vole
        check(f'strict_level_all_boxes_{m}',vol==sum(v) and vole==sum(ve))
        check(f'weak_ratio_exact_{m}',ratio==(1-Q(1,2**m))*D)
        check(f'enhanced_weak_ratio_exact_{m}',ratioe==(1-Q(1,2**m))*sup_enh)
        check(f'monotone_strict_approach_{m}',last<ratio<D and laste<ratioe<sup_enh)
        check(f'exact_gap_{m}',D-ratio==D/2**m and sup_enh-ratioe==sup_enh/2**m)
        last,laste=ratio,ratioe
        thresholds.append({'m':m,'threshold':qr(t),'baseline_ratio':qr(ratio),
                           'enhanced_ratio':qr(ratioe)})
    vol_lam=sum(v[k] for k in range(n+1) if amp[k]*b[k]>lam)
    vole_lam=sum(ve[k] for k in range(n+1) if ae[k]*b[k]>lam)
    vole_H=sum(ve[k] for k in range(n+1) if ae[k]*b[k]>H)
    check('strict_at_lambda_baseline_empty',vol_lam==0)
    check('strict_at_lambda_enhanced_E0_only',vole_lam==ve[0])
    check('strict_at_H_enhanced_empty',vole_H==0)
    check('above_lambda_branch_sup_one',sup_enh>1 and H*ve[0]==1)
    rows=[{'k':k,'p':qr(p[k]),'b':qr(b[k]),'volume':qr(v[k]),
           'amplitude':qr(amp[k]),'enhanced_volume':qr(ve[k]),
           'enhanced_amplitude':qr(ae[k])} for k in range(n+1)]
    rounds.append({'n':n,'lambda':qr(lam),'H':receipt(H),
                   'baseline_exact_weak_sup':receipt(D),
                   'enhanced_exact_weak_sup':receipt(sup_enh),
                   'r_values':[qr(r) for r in rv],
                   'strict_thresholds':thresholds,'layer_rows':rows,
                   'passed_checks':len(tests),'test_names':tests})
    total+=len(tests)
res={'status':'PASS_EXACT','checks':total,'rounds':rounds,
     'registration_sha256':sha256(REG.read_bytes()).hexdigest(),
     'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
     'scope':'Posterior-only weak stopping relaxation obstruction, not an original G_c counterexample.'}
OUT.write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':res['status'],'checks':total,
                  'round_checks':[(r['n'],r['passed_checks']) for r in rounds],
                  'weak_sup_displays_only':[(r['n'],float(Q(r['baseline_exact_weak_sup']['rational'])),
                                            float(Q(r['enhanced_exact_weak_sup']['rational']))) for r in rounds]}))
