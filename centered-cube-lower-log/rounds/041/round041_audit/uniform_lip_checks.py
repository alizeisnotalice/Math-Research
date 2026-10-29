"""Independent exhaustive piecewise-linear uniform MGF/stoploss controls."""
from pathlib import Path
from fractions import Fraction as F
import itertools,math,json,hashlib
HERE=Path(__file__).resolve().parent
stop=mgf=0;maxratio=0.;affine_equal=0
for slopes in itertools.product([-1,0,1],repeat=4):
    values=[F(0)]
    for s in slopes:values.append(values[-1]+F(s,4))
    mean=sum((F(1,4)*(a+b)/2 for a,b in zip(values,values[1:])),F(0))
    vals=[x-mean for x in values]
    assert min(vals)>=F(-1,2)and max(vals)<=F(1,2)
    for t in [F(-1,2),F(-1,4),F(0),F(1,4),F(1,2)]:
        v=F(0)
        for a,b in zip(vals,vals[1:]):
            a-=t;b-=t
            if min(a,b)>=0:v+=F(1,4)*(a+b)/2
            elif max(a,b)>0:v+=F(1,4)*max(a,b)**2/(2*abs(a-b))
        upper=(F(1,2)-t)**2/2
        assert v<=upper;stop+=1
    for t in [-8.,-2.,-1.,.5,1.,2.,8.]:
        total=0.
        for a,b in zip(vals,vals[1:]):
            aa,bb=float(a),float(b);d=t*(bb-aa)
            total+=.25*math.exp(t*aa)*(math.expm1(d)/d if d else 1.)
        exact_uniform=math.sinh(t/2)/(t/2);gauss=math.exp(t*t/24)
        assert total<=exact_uniform*(1+1e-13)and exact_uniform<=gauss*(1+1e-13)
        maxratio=max(maxratio,total/exact_uniform);mgf+=1
        if all(s==1 for s in slopes)or all(s==-1 for s in slopes):
            assert abs(total/exact_uniform-1)<1e-13;affine_equal+=1
out={'exhaustive_functions':81,'slopes':[-1,0,1],'equal_piece_width':'1/4',
     'exact_Fraction_stoploss_checks':stop,'analytic_piecewise_integral_MGF_checks':mgf,
     'nonzero_t_affine_equalities':affine_equal,'max_MGF_over_uniform_MGF':maxratio,
     'reviewed_proof_SHA256':hashlib.sha256((HERE.parent/'round041_mcd/uniform_PROOF.md').read_bytes()).hexdigest(),
     'scope':'No random samples; exact stoploss controls and ordinary-FP exponential integrals; not a source FP certificate'}
p=HERE/'uniform_lip_checks.json'
if p.exists():raise RuntimeError('Refusing replacement')
p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
