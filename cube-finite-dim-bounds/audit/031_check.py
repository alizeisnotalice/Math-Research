"""Independent Gram/symmetry/positive-power target checks for031.

No Green formula, sampled minimum, cubature, LP, or primary module imports.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse,hashlib,json,time
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
def ball(x):
    x=F(x);return arb(x.numerator)/x.denominator
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def normalize(v):
    first=next(t for t in v if t)
    return tuple(v) if first>0 else tuple(-t for t in v)

def run(root):
    started=time.monotonic();snapshot=json.loads((HERE/'031_source_snapshot.json').read_text())
    for name,h in snapshot.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    target=json.loads((root/'results/round011_candidate_targets.json').read_text())
    primary=json.loads((root/'results/round011_weighted_frame_checks.json').read_text())
    assert primary['proof_sha256']==snapshot['proofs/round011_weighted_frame_transfer.md']
    tau=s.symbols('tau',positive=True);r=s.symbols('r',positive=True)
    epsilon=s.symbols('epsilon',positive=True);g,I,J=s.symbols('g I J',positive=True)
    cases=[]
    for n in (3,4):
        body={(1,)+tuple(2*((k>>j)&1)-1 for j in range(n-1)) for k in range(2**(n-1))}
        face=set()
        for i,j in combinations(range(n),2):
            for sign in (-1,1):
                v=[0]*n;v[i]=1;v[j]=sign;face.add(tuple(v))
        mb,mf=len(body),len(face)
        gb=sum((s.Matrix(v)*s.Matrix(v).T for v in body),s.zeros(n))
        gf=sum((s.Matrix(v)*s.Matrix(v).T for v in face),s.zeros(n))
        assert gb==mb*s.eye(n) and gf==2*(n-1)*s.eye(n)
        assert mb==2**(n-1) and mf==n*(n-1)
        for vectors in (body,face):
            # Generators of the full cube group, applied as unoriented lines.
            for axis in range(n):
                transformed={normalize(tuple(-a if j==axis else a for j,a in enumerate(v))) for v in vectors}
                assert transformed==vectors
            for axis in range(n-1):
                transformed=set()
                for v in vectors:
                    w=list(v);w[axis],w[axis+1]=w[axis+1],w[axis]
                    transformed.add(normalize(w))
                assert transformed==vectors
        xi=s.symbols('x0:'+str(n),real=True)
        symbol=sum(s.Abs(sum(a*b for a,b in zip(v,xi))) for v in body)+tau*sum(s.Abs(sum(a*b for a,b in zip(v,xi))) for v in face)
        scaled=symbol.subs({x:r*x for x in xi},simultaneous=True)
        assert s.simplify(scaled-r*symbol)==0
        basis=(s.Matrix([[1,1,1],[-1,1,-1],[-1,-1,1]]) if n==3
               else s.Matrix([[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]]))
        assert all(tuple(basis[:,j]) in body for j in range(n))
        det=abs(int(basis.det()));assert det=={3:4,4:16}[n]
        lower=s.sqrt(mb);upper=mb+tau*s.sqrt(mf*2*(n-1))
        assert s.simplify(upper-(mb+tau*(n-1)*s.sqrt(2*n)))==0
        weighted_lengths=sum(s.sqrt(sum(a*a for a in v)) for v in body)+tau*sum(s.sqrt(sum(a*a for a in v)) for v in face)
        jump_rate=2*(mb+tau*mf)/(s.pi*epsilon)
        q=s.Rational(n,n-1);c=g/n
        rawmass=s.Rational(n-1,n)*J*c**(-s.Rational(1,n-1))
        cost=s.simplify(rawmass/(2**n*(g-c)))
        stated=n**s.Rational(1,n-1)*J/(2**n*g**q)
        assert s.simplify(cost-stated)==0
        byI=s.simplify(cost.subs(J,n*2**n*I))
        powered=s.simplify(byI**(n-1))
        assert s.simplify(powered-n**n*I**(n-1)/g**n)==0
        U={3:8,4:16}[n]
        coeff=F(U**(n-1),n**n)
        row=next(v for v in primary['cases'] if v['n']==n)
        assert F(row['coefficient_exact'])==coeff
        cases.append({'n':n,'body_count':mb,'face_count':mf,'body_gram':str(gb),'face_gram':str(gf),
                      'explicit_body_basis_det':det,'symbol_lower_norm_factor':str(lower),
                      'symbol_upper_norm_factor':str(upper),'small_cutoff_weighted_length_sum':str(weighted_lengths),
                      'all_direction_truncation_rate':str(jump_rate),'cube_group_preserves_weight_classes':True,
                      'positive_homogeneity_exact':True,'c_choice':str(c),'H_mass':str(rawmass),
                      'cost_from_J':str(cost),'cost_from_I':str(byI),'cost_power_n_minus_1':str(powered),
                      'breakthrough_I_power_less_than_coefficient_times_g_power':str(coeff),
                      'current_U':U})
    gl=F(target['g_lower_target']);iu=F(target['I_upper_target']);assert gl>0 and iu>0 and F(target['tau'])==F(1,3)
    slack=64*gl**3-27*iu**2
    assert slack>0 and slack==F(56653,250000000000)==F(target['exact_slack_64g3_minus27I2'])
    squared=27*iu*iu/(gl**3);assert squared<64
    ctx.prec=256;conditional=ball(squared).sqrt()
    lo,hi=map(F,ends(conditional));a,b=map(F,target['conditional_cost_interval'])
    assert a<lo<=hi<b and conditional<8
    return {'status':'PASS exact weighted Gram/symmetry/scales/cost targets only; g and I target inequalities are not certified',
            'source_sha256':snapshot,'cases':cases,
            'target':{'tau':'1/3','g_lower_target':str(gl),'I_upper_target':str(iu),
                      'exact_slack_64g3_minus27I2':str(slack),'conditional_cost_squared':str(squared),
                      'conditional_cost_interval':ends(conditional),'conditional_cost_ball':conditional.str(55),
                      'strictly_inside_primary_192bit_interval':True,'conditional_cost_below_8':True,
                      'g_and_I_are_targets_only':True},
            'proof_steps_not_numerically_substituted':'weighted cutoff, obstacle/infinite-source closure, global positive-part source, one-D/all-common-r/general-L1 migration are reviewed in031_audit.md using026 dependencies',
            'precision_bits':ctx.prec,'shared_libraries':{'SymPy':s.__version__,'Fraction':'stdlib','Arb':'python-flint256bit'},
            'scope':'n3/n4 finite body+face with fixed positive tau (tau0 deletes face), true semigroup migration and conditional target algebra, no cb11/Green values/g/I integration or new bound',
            'seconds':time.monotonic()-started}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    result=run(args.project);result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:(HERE/'031_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
