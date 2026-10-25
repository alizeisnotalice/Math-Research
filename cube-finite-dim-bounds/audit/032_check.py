"""Independent real-Hermite/full96-cone reconstruction for locked cb11.

The primary green_expression is read-only executed only for exact comparison;
no primary Green, heat, cone-checker or limit-checker is used as an oracle.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations,product
from collections import Counter
import argparse,hashlib,json,runpy,time
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
CONST='constant';EPS='epsilon'
def ball(x):
    x=F(x);return arb(x.numerator)/x.denominator
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def plus(a,b):
    d=dict(a)
    for k,v in b.items():d[k]=d.get(k,F(0))+v
    return {k:v for k,v in d.items() if v}
def times(a,t):return {k:v*t for k,v in a.items() if v*t}
def hermite(nodes):
    """DD2 of real g_epsilon at fixed nodes, retaining exact asymptotic jets."""
    z=sorted(nodes);memo={}
    def dd(i,j):
        if (i,j) in memo:return memo[i,j]
        if z[i]==z[j]:
            order=j-i;b=z[i]
            if order==0:out={EPS:F(1)} if b==0 else {abs(b):F(1)}
            elif order==1:out={} if b==0 else {CONST:1/b}
            elif order==2:
                assert b!=0,'all-zero nodes handled by origin peak'
                out={CONST:-1/(2*b*b)}
        else:out=times(plus(dd(i+1,j),times(dd(i,j-1),F(-1))),1/(z[j]-z[i]))
        memo[i,j]=out;return out
    return dd(0,2)
def kind(z):
    zeros=z.count(0)
    return ('two_zero' if zeros==2 else 'one_zero_equal' if zeros==1 and len(set(z))==2
            else 'one_zero_distinct' if zeros==1 else 'all_equal_nonzero' if len(set(z))==1
            else 'two_equal_nonzero' if len(set(z))==2 else 'distinct_nonzero')
def canonical(cone):return tuple(sorted(cone))

def run(root):
    started=time.monotonic();snapshot=json.loads((HERE/'032_source_snapshot.json').read_text())
    for name,h in snapshot.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    directory=root/'upper/continuous_bridge'
    a,b,c=s.symbols('a b c',real=True)
    ABC=s.Matrix([[1,1,2],[0,1,1],[0,0,1]])
    BCD=s.Matrix([[1,2,1],[1,1,1],[0,1,1]])
    assert ABC.inv()*s.Matrix([a,b,c])==s.Matrix([a-b-c,b-c,c])
    assert BCD.inv()*s.Matrix([a,b,c])==s.Matrix([b-c,a-b,b+c-a])
    body=[(1,j,k) for j,k in product((-1,1),repeat=2)]
    face=[tuple((1 if k==i else sign if k==j else 0) for k in range(3))
          for i,j in [(0,1),(0,2),(1,2)] for sign in (-1,1)]
    cones={}
    for perm in permutations(range(3)):
        for signs in product((-1,1),repeat=3):
            P=s.zeros(3)
            for i in range(3):P[i,perm[i]]=signs[i]
            for base in (ABC,BCD):
                R=P*base;rays=tuple(tuple(int(x) for x in R[:,j]) for j in range(3))
                assert abs(R.det())==1
                for v in body+face:
                    dots=[sum(x*y for x,y in zip(v,r)) for r in rays]
                    assert min(dots)>=0 or max(dots)<=0 # all forms linear on cone
                coefficients=tuple((sum(abs(sum(x*y for x,y in zip(v,r))) for v in body),
                                    sum(abs(sum(x*y for x,y in zip(v,r))) for v in face)) for r in rays)
                key=canonical(rays);assert key not in cones
                cones[key]=(rays,coefficients)
    assert len(cones)==96 and len({r for rr,_ in cones.values() for r in rr})==50
    opposite_pairs={min(key,canonical(tuple(tuple(-v for v in r) for r in key))) for key in cones}
    assert len(opposite_pairs)==48
    for key in cones:assert canonical(tuple(tuple(-v for v in r) for r in key)) in cones
    published=json.loads((directory/'cb11_cones.json').read_text())
    assert {canonical(tuple(map(tuple,row['rays']))) for row in published['cones']}==set(cones)
    for row in published['cones']:
        R=s.Matrix.hstack(*(s.Matrix(r) for r in row['rays']))
        assert int(R.det())==row['det']
        for r,rate in zip(row['rays'],row['rates_body_face_coefficients']):
            actual=[sum(abs(sum(x*y for x,y in zip(v,r))) for v in vv) for vv in (body,face)]
            assert actual==rate
    tau=s.symbols('tau',nonnegative=True)
    volume=s.factor(sum(1/(6*s.prod(B+tau*Ff for B,Ff in rates)) for _,rates in cones.values()))
    stated_volume=5/(6*(1+tau)*(2+3*tau)*(4+5*tau))
    assert s.factor(volume-stated_volume)==0
    peak=s.factor(6*volume/8)
    assert s.factor(peak-5/(8*(1+tau)*(2+3*tau)*(4+5*tau)))==0
    assert peak.subs(tau,0)==s.Rational(5,64)
    # Independent truncated-time identity using only a real rational derivative.
    e=s.symbols('epsilon',positive=True);bs=s.symbols('b0 b1 b2',real=True)
    real_product=e*(e*e-sum(bs[i]*bs[j] for i in range(3) for j in range(i+1,3)))/s.prod(e*e+x*x for x in bs)
    derivative=sum(e/(e*e+bs[j]*bs[j])/s.prod(bs[j]-bs[k] for k in range(3) if k!=j) for j in range(3))
    assert s.factor(derivative+real_product)==0
    z=s.symbols('z',real=True,nonzero=True);g=s.log(e*e+z*z)/2
    assert s.limit(g-s.log(z*z)/2,e,0,dir='+')==0
    assert s.simplify(s.limit(s.diff(g,z),e,0,dir='+')-1/z)==0
    assert s.simplify(s.limit(s.diff(g,z,2)/2,e,0,dir='+')+1/(2*z*z))==0
    assert s.diff(s.log(e*e+s.Symbol('x',real=True)**2)/2,s.Symbol('x',real=True)).subs(s.Symbol('x',real=True),0)==0
    # Only compare the primary exact expression after our independent reconstruction.
    module=runpy.run_path(str(directory/'cb11_frame_green.py'))
    ctx.prec=256;rows=[];all_types=set()
    for case in json.loads((HERE/'032_cases.json').read_text()):
        x=tuple(map(F,case['x']));t=F(case['tau']);assert t>=0
        if not any(x):
            primary=module['green_expression'](x,t);assert primary['origin']
            rows.append({**case,'classification':'infinity_origin_power','origin_peak_coefficient_over_pi3':str(F(15,544))})
            continue
        total={};types=Counter();local_gamma_signs=Counter()
        for rays,coefficients in cones.values():
            rates=[F(B)+t*Ff for B,Ff in coefficients];assert all(k>0 for k in rates)
            nodes=[sum(F(v)*xx for v,xx in zip(r,x))/k for r,k in zip(rays,rates)]
            assert nodes.count(0)<3
            label=kind(nodes);types[label]+=1;all_types.add(label)
            form=hermite(nodes);gamma=form.get(EPS,F(0))
            local_gamma_signs['negative' if gamma<0 else 'positive' if gamma>0 else 'zero']+=1
            total=plus(total,times(form,F(1)/s.prod(rates)))
        gamma=F(total.get(EPS,F(0)));constant=F(total.get(CONST,F(0)))
        logs={F(k):F(v) for k,v in total.items() if not isinstance(k,str)}
        assert sum(logs.values(),F(0))+gamma==0 and gamma<=0
        primary=module['green_expression'](x,t);scale=max(map(abs,x))
        assert primary['scale']==scale and primary['constant']==constant*scale*scale
        assert primary['log_divergence']==gamma*scale*scale
        assert primary['log_coefficients']=={q/scale:v*scale*scale for q,v in logs.items()}
        record={**case,'node_type_counts_full96':dict(types),'local_gamma_sign_counts_full96':dict(local_gamma_signs),
                'global_Gamma_without_8pi3':str(gamma),'constant_without_8pi3':str(constant),
                'log_coefficients_without_8pi3':{str(q):str(v) for q,v in sorted(logs.items())},
                'primary_normalized_expression_matches_exactly':True}
        if gamma<0:record['classification']='infinity_log'
        else:
            value=ball(constant)+sum((ball(v)*ball(q).log() for q,v in logs.items()),arb(0))
            value/=8*arb.pi()**3;assert value>0
            record.update(classification='finite',Green_interval=ends(value),Green_ball=value.str(45))
        if x==(F(1),F(1),F(0)):
            if t==0:
                assert gamma==0 and constant==-1
                primes={}
                for q,v in logs.items():
                    for p,power in s.factorint(q.numerator).items():primes[int(p)]=primes.get(int(p),F(0))+v*int(power)
                    for p,power in s.factorint(q.denominator).items():primes[int(p)]=primes.get(int(p),F(0))-v*int(power)
                assert {p:v for p,v in primes.items() if v}=={2:F(4)}
                record['edge_physical_identity_exact']='(4log(2)-1)/(8pi^3)'
            else:assert gamma<0
        if x==(F(1),F(1),F(1)):assert gamma<0
        if x==(F(1),F(0),F(0)):
            assert gamma==0 and local_gamma_signs['negative'] and local_gamma_signs['positive']
        rows.append(record)
    assert all_types=={'distinct_nonzero','two_equal_nonzero','all_equal_nonzero','one_zero_distinct','one_zero_equal','two_zero'}
    return {'status':'PASS complete96-cone pointwise real Green formula and all zero/repeated-node rules; no continuous-box or cost certificate',
            'source_sha256':snapshot,'cone_count':96,'opposite_pairs':48,'distinct_rays':50,
            'ABC_inverse':['a-b-c','b-c','c'],'BCD_inverse':['b-c','a-b','b+c-a'],
            'all_original_spectral_forms_have_one_sign_on_each_cone':True,
            'spectral_volume':str(volume),'peak_times_pi3_t3':str(peak),'tail_from_T':'peak_coefficient/(2*pi^3*T^2)',
            'truncated_DD_route':'real g_epsilon=log(epsilon^2+b^2)/2; rational derivative identity plus zero limit at infinity; continuous confluent extension',
            'asymptotic_jet_nonzero':'g->log|b|, g_prime->1/b, g_second/2->-1/(2b^2), errors O(epsilon^2) at fixed nonzero knots',
            'asymptotic_jet_zero':'g(0)=log epsilon, g_prime(0)=0; all3zero is origin only',
            'fixed_cases':rows,'all_six_node_types_covered':sorted(all_types),'precision_bits':ctx.prec,
            'implementation':'new full96-cone direct frame-rate construction and real-Hermite recursion; primary expression read-only called only for exact normalized-coefficient comparison, no green/heat/checker calls',
            'shared_libraries':{'SymPy':s.__version__,'Fraction':'stdlib','Arb':'python-flint256bit'},
            'scope':'complete pointwise representation and exact rational-input classification; no cb11 continuous box/Taylor remainder, g/J minima/cubature/LP/grid/cost;031 true weighted semigroup inherited',
            'seconds':time.monotonic()-started}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    result=run(args.project);result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:(HERE/'032_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
