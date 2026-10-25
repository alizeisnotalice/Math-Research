"""Independent finite frame/cap budget checks; no primary module imports."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
import argparse,hashlib,json,time
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
def ball(q):
    q=F(q);return arb(q.numerator)/q.denominator
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def strictly_inside(x,published):
    lo,hi=map(F,ends(x));a,b=map(F,published);return a<lo<=hi<b

def frames(n,mixed):
    # Construct hypercube vertices by bits, face diagonals as an unoriented set.
    body=[(1,)+tuple(2*((bits>>j)&1)-1 for j in range(n-1))
          for bits in range(2**(n-1))]
    faces=set()
    for i in range(n):
        for j in range(i+1,n):
            for sign in (-1,1):
                v=[0]*n;v[i]=1;v[j]=sign;faces.add(tuple(v))
    return body,body+sorted(faces) if mixed else body

def rebuild_geometry(n,mixed):
    body,vectors=frames(n,mixed)
    determinants={}
    for i,j in combinations(range(len(vectors)),2):
        values=[abs(int(s.Matrix.hstack(s.Matrix(vectors[i]),s.Matrix(vectors[j]),
                         *(s.Matrix(v) for v in rest)).det()))
                for rest in combinations(body,n-2)]
        best=max(values);assert best>0
        determinants[i,j]=best
    dstar=max(abs(int(s.Matrix.hstack(*(s.Matrix(v) for v in basis)).det()))
              for basis in combinations(body,n))
    assert dstar=={3:4,4:16}[n]
    face_counts=[];q_norm2=set();psd_types=set()
    for k in range(n):
        for sign in (-1,1):
            Q={tuple(F(sign*a,v[k]) for a in v) for v in vectors if v[k]}
            face_counts.append(len(Q))
            for q in Q:
                norm2=sum(a*a for a in q);q_norm2.add(norm2)
                assert q[k]==sign and norm2<=n
                tangent=s.Matrix([s.Rational(a.numerator,a.denominator) for j,a in enumerate(q) if j!=k])
                # dist(q+z,span(q))^2 >= |z|^2/|q|^2 by PSD.
                M=(s.Rational(norm2.numerator,norm2.denominator)-1)*s.eye(n-1)-tangent*tangent.T
                eig=M.eigenvals();assert all(a>=0 for a in eig)
                psd_types.add(tuple(sorted((str(a),b) for a,b in eig.items())))
    assert face_counts==[8]*(2*n)
    return body,vectors,determinants,dstar,face_counts,q_norm2,psd_types

def run(root):
    started=time.monotonic();snapshot=json.loads((HERE/'030_source_snapshot.json').read_text())
    for name,h in snapshot.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    author=json.loads((root/'results/round011_face_caps.json').read_text())
    assert len(author['cases'])==5
    # Check all time and cap primitive identities, including the endpoint at zero.
    t,d,c,r=s.symbols('t d c r',positive=True)
    time3=s.integrate(t/(d*d+c*c*t*t),(t,0,1))
    time4=s.integrate(1/(d*d+c*c*t*t),(t,0,1))
    assert s.simplify(time3-s.log(1+c*c/(d*d))/(2*c*c))==0
    assert s.simplify(time4-s.atan(c/d)/(c*d))==0
    primitive=r*r*(s.log(r)**2-s.log(r)+s.Rational(1,2))/2
    assert s.simplify(s.diff(primitive,r)-r*s.log(r)**2)==0
    assert s.limit(primitive,r,0,dir='+')==0
    assert s.diff(s.Rational(3,5)*r**s.Rational(5,3),r)==r**s.Rational(2,3)
    # Use unordered pairs and double the sum; primary uses ordered pair loops.
    geometries={3:rebuild_geometry(3,True),4:rebuild_geometry(4,False)}
    ctx.prec=256;cases=[]
    expected=[F(5892),F(36800,3),F(38400),F(187200),F(5632)]
    for index,(n,tau) in enumerate([(3,F(1,10)),(3,F(1,3)),(3,F(1)),(3,F(3)),(4,F(0))]):
        body,vectors,dets,dstar,face_counts,q_norm2,psd_types=geometries[n]
        m=len(vectors);rates=[F(1) if v in body else tau for v in vectors]
        assert min(rates)>0
        S=2*sum((rates[i]*rates[j]/det for (i,j),det in dets.items()),F(0))
        Cnum=4*m*m*n*S;assert Cnum==expected[index]
        hist=Counter()
        for v in dets.values():hist[v]+=2
        primary=author['cases'][index]
        assert primary['n']==n and F(primary['tau'])==tau
        assert F(primary['pair_sum_exact'])==S
        assert F(primary['C_numerator_over_pi_n'])==Cnum
        assert primary['m']==m and primary['ordered_pairs']==m*(m-1)
        assert {int(k):v for k,v in primary['det_histogram'].items()}==dict(hist)
        assert set(map(tuple,primary['directions']))==set(vectors)
        # Verify each author's explicit completion as well as our independent max.
        for witness in primary['basis_certificate']:
            i,j=witness['i'],witness['j']
            v,w=map(tuple,(primary['directions'][i],primary['directions'][j]))
            rest=list(map(tuple,witness['completion']))
            assert len(rest)==n-2 and all(q in body for q in rest)
            actual=abs(int(s.Matrix.hstack(s.Matrix(v),s.Matrix(w),*(s.Matrix(q) for q in rest)).det()))
            own_i,own_j=sorted((vectors.index(v),vectors.index(w)))
            assert actual==witness['abs_det']==dets[own_i,own_j]
        pi=arb.pi();C=ball(Cnum)/pi**n;cm=ball(min(rates));D=arb(1)/(dstar*(n-1)*pi**n)
        budgets=[]
        for epsilon in (F(1,10**6),F(1,10**8)):
            e=ball(epsilon);assert e<(-arb(1)).exp()
            if n==3:
                A=D+C/(2*cm*cm)*(1+4*m*m*n*n*cm*cm).log()
                B=C/(cm*cm);L=-e.log()
                # Independent positive sqrt route, no exp(log(x)*p).
                cap=arb(2).sqrt()*pi*e*e*(A*A.sqrt()+B*B.sqrt()*(L*L+L+ball(F(1,2))))
            else:
                B=C*pi*m*n/cm
                cap=arb(2).root(3)*4*pi*(D*D.root(3)*e**3/3
                      +ball(F(3,5))*B*B.root(3)*e*e.root(3)**2)
            published=next(q for q in primary['cap_displays'] if F(q['epsilon_exact'])==epsilon)
            assert strictly_inside(cap,published['analytic_budget_interval'])
            total=cap*sum(face_counts)
            budgets.append({'epsilon_exact':str(epsilon),'one_face_one_cap_interval':ends(cap),
                            'one_face_one_cap_ball':cap.str(45),'all_faces_Q_count':sum(face_counts),
                            'all_face_caps_total_interval':ends(total),'all_face_caps_total_ball':total.str(45),
                            'strictly_inside_primary_192bit_interval':True})
        cases.append({'n':n,'tau':str(tau),'m':m,'cmin':str(min(rates)),
                      'pair_sum':str(S),'C_numerator_over_pi_n':str(Cnum),
                      'ordered_pair_det_histogram':dict(hist),'body_max_det':dstar,
                      'face_Q_counts':face_counts,'face_Q_norm_squared':list(map(str,sorted(q_norm2))),
                      'projection_PSD_eigenvalue_types':list(map(str,sorted(psd_types))),
                      'budgets':budgets})
    return {'status':'PASS finite exact pair densities/projection/cap budgets for five locked frames; no global cost certified',
            'source_sha256':snapshot,'precision_bits':ctx.prec,'cases':cases,
            'primitive_identities':{'n3_time':str(time3),'n4_time':str(time4),'log_cap_primitive':str(primitive),'n4_power_primitive':'3*r^(5/3)/5'},
            'implementation':'new SymPy matrix determinants/unordered-pair aggregation, all face tangent PSD matrices, exact rational rates/epsilons, positive sqrt/cuberoot Arb evaluation; no primary imports',
            'shared_libraries':{'SymPy':s.__version__,'Fraction':'stdlib','Arb':'python-flint256bit'},
            'scope':'five frames and ten analytic cap expressions, not Green values/g/fullJ/positive-source certification/new Cn bound; no spectrum cones/cost search',
            'seconds':time.monotonic()-started}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    result=run(args.project);result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:(HERE/'030_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
