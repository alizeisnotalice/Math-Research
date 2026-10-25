#!/usr/bin/env python3
"""Independent finite checks for the second review's shared theory.
No primary/audit program imports, no Green evaluations or optimization.
"""
from pathlib import Path
from fractions import Fraction
from itertools import product, permutations, combinations
from hashlib import sha256
import argparse,json,time,platform
import sympy as S

def canonical(v):
    v=tuple(v)
    k=next(x for x in v if x)
    return v if k>0 else tuple(-x for x in v)

def check():
    project=Path(__file__).resolve().parents[3]
    snapshot=json.loads((Path(__file__).with_name('source_snapshot.json')).read_text())
    for rel,expected in snapshot['inputs'].items():
        actual=sha256((project/rel).read_bytes()).hexdigest()
        assert actual==expected,(rel,actual,expected)
    cases=[]
    for n,tau in [(3,S.Rational(1,3)),(4,S.Integer(0))]:
        body=[(1,)+x for x in product((-1,1),repeat=n-1)]
        face=[]
        for i,j in combinations(range(n),2):
            for s in [-1,1]:
                v=[0]*n;v[i]=1;v[j]=s;face.append(tuple(v))
        gram=lambda vv:sum((S.Matrix(v)*S.Matrix(v).T for v in vv),S.zeros(n))
        assert gram(body)==2**(n-1)*S.eye(n)
        assert gram(face)==2*(n-1)*S.eye(n)
        weights={canonical(v):S.Integer(1) for v in body}
        if tau:
            weights.update({canonical(v):tau for v in face})
        group_count=0
        for p in permutations(range(n)):
            for signs in product((-1,1),repeat=n):
                transformed={canonical(tuple(signs[i]*v[p[i]] for i in range(n))):w for v,w in weights.items()}
                assert transformed==weights
                group_count+=1
        if n==3:
            basis=[(1,1,1),(1,1,-1),(1,-1,1)]
        else:
            basis=[(1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1)]
        det=abs(S.Matrix.hstack(*(S.Matrix(v) for v in basis)).det())
        assert det=={3:4,4:16}[n]
        pair_count=0
        # Each pair has a basis extension: used by the line-exterior continuity proof.
        directions=list(weights)
        for i,j in combinations(range(len(directions)),2):
            found=False
            for tail in combinations([k for k in range(len(directions)) if k not in (i,j)],n-2):
                d=S.Matrix.hstack(*(S.Matrix(directions[k]) for k in (i,j)+tail)).det()
                if d:
                    found=True;break
            assert found
            pair_count+=1
        # Exact small-cutoff diagnostic g_e=(e1.x)^2 (1-|x|^2/e^2)_+^3.
        z=S.symbols('z',real=True)
        integral=S.integrate((1-z*z)**3,(z,-1,1))
        assert integral==S.Rational(32,35)
        cutoff_sum=S.simplify(sum(w*S.Rational(v[0]**2)/S.sqrt(sum(x*x for x in v)) for v,w in weights.items()))
        assert S.simplify(cutoff_sum-(2**(n-1)/S.sqrt(n)+tau*S.sqrt(2)*(n-1)))==0
        # At x=2*(1,...,1), only that direction reaches support of radial unit-ball bump.
        v0=S.ones(n,1);x0=2*v0
        distances=[]
        for v in directions:
            if canonical(v)==canonical(v0):continue
            w=S.Matrix(v);dsq=S.simplify(x0.dot(x0)-x0.dot(w)**2/w.dot(w))
            assert dsq>1
            distances.append(dsq)
        contact_lb=S.Rational(27,64)/(S.pi*S.sqrt(n)*(2+1/(2*S.sqrt(n)))**2)
        assert contact_lb>0
        # Full radial mass and scalar optimum, fixed dimensions; positive symbols.
        a,c,g,J,I,r=S.symbols('a c g J I r',positive=True)
        alpha=n-1;q=S.Rational(n,n-1);R=(a/c)**S.Rational(1,n-1)
        radial=S.simplify(a*R-c*R**n/n)
        expected=S.Rational(n-1,n)*a**q*c**(-S.Rational(1,n-1))
        assert S.simplify(radial-expected)==0
        B=S.Rational(n-1,n)*J*c**(-S.Rational(1,n-1))/(2**n*(g-c))
        logarithmic=S.simplify(S.diff(B,c)/B)
        assert S.simplify(logarithmic-(n*c-g)/((n-1)*c*(g-c)))==0
        optimum=S.simplify(B.subs(c,g/n))
        assert S.simplify(optimum-n**S.Rational(1,n-1)*J/(2**n*g**q))==0
        face_cost=S.simplify(optimum.subs(J,n*2**n*I))
        assert S.simplify(face_cost**(n-1)-n**n*I**(n-1)/g**n)==0
        cases.append({'n':n,'tau':str(tau),'body_gram':str(gram(body)),'face_gram':str(gram(face)),
          'basis_columns':basis,'absolute_basis_det':int(det),'cube_group_elements_checked':group_count,
          'pair_basis_extensions_checked':pair_count,'small_cutoff_L_at_zero_divided_by_epsilon':str(S.simplify(integral*cutoff_sum/S.pi)),
          'contact_bump_point':[2]*n,'other_line_distance_squared_min':str(min(distances)),
          'strict_positive_contact_lower_bound':str(contact_lb),'mass_radial_formula':str(expected),
          'cost_log_derivative':str(logarithmic),'optimized_cost':str(face_cost),'cost_power':str(S.simplify(face_cost**(n-1)))})
    return {'status':'PASS finite algebra only; continuous proof is audit.md','python':platform.python_version(),
      'sympy':S.__version__,'arithmetic':'SymPy exact rationals and radicals; no floating tolerance',
      'source_hashes_checked':len(snapshot['inputs']),'cases':cases,
      'no_primary_code_imported':True,'no_Green_or_g_I_evaluations':True}
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--save',action='store_true');args=parser.parse_args()
    start=time.monotonic();result=check();result['seconds']=time.monotonic()-start
    text=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.save:Path(__file__).with_name('results.json').write_text(text)
    print(text,end='')
