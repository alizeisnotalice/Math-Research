"""One preregistered exact symbolic batch; no cost/parameter optimization."""
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
import json,hashlib
OUT=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds')

def sphere_values(n,B4,B6,R=Q(1)):
    alpha=Q(3,n+2);beta=Q(15,n+8);gamma=Q(30,(n+8)*(n+4))
    d=B6*R**6;e=B4*R**4-beta*d;const=-alpha*B4*R**4+gamma*d
    vals=[];degenerate=[]
    for support in range(1,n+1):
        vals.append(d/Q(support**2)+e/Q(support)+const)
        if d==0:continue
        T=-2*e/(3*d)
        for p in range(1,support):
            q=support-p
            if p==q:
                if p*T==1:degenerate.append({'support':support,'p':p,'q':q,'sum_u_v':str(T),'representative':'uniform support'})
                continue
            u=(1-q*T)/(p-q);v=T-u
            if u>0 and v>0:
                vals.append(d*(p*u**3+q*v**3)+e*(p*u**2+q*v**2)+const)
    return min(vals),max(vals),degenerate

def run():
    sphere=[]
    for n,B4,B6 in ((3,Q(-3,22),Q(1)),(4,Q(1,2),Q(1))):
        lo,hi,degen=sphere_values(n,B4,B6)
        assert degen
        p=1 if n==3 else 2
        T=Q(1,p)
        u=T/4;v=3*T/4
        d=B6;e=B4-Q(15,n+8)*B6;const=-Q(3,n+2)*B4+Q(30,(n+8)*(n+4))*B6
        nonuniform=d*p*(u**3+v**3)+e*p*(u**2+v**2)+const
        uniform=d/Q((2*p)**2)+e/Q(2*p)+const
        assert nonuniform==uniform
        sphere.append({'n':n,'B4':str(B4),'B6':str(B6),'exact_global_range':[str(lo),str(hi)],'degenerate_stationary_families':degen,'nonuniform_equals_uniform':str(uniform)})
    # Equal multiplicities in cube remain finite, unlike the sphere degeneracy.
    # n4, all four y free: two y=1/4 and two y=3/4, B4=B6=2/3.
    n=4;B4=s.Rational(2,3);B6=s.Rational(2,3)
    ys=s.symbols('y0:4',positive=True);S=sum(ys);P2=sum(y*y for y in ys);P3=sum(y**3 for y in ys)
    f=1/S-B4*(P2-s.Rational(1,2)*S*S)-B6*(P3-s.Rational(5,4)*S*P2+s.Rational(5,16)*S**3)
    point=dict(zip(ys,[s.Rational(1,4)]*2+[s.Rational(3,4)]*2))
    grad=[s.factor(s.diff(f,y).subs(point)) for y in ys]
    assert grad==[0]*4
    t=s.symbols('t',real=True)
    # On equal-multiplicity linear reduction u+v=1 and S=2,
    # the remaining stationarity is 1/4 - (2/3)u(1-u) = 0? Derive it directly.
    grad_on_line=s.factor(s.diff(f,ys[0]).subs(dict(zip(ys,[t,t,1-t,1-t]))))
    assert grad_on_line!=0
    roots=s.solve(grad_on_line,t)
    assert roots==[s.Rational(1,4),s.Rational(3,4)]
    cube={'n':4,'B4':'2/3','B6':'2/3','equal_multiplicity':'2+2','stationary_y':['1/4','1/4','3/4','3/4'],'gradient':list(map(str,grad)),'reduced_polynomial':str(grad_on_line),'admissible_roots':list(map(str,roots)),'scope':'stationary candidate completeness only; not a global cube minimum or kernel certificate'}
    radial=[]
    t=s.symbols('t',real=True);A,D,E=s.symbols('A D E')
    for n,expr in ((3,1-A*t-D*t**5-E*t**7),(4,1-A*t-D*t**3-E*t**4)):
        sol=s.solve([expr.subs(t,v) for v in (1,2,3)],(A,D,E),dict=True)[0]
        poly=s.factor(expr.subs(sol));factors=s.factor_list(poly)
        intervals=s.Poly(poly,t).intervals(eps=s.Rational(1,10**12))
        positives=[(a,b,m) for (a,b),m in intervals if a>0]
        assert len(positives)==3 and all(m==1 for a,b,m in positives)
        signs=[s.sign(poly.subs(t,v)) for v in (s.Rational(1,2),s.Rational(3,2),s.Rational(5,2),s.Rational(7,2))]
        assert signs==[1,-1,1,-1]
        radial.append({'n':n,'variable':'r' if n==3 else 'z=r^2','A':str(sol[A]),'D=B4*h4':str(sol[D]),'E=B6*h6':str(sol[E]),'factorization':str(poly),'positive_roots_with_multiplicity':[{'interval':[str(a),str(b)],'multiplicity':m} for a,b,m in positives],'interval_signs_at_half_points':list(map(str,signs)),'scope':'univariate radial mechanism example only; no global outer-ball negativity or cube certificate claimed'})
    rec={'status':'PASS exact finite symbolic batch','sympy':s.__version__,'sphere':sphere,'cube_equal_multiplicity':cube,'radial_three_root_examples':radial,'new_upper_bound':False,'shared_arithmetic':'derivations and implementation share normalized formulas/SymPy; independent audit required'}
    (OUT/'cb7_results.json').write_text(json.dumps(rec,indent=2)+'\n')
    paths=[MAIN/'rounds/007-registration.md',MAIN/'STATE.md',MAIN/'results/round007_h6_initial.json',MAIN/'src/check_h6_identities.py',MAIN/'upper/continuous_bridge/cb6_legality.md']
    (OUT/'cb7_sources.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2)+'\n')
    print(json.dumps(rec,indent=2))
if __name__=='__main__':run()
