"""One finite exact symbolic batch for the cb6 continuous legality interface.
No cost optimization, numerical mesh, old Q1 test, or LP solve.
"""
from pathlib import Path
import sympy as sp
import hashlib,json
OUT=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds')
AUDIT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-audit/cube-finite-dim-bounds')

def run():
    records=[]
    r=sp.symbols('r',positive=True)
    for n in (3,4):
        xs=sp.symbols('x0:'+str(n),real=True)
        s=sum(x*x for x in xs)
        H=sum(x**4 for x in xs)-sp.Rational(3,n+2)*s*s
        lap=sp.expand(sum(sp.diff(H,x,2) for x in xs))
        assert lap==0
        low=-sp.Rational(2*(n-1),n*(n+2));high=sp.Rational(n-1,n+2)
        axis=sp.simplify(H.subs(dict(zip(xs,[1]+[0]*(n-1)))))
        # Homogeneous sphere diagonal: x_i=1/sqrt(n).
        diagonal=sp.simplify(H.subs({x:1/sp.sqrt(n) for x in xs}))
        assert axis==high and diagonal==low
        sphere_max=sp.Rational(2)**(2-n)-1-sp.Rational(1,100)*low*2**4
        assert sphere_max<0
        expected={3:-sp.Rational(343,750),4:-sp.Rational(71,100)}[n]
        assert sphere_max==expected
        # Normalized W=1 radial small example, A=1,B=0,R=2.
        sigma=2*sp.pi**(sp.Rational(n,2))/sp.gamma(sp.Rational(n,2))
        cn=sp.simplify(1/((n-2)*sigma))
        mass=sp.simplify(cn*sigma*sp.integrate(r-r**(n-1),(r,0,1)))
        assert mass==sp.Rational(1,2*n)
        source_density=cn*(n-2)
        assert sp.simplify(source_density*sigma)==1
        # Core mass at rho=1/2; harmonic quartic integrates to zero.
        rho=sp.Rational(1,2)
        core_mass=sp.simplify(cn*sigma*(rho**2/2-rho**n/n))
        assert core_mass>0
        # Direct fundamental radial Laplacian off origin.
        phi=cn*r**(2-n)
        assert sp.simplify(sp.diff(phi,r,2)+(n-1)*sp.diff(phi,r)/r)==0
        records.append({'n':n,'Delta_H4':str(lap),'unit_sphere_H4_range':[str(low),str(high)],
            'c_n':str(cn),'nonradial_legality_example':{'normalized_c':'1','a':'1','b':'1/100','R':'2','sphere_max_F':str(sphere_max),'cube_certificate':False},
            'W1_radial_example':{'A':'1','B':'0','R':'2','free_radius':'1','mass':str(mass),'source_mass':'1','core_mass_rho_half':str(core_mass),'cube_certificate':False}})
    x,y,z,e=sp.symbols('x y z e',real=True,positive=None)
    P=x*x-y*y
    beta=(P+sp.sqrt(P*P+e*e))/2
    actual=sp.simplify(sum(sp.diff(beta,t,2) for t in (x,y,z)))
    expected=2*e*e*(x*x+y*y)/(P*P+e*e)**sp.Rational(3,2)
    assert sp.simplify(actual-expected)==0
    # Pole, outer radius, and cube interfaces are separate from this local model.
    rec={'status':'PASS exact symbolic finite batch','sympy':sp.__version__,'dimensions':records,
         'critical_zero_local_model':{'P':'x^2-y^2','critical_set':'x=y=0','Delta_smoothed_positive_part':str(expected),'sign':'nonnegative for real x,y,z and epsilon>0'},
         'scope':'legality and normalization cross-check only; no candidate cube domination or new U3/U4','shared_arithmetic':'formulas shared with analytic proof; not independent mathematical audit'}
    (OUT/'cb6_results.json').write_text(json.dumps(rec,indent=2)+'\n')
    # Locate receipt names rather than silently guessing audit filenames.
    audit_paths=sorted(AUDIT.glob('audit/008*'))
    upper=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-upper/cube-finite-dim-bounds/upper')
    paths=[MAIN/'STATE.md',MAIN/'results/bounds_ledger.json']+[p for p in audit_paths if p.is_file()]
    paths += [upper/name for name in ('u05_results.json','u05_h4_scan.py') if (upper/name).is_file()]
    (OUT/'cb6_sources.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2)+'\n')
    print(json.dumps(rec,indent=2))
if __name__=='__main__':run()
