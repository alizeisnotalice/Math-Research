"""Finite independent checks; only Python Fraction and SymPy, no author imports."""
from fractions import Fraction as F
from pathlib import Path
import json, hashlib
import sympy as s
HERE=Path(__file__).resolve().parent
SOURCES={
 'cb7':'/Users/zhengzhihao/Desktop/T/Math-Research-finite-bridge/cube-finite-dim-bounds/upper/continuous_bridge/cb7_stratification.md',
 'vg7':'/Users/zhengzhihao/Desktop/T/Math-Research-finite-variation/cube-finite-dim-bounds/lower/variation_geometry/vg7_joint_witness.md'}
EXPECTED={'cb7':'e55314369108302cc795f117233c32161361711bf6aab14d222402b166542caf','vg7':'034c633ef1ab779ad669d3eba1d6239eff4d924def91d85fddaa85948c4ce220'}
def sphere(n,R,B4,B6):
    alpha=F(3,n+2);beta=F(15,n+8);gamma=F(30,(n+8)*(n+4))
    d=B6*R**6;e=B4*R**4-beta*d;C=-alpha*B4*R**4+gamma*d
    out=[d/F(k*k)+e/F(k)+C for k in range(1,n+1)]
    if d:
        T=-2*e/(3*d)
        for p in range(1,n):
            for q in range(1,n-p+1):
                if p==q: continue
                a=(1-q*T)/(p-q);b=T-a
                if a>0 and b>0:out.append(d*(p*a**3+q*b**3)+e*(p*a*a+q*b*b)+C)
    return min(out),max(out)
# Direct endpoint partition, rather than author normal-interval formulas.
def interval_length(intervals):
    cuts=sorted(set(z for lo,hi in intervals if lo<hi for z in (lo,hi)))
    return sum((b-a for a,b in zip(cuts,cuts[1:]) if any(lo<(a+b)/2<hi for lo,hi in intervals)),F(0))
def run():
    hashes={k:hashlib.sha256(Path(v).read_bytes()).hexdigest() for k,v in SOURCES.items()}; assert hashes==EXPECTED
    guards=[]
    for n,R,A,B4,B6,expected in [
      (3,F(16,5),F('.4089479272201547'),F('.005004477192161995'),F('.0009386738283819751'),F(2152431878903217985131,1503906250000000000000000)),
      (4,F(18,5),F('.18157763153667977'),F('.004'),F('.0005'),F(8178647106137,8100000000000000000))]:
        lo,hi=sphere(n,R,B4,B6);guard=A+lo-R**(2-n);assert guard==expected and guard>0
        guards.append({'n':n,'margin':str(guard)})
    assert sphere(3,F(1),F(-3,22),F(1))==(F(-1,35),F(26,315))
    assert sphere(4,F(1),F(1,2),F(1))==(F(-11,144),F(5,16))
    symbolic=[]
    for n in (3,4):
        ys=s.symbols('y:'+str(n), positive=True);B4,B6=s.symbols('B4 B6');S=sum(ys);P2=sum(t*t for t in ys);P3=sum(t**3 for t in ys)
        alpha=s.Rational(3,n+2);beta=s.Rational(15,n+8);gamma=s.Rational(30,(n+8)*(n+4));nu=s.Rational(2-n,2)
        fun=S**nu-B4*(P2-alpha*S*S)-B6*(P3-beta*S*P2+gamma*S**3)
        diff=s.diff(fun,ys[0])-s.diff(fun,ys[1]);claimed=(ys[0]-ys[1])*(-2*B4-3*B6*(ys[0]+ys[1])+2*beta*B6*S)
        assert s.simplify(diff-claimed)==0
        for ell in range(1,n//2+1):assert 1-s.Rational(10*ell,n+8)!=0
        r,R=s.symbols('r R',positive=True)
        for degree in (0,4,6):assert s.integrate(r**(n-1+degree),(r,0,R))==R**(n+degree)/(n+degree)
        symbolic.append(n)
    # Full local normal length including both old fronts; joint only adds its true interval.
    left,right=F(-3,32),F(3,32);born=[None,None,None,(F(-1,32),F(1,32))];Ds=[]
    for z in born:
        intervals=[(F(-2),left),(right,F(2))]+([] if z is None else [z])
        Ds.append(interval_length(intervals)-4)
    assert Ds==[F(-3,16),F(-3,16),F(-3,16),F(-1,8)]
    mixed=Ds[3]-Ds[1]-Ds[2]+Ds[0];assert mixed==F(1,16)
    P=F(1,8);V0=F(1,2);jac=F(3,16);cost=[F(0),F(1,2),F(1,2),F(1)]
    G=[P*d+V0*(jac-w) for d,w in zip(Ds,cost)];assert G==[F(9,128),F(-23,128),F(-23,128),F(-27,64)]
    assert G[3]-G[1]-G[2]+G[0]==F(1,128)
    assert F(1,8)-F(1,32)>F(1,32) # moving tied maximum makes born empty.
    q=F(2049,2048);sigma=q**4-1;h=1+F(3,16)*sigma;l=1/(2*h)-F(1,2);b=q/(2*h)-F(1,2)
    assert 0<b<-l and h**4>1+sigma/2 and sigma<F(1,256)
    length=interval_length([(F(-2),l),(-l,F(2)),(-b,b)]);V=P*length
    assert V==F(563121827635203,1126312525791244)
    assert V-P*interval_length([(F(-2),l),(-l,F(2))])==F(34284216317,2252625051582488)
    assert h*V/(1+sigma)==F(20856363986489,41781501527232)
    return {'status':'PASS','source_sha256':hashes,'independent_guard_margins':guards,'symbolic_dimensions':symbolic,'normal_derivatives':list(map(str,Ds)),'mixed_normal':str(mixed),'complete_quotient_derivatives':list(map(str,G)),'mixed_quotient':'1/128','finite_sigma':str(sigma),'finite_volume':str(V),'scope':'bounded exact identities and endpoint partition, no author imports, no search; analytical completeness separately audited','libraries':{'python':'int/Fraction','sympy':s.__version__}}
if __name__=='__main__':
    out=run();out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();HERE.joinpath('016_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
