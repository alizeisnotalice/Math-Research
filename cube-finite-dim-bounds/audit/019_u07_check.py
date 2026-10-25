"""New locked U07 parameters only; no author imports or optimization."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,argparse
import sympy as s
HERE=Path(__file__).resolve().parent
HASHES={'u07_certificate.md':'edbc2a59331019783dbd6fcfe5fced106a119d03081e4c667ca71a538aff3eb2','u07_parameter_check.py':'bcc1847b9f3a328d8b5b44ddaed7511ca44920ee725b31eec069c3308ed39051'}
def angular_min(n,R,b4,b6):
    d=b6*R**6;e=b4*R**4-F(15,n+8)*d;C=-F(3,n+2)*b4*R**4+F(30,(n+8)*(n+4))*d
    vals=[d/F(k*k)+e/F(k)+C for k in range(1,n+1)]
    T=-2*e/(3*d)
    for p in range(1,n):
      for q in range(1,n-p+1):
        if p!=q:
          u=(1-q*T)/(p-q);v=T-u
          if u>0 and v>0:vals.append(d*(p*u**3+q*v**3)+e*(p*u*u+q*v*v)+C)
    return min(vals),len(vals)
def run(root):
    for n,h in HASHES.items():assert hashlib.sha256((root/n).read_bytes()).hexdigest()==h
    rows=[]
    for n,A,b4,b6,R,core in [(3,F('.41061707799980696'),F('.0052535036230922875'),F('.0009893885508363638'),F('3.1643099997673794'),F(21,10)),(4,F('.14344247730484339'),F('.002778738033100787'),F('.0003790027706453'),F('3.3162574578337103'),F(23,10))]:
      yy=s.symbols('y:'+str(n),positive=True);p=s.Rational(n-2,2);a=s.Rational(3,n+2);beta=s.Rational(15,n+8);gamma=s.Rational(30,(n+8)*(n+4));S=sum(yy);P2=sum(y*y for y in yy);P3=sum(y**3 for y in yy);B4,B6=s.symbols('B4 B6');f=S**(-p)-A-B4*(P2-a*S*S)-B6*(P3-beta*S*P2+gamma*S**3)
      for y in yy:
        claimed=-p*S**(-s.Rational(n,2))+2*a*B4*S+beta*B6*P2-3*gamma*B6*S*S+(-2*B4+2*beta*B6*S)*y-3*B6*y*y
        assert s.simplify(s.diff(f,y)-claimed)==0
      assert F(3)<F(49,16) and 3*F(577351,10**6)**2>1
      radial=F(2,21) if n==3 else F(1,16)
      margin=radial-2*F(3,n+2)*b4*n-3*F(15,n+8)*b6*n;assert margin>0
      low,count=angular_min(n,R,b4,b6);eta=A+low-R**(2-n);assert eta>0 and R*R>n
      p4=F(n-1,n+2);p6=F(32,693) if n==3 else F(1,16)
      coremargin=core**(2-n)-A-b4*p4*core**4-b6*p6*core**6;assert coremargin>0 and core*core>n and core<R
      cornerh4=n-F(3,n+2)*n*n;cornerh6=n-F(15,n+8)*n*n+F(30,(n+8)*(n+4))*n**3
      upper=(F(577351,10**6) if n==3 else F(1,4))-A-b4*cornerh4-b6*cornerh6;assert upper>0
      mass=(4*F(157,50) if n==3 else 2*F(157,50)**2)*(core*core/2-A*core**n/n);cost=mass/(2**n*upper)
      exact=F(8332743966701654419227,1008410159020737344500) if n==3 else F(161842656400168712483969749,9292517099598284640000000)
      assert cost==exact and cost>2**n
      rows.append({'n':n,'A':str(A),'B4':str(b4),'B6':str(b6),'R':str(R),'sphere_candidates':count,'outer_negative_margin_exact':str(eta),'outer_negative_margin_display_only':float(eta),'cube_negative_derivative_margin_exact':str(margin),'cube_negative_derivative_margin_display_only':float(margin),'core_radius':str(core),'core_margin_exact':str(coremargin),'mass_lower_exact':str(mass),'corner_upper_exact':str(upper),'specific_kernel_cost_lower_exact':str(cost),'specific_kernel_cost_lower_display_only':float(cost)})
    return {'status':'PASS new U07 fixed-parameter independent audit','source_sha256':HASHES,'runs':rows,'scope':'two exact kernels only; complete cube monotonicity sufficient condition, global sphere guard and core-mass cost rejection; no C_n or whole-family lower bound','libraries':'Python int/Fraction; SymPy '+s.__version__}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--upper-dir',type=Path,default=HERE.parent/'upper');ap.add_argument('--save',action='store_true');args=ap.parse_args();out=run(args.upper_dir);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('019_u07_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
