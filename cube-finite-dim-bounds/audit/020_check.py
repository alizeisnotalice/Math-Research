"""Round8 finite interface prechecks only: no old certificate/LP imports."""
from fractions import Fraction as F
from math import comb
from itertools import product
from pathlib import Path
import hashlib,json,argparse
import sympy as s
HERE=Path(__file__).resolve().parent

def origin_rows(n):
    # Count actual binary signed centers independently, not a stored row.
    row=[F(0)]*(n+1);hessian=[F(0)]*(n+1)
    for k in product((-1,0,1),repeat=n):
        j=sum(z!=0 for z in k);vals=[F(2,3) if z==0 else F(1,6) for z in k]
        row[j]+=__import__('math').prod(vals)
        hessian[j]+=(F(-2) if k[0]==0 else F(1))*__import__('math').prod(vals[1:])
    return row,hessian

def split(lo,hi,layer):
    L=lo[0]+layer;U=hi[0]+layer;out=[]
    for parity in (0,1):
        a=-((-(L-parity))//2);b=(U-parity)//2
        for m in product(range(a,b+1),*(range(lo[j],hi[j]+1) for j in range(1,4))):out.append((parity,m))
    return out

def run():
    rows=[]
    for n,expected in [(3,[8,12,6,1]),(4,[16,32,24,8,1])]:
        v,H=origin_rows(n);assert [int(z*3**n) for z in v]==expected
        assert all(H[j]==v[j]*(F(9*j,n)-3) for j in range(n+1))
        # Accurate eliminate all-ones coefficient (integer row's final weight1).
        aa=[F(7*j-4,11+j) for j in range(n)];aa.append(-sum(expected[j]*aa[j] for j in range(n)))
        assert sum(v[j]*aa[j] for j in range(n+1))==0
        diag=sum(H[j]*aa[j] for j in range(n+1));assert diag==F(9,n*3**n)*sum(j*expected[j]*aa[j] for j in range(n+1))
        rows.append({'n':n,'integer_origin_row':expected,'origin_denominator':3**n,'unit_grid_Hessian_diagonal_coefficients':list(map(str,H)),'general_grid_scaling':'multiply Hessian by q^2 for h=1/q','all_ones_elimination_exact_zero':True})
    bad_origin=(8*F(-1)+12*F('0.666666666667'))/27;assert bad_origin==F(1,6750000000000)>0 and bad_origin/2>0
    intervals=[(-5,-1),(-4,-1),(-3,2),(-2,2),(-1,-1),(0,0),(1,1),(1,4),(2,6)];period=[F(3,2),F(5,4),F(7,6),F(9,8)];cases=0
    for L,U in intervals:
      lo=(L,-1,0,1);hi=(U,1,1,1);other=6
      for layer in (0,1):
        direct=set(product(range(L+layer,U+layer+1),range(-1,2),range(0,2),range(1,2)))
        new=split(lo,hi,layer);recovered={((2*m[0]+p),)+m[1:] for p,m in new};assert recovered==direct and len(new)==len(direct)
        E=sum(p==0 for p,m in new);O=sum(p==1 for p,m in new);N=(U-L+1)*other;delta=E-O
        assert E+O==N and delta==(0 if (U-L+1)%2==0 else (-1)**(L+layer)*other)
        for p,m in new:
          physical=[period[0]*2*(m[0]+F(p,2))]+[period[j]*m[j] for j in range(1,4)]
          old=[period[0]*(2*m[0]+p)]+[period[j]*m[j] for j in range(1,4)];assert physical==old
        for eps in (F(-1,3),F(0),F(1,4)):
          assert 1+eps>0 and 1-eps>0 and (1+eps)+(1-eps)==2
          assert sum((1+eps if p==0 else 1-eps for p,m in new),F(0))==N+eps*delta;cases+=1
    # Explicit failure of the positive-slope Taylor bound at negative slope.
    eps=F(1,2);r0=F(1,2);tangent=r0*(1-eps/4);old_error=3*r0*eps**2/32;old_lower=tangent-old_error
    assert old_lower==F(109,256) and old_lower**4>F(1,32)
    # Safe rational replacement of (1-rho)^(-7/4) by (1-rho)^(-2).
    rho=F(1,2);B=3*r0/(32*(1-rho)**2)
    for delta in (-1,1):
      for eps in (F(-1,2),F(-1,4),F(0),F(1,4),F(1,2)):
        tan=r0*(1+F(delta,4)*eps);lower=tan-B*eps**2;upper=tan+B*eps**2;mass_radius_fourth=(1+eps*delta)/16
        assert lower>=0 and lower**4<=mass_radius_fourth<=upper**4
    return {'status':'PASS finite independent round8 interface prechecks','origin':rows,'independent_rounding_error_origin':str(bad_origin),'parity_mass_cases':cases,'integer_blocks':len(intervals),'receiver_layers':[0,1],'supercell_mass_exact':'2','negative_slope_counterexample':{'N':1,'delta':-1,'epsilon':'1/2','true_radius_fourth':'1/32','incorrect_positive_slope_lower':'109/256','incorrect_lower_fourth_exceeds_true':True},'safe_Taylor_condition':'|epsilon*delta/N|<=rho<1, remainder <=(3*r0/32)*(delta/N)^2*(1-rho)^(-7/4)*epsilon^2; rational factor (1-rho)^(-2) also sufficient','scope':'new finite identities and parity embedding only; no old R4 volume or LP run; author source audit pending','libraries':'Python int/Fraction; SymPy '+s.__version__}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--save',action='store_true');args=ap.parse_args();out=run();out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('020_precheck_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
