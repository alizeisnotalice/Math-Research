"""Exact finite-frame spectral geometry, not a kernel-cost calculation."""
from pathlib import Path
from itertools import product
import hashlib,json
import sympy as s

def run():
 rows=[]
 for n in (3,4):
  V=s.Matrix([(1,)+v for v in product((-1,1),repeat=n-1)])
  assert V.T*V==2**(n-1)*s.eye(n)
  if n==3:
   B=s.Matrix([[1,1,1],[1,1,-1],[1,-1,1]]).T
  else:
   B=s.Matrix([[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]]).T
  det=abs(int(B.det()));assert det==(4 if n==3 else 16)
  checks=0
  for x in product(range(-3,4),repeat=n):
   ordered=sorted(map(abs,x),reverse=True);a,b,c=ordered[:3]
   exact=sum(abs(int(v.dot(s.Matrix(x)))) for v in [V.row(i) for i in range(V.rows)])
   predicted=max(4*a,2*(a+b+c)) if n==3 else max(8*a,4*(a+b+c),2*(3*a+b+c+ordered[3]))
   assert exact==predicted,(n,x,exact,predicted)
   checks+=1
  # Chamber identities: nonnegative terms plus the remaining signed term(s).
  a,b,c,d=s.symbols('a b c d',nonnegative=True)
  if n==3:
   first=(a+b+c)+(a+b-c)+(a-b+c)
   assert s.expand(first+(a-b-c)-4*a)==0
   assert s.expand(first-(a-b-c)-2*(a+b+c))==0
  else:
   first=(a+b+c+d)+(a+b+c-d)+(a+b-c+d)+(a+b-c-d)+(a-b+c+d)+(a-b+c-d)
   assert s.expand(first-(6*a+2*b+2*c))==0
   assert s.expand(first+2*(a-b-c)-8*a)==0
   assert s.expand(first-2*(a-b-c)-4*(a+b+c))==0
   assert s.expand(first+2*d-2*(3*a+b+c+d))==0
  if n==3:
   volume=s.Rational(1,2)**3-8*s.Rational(1,4)**3/s.factorial(3)
   assert volume==s.Rational(5,48)
   assert s.factorial(3)*volume/8==s.Rational(5,64)
  rows.append({'n':n,
               'Fourier_unit_ball_volume':str(volume) if n==3 else 'not computed',
               'exact_density_peak':'5/(64*pi^3*t^3)' if n==3 else 'not computed','directions':[list(map(int,V.row(i))) for i in range(V.rows)],'Gram_diagonal':2**(n-1),'basis_columns':[[int(B[i,j]) for i in range(n)] for j in range(n)],'absolute_basis_determinant':det,'density_sup_sufficient_bound':'1/('+str(det)+'*pi^'+str(n)+'*t^'+str(n)+')','ordered_chamber_symbolic_checks':'PASS','integer_points_diagnostic_only':checks,'psi_formula':'max(4a,2(a+b+c))' if n==3 else 'max(8a,4(a+b+c),2(3a+b+c+d))'})
 return {'status':'PASS exact spectral geometry and spanning; not kernel positivity or cost certification','rows':rows,'sympy':s.__version__,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':
 out=run();p=Path(__file__).resolve().parents[1]/'results/round008_frame_interface.json';p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
