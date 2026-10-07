"""Exact finite checks of a general source-cell covering implication.
Each record has its stated threshold; not a common-threshold extremal example.
"""
from fractions import Fraction as F
from pathlib import Path
import random,json,hashlib
rows=[]
for stage,n in enumerate((8,32,128),1):
 rng=random.Random(202610070900+stage)
 N=8*stage
 for shape in ('one_cell','several_cells','spread'):
  atoms=[]
  for j in range(N):
   if shape=='one_cell': z=[F(rng.randrange(1,100),100*n) for k in range(n)]
   elif shape=='several_cells': z=[F((j%4)*100+rng.randrange(1,100),100*n) for k in range(n)]
   else:z=[F(rng.randrange(-150,151),100) for k in range(n)]
   atoms.append(z)
  raw=[rng.randrange(1,11) for j in range(N)];w=[F(v,sum(raw)) for v in raw]
  cells=[tuple((v*n).__floor__() for v in z) for z in atoms]
  mass={}
  for cell,weight in zip(cells,w):mass[cell]=mass.get(cell,F(0))+weight
  for row in range(3):
   R=(F(1),F(3,2),F(2))[row]
   x=[v+F(rng.randrange(-10,11),100)*R for v in atoms[row]]
   captured=[max(abs(a-b) for a,b in zip(z,x))<=R/2 for z in atoms]
   m=sum((weight for weight,hit in zip(w,captured) if hit),F(0));assert m>0
   lam=m/(3*R**n) # This record's hard band: u=3lambda.
   caps={}
   for cell,weight,hit in zip(cells,w,captured):
    if hit:caps[cell]=caps.get(cell,F(0))+weight
   for eta in (F(1,2),F(1,4)):
    dominant=[cell for cell,v in caps.items() if v>=eta*m]
    for cell in dominant:
     M=mass[cell];center=[F(2*k+1,2*n) for k in cell]
     assert R**n < M/(2*eta*lam)
     d=max(2*abs(a-b) for a,b in zip(x,center))
     assert d<=R+F(1,n)
     assert d**n <= (1+F(1,n))**n*M/(2*eta*lam)
    rows.append(dict(n=n,N=N,shape=shape,row=row,R=str(R),eta=str(eta),captured_cells=len(caps),dominant_cells=len(dominant),captured_mass=str(m),threshold=str(lam)))
  assert sum(mass.values())==1
p=Path(__file__);out=p.with_name(p.stem+'_results.json');assert not out.exists()
out.write_text(json.dumps(dict(status='PASS_EXACT_CELL_COVER_IMPLICATION',rows=rows,rounds=3,scope='Per-record rational geometric/mass implication for chosen cube, not a whole maximal/FIRST experiment. Analytic proof handles common threshold and arbitrary finite measure.',script_sha256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2)+'\n')
print(len(rows),'exact cell geometry records passed')
