from pathlib import Path
import json,numpy as np
P=Path(__file__).resolve().parent;rows=[]
for r in json.loads((P/'refinement.json').read_text()):
 b=np.array(r['coefficients']);m=np.arange(1,len(b)+1);a=10/9;g=27/29
 bound=(2*np.pi**2*np.sum(m*m*b)/a**2+4*np.pi*np.sum(m*b)/a**3+6*(1+b.sum())/a**4)/g
 gridmax=0
 for l,u in [(10/9,11/9),(1.3,1.4)]:
  x=np.linspace(l,u,20001);vals=[]
  for batch in np.array_split(x,50):vals.extend(abs((1+np.sin(np.pi*batch[:,None]*m)**2@b)/(g*batch*batch)-1))
  gridmax=max(gridmax,max(vals))
 step=(1/9)/20000
 rows.append(dict(M=len(b),raw_min=r['rawMinimum'],nonnegative=bool(b.min()>=0),mass=float(b.sum()),grid_error=float(gridmax),second_derivative_bound=float(bound),interpolation_allowance=float(bound*step**2/8),continuous_estimate=float(gridmax+bound*step**2/8)))
(P/'check.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
