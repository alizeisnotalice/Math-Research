"""Exact all-macroscopic-length finite controls; no simulation or floating CI."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import math,json
out=Path(__file__).resolve().parent
K=4;n=3;labels=[2,0,3,1];x=[F(3,17),F(5,19),F(7,23)]
rows=[]
for R in [2,4]:
 delta=F(1,64*K*R)
 for L in [F(10,9),F(7,4),F(5,2),F(13,4),F(23,4)]:
  w=[]
  for b in range(R):
   wi=[]
   for xx in x:
    wr=[]
    for j in range(K):
     a=(F(j)+F(2*((labels[j]+b)%R)+1,2*R))/K
     weight=F(0)
     for sh in range(math.floor(xx-L/2)-1,math.ceil(xx+L/2)+2):
      weight+=max(F(0),min(a+sh+delta/2,xx+L/2)-max(a+sh-delta/2,xx-L/2))/delta/(K*L)
     wr.append(weight)
    wi.append(wr)
   w.append(wi)
  direct=sum(math.prod(w[sum(cells)%R][i][j]for i,j in enumerate(cells))for cells in product(range(K),repeat=n))
  fixed=[math.prod(sum(w[b][i],F(0))for i in range(n))for b in range(R)]
  prior=sum(fixed,F(0))/R;eps=F(R-1,2**n)
  assert abs(direct/prior-1)<=eps
  rows.append({'R':R,'L':str(L),'joint':str(direct),'prior':str(prior),
   'relative_difference':str(abs(direct/prior-1)),'universal_candidate':str(eps)})
res={'K':K,'n':n,'x':list(map(str,x)),'rows':rows,'scope':'Exact rational finite controls; general theorem separately proved.'}
(out/'rational_results.json').write_text(json.dumps(res,indent=2)+'\n')
print('Exact rational checks:',len(rows),'all passed')
