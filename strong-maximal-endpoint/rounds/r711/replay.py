"""Solver-free Fraction replay of saved full-product and layer certificates."""
from fractions import Fraction as F
from pathlib import Path
import json
from product_lp import star_catalogue,build_full_lp,build_layer_lp
HERE=Path(__file__).resolve().parent
def verify(A,b,c,cert):
 x=list(map(F,cert['primal_rational']));y=list(map(F,cert['dual_rational']))
 assert len(x)==len(c) and len(y)==len(b)
 assert all(0<=v<=1 for v in x) and all(v>=0 for v in y)
 assert all(sum((a*v for a,v in zip(row,x)),F(0))>=rhs for row,rhs in zip(A,b))
 primal=sum((a*v for a,v in zip(c,x)),F(0))
 dual=sum((rhs*v for rhs,v in zip(b,y)),F(0))+sum((min(F(0),c[j]-sum((row[j]*v for row,v in zip(A,y)),F(0))) for j in range(len(c))),F(0))
 assert [dual,primal]==list(map(F,cert['objective_interval']))
 assert 0<dual<=primal
 return dual,primal
if __name__=='__main__':
 data=json.loads((HERE/'product_lp.json').read_text());count=0
 for rec in data['cases']:
  boxes=[[[F(z) for z in axis] for axis in box] for box in rec['base_boxes']];k=rec['k'];h=rec['h'];a=F(rec['a']);p=F(rec['p'])
  _,_,A,b,c=build_full_lp(boxes,star_catalogue(k,h,a),p);full=verify(A,b,c,rec['full'])
  *_,A,b,c=build_layer_lp(boxes,k,h,a,p);layer=verify(A,b,c,rec['layers'])
  assert max(full[0],layer[0])<=min(full[1],layer[1]);count+=1
 out={'status':'verified','cases':count,'exact_certificates':2*count,'method':'Fraction; geometry reconstructed; no optimizer invoked'}
 (HERE/'replay.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
