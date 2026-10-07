"""Exact finite LCA refinement arithmetic; not a geometric FIRST fixture."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json

def normalize(xs):
 s=sum(xs);return [F(x,s) for x in xs]
def run():
 rows=[]
 for stage,J in enumerate((2,4,8),1):
  N=2**J;delta=F(49,32768*J)
  for name in ('uniform','unequal','zeros'):
   p=normalize([1]*N if name=='uniform' else [1+(i%7) for i in range(N)])
   r=normalize([1]*N if name=='uniform' else ([2**(i%20) for i in range(N)] if name=='unequal' else [int(i%3==0) for i in range(N)]))
   light=F(0);total=F(0);retained=F(0);sum_s=F(0);channels=0;zeros=0
   for depth in range(J):
    width=2**(J-depth)
    for left in range(0,N,width):
     ids=list(range(left,left+width));s=sum((p[i] for i in ids),F(0));sum_s+=s
     half=width//2;cl=sum((p[i] for i in ids[:half]),F(0));cr=s-cl
     for i in ids:
      w=p[i]/s;t=r[i];c=t*(cr if i<left+half else cl)
      assert c<=s*t
      # Deterministic additional gate, unchanged through the split.
      c*=int((i+depth)%5!=0)
      total+=c;channels+=1;zeros+=int(t==0)
      if t<=delta*w:light+=c
      else:retained+=c
   assert sum_s==J and total==light+retained
   assert light<=delta*sum_s and total<=1
   assert 3*light<=F(49,8192)
   rows.append(dict(stage=stage,J=J,leaves=N,pattern=name,channels=channels,zero_capture_channels=zeros,delta=str(delta),total=str(total),light=str(light),retained=str(retained),sum_soft_parent_mass=str(sum_s),traffic_u_over_lambda='3',light_traffic_over_lambda=str(3*light),universal_absorption_bound='49/8192'))
 return rows
if __name__=='__main__':
 out=Path(__file__).with_name('refinement_pressure_results.json')
 if out.exists():raise SystemExit('Refuse overwrite')
 rows=run();r={'status':'PASS_EXACT_TREE_ARITHMETIC','records':len(rows),'rounds':3,'rows':rows,'seed':None,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Normalized posterior tree/LCA and gates only; no claim these rows are original cube FIRST or geom traffic.'};out.write_text(json.dumps(r,indent=2));print(json.dumps({k:v for k,v in r.items() if k!='rows'}))
