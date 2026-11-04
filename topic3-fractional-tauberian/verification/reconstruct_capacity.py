from fractions import Fraction as F
from pathlib import Path
import json
from witness_capacity import atoms

def overlap(R,Q):
 return max(F(0),min(R[1],Q[1])-max(R[0],Q[0]))*max(F(0),min(R[3],Q[3])-max(R[2],Q[2]))
def reconstruct(row):
 al=F(row['alpha']);a=F(row['radius']);ex,ey,ed=map(F,row['scales']);ells=[ex,ex,ey,ey,ed,ed,ed,ed];C=[(a,0),(-a,0),(0,a),(0,-a),(a/2,a/2),(-a/2,a/2),(a/2,-a/2),(-a/2,-a/2)];boxes=[(x-l/2,x+l/2,y-l/2,y+l/2) for (x,y),l in zip(C,ells)];caps=atoms(boxes);m={int(s):F(v) for s,v in row['primal_masses'].items()};assert sum(m.values())<=1 and all(0<=v<=caps[s] for s,v in m.items())
 X=sorted({x for b in boxes for x in b[:2]});Y=sorted({y for b in boxes for y in b[2:]});left=dict(m);E=[]
 for x0,x1 in zip(X,X[1:]):
  for y0,y1 in zip(Y,Y[1:]):
   x=(x0+x1)/2;y=(y0+y1)/2;mask=sum(1<<i for i,(l,r,b,t) in enumerate(boxes) if l<x<r and b<y<t);take=min(left.get(mask,F(0)),(x1-x0)*(y1-y0))
   if take:E.append((x0,x0+take/(y1-y0),y0,y1));left[mask]-=take
 assert all(v==0 for v in left.values())
 fill=1-sum((r-l)*(t-b) for l,r,b,t in E)
 if fill:E.append((F(10),10+fill,F(10),F(11)))
 assert sum((r-l)*(t-b) for l,r,b,t in E)==1
 for i,R in enumerate(E):
  for Q in E[:i]:assert overlap(R,Q)==0
 captures=[sum(overlap(R,Q) for R in E) for Q in boxes];slack=min(s-al*l for s,l in zip(captures,ells));assert slack>0 and slack==F(row['primal_strict_slack'])
 return {'alpha':str(al),'radius':str(a),'scales':row['scales'],'mass':'1','rectangles':[[str(x) for x in R] for R in E],'captures':[str(x) for x in captures],'strict_slack':str(slack),'rectangle_count':len(E)}
if __name__=='__main__':
 p=Path(__file__).resolve().parents[1]/'results';data=json.loads((p/'varying_scales.json').read_text());out={'scope':'exact rational reconstruction and verification of five actual measurable indicator inputs','inputs':[reconstruct(r) for r in data['selected_certificates']+data['refined_cases']]};(p/'reconstructed_inputs.json').write_text(json.dumps(out,indent=2));print([(r['alpha'],r['rectangle_count'],float(F(r['strict_slack']))) for r in out['inputs']])
