"""Independent integer-grid inner/outer union check for the satellite certificate.
Uses SymPy integer_nthroot, independent x-slab/2D difference grid union.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np
from sympy import integer_nthroot
ROOT=Path(__file__).resolve().parent
D=10**24;G=10**12

def flo(x):return x.numerator//x.denominator

def cei(x):return -((-x.numerator)//x.denominator)

def build(inp):
 p=list(map(F,inp['periods']));ph=[list(map(F,q)) for q in inp['phases']];w=list(map(F,inp['weights']));inner=[];outer=[]
 for key in inp['witness_keys']:
  lower=[key[:3],key[3:6]];upper=[key[6:9],key[9:12]];count=[]
  for b in range(2):
   n=1
   for j in range(3):n*=max(0,upper[b][j]-lower[b][j]+1)
   count.append(n)
  active=[b for b in range(2) if count[b]];mass=sum((count[b]*w[b] for b in active),F(0));q=int(integer_nthroot(mass.numerator*D**3//mass.denominator,3)[0]);assert q**3*mass.denominator<=mass.numerator*D**3<(q+1)**3*mass.denominator
  for target,r,lf,hf in [(inner,F(q,2*D),cei,flo),(outer,F(q+1,2*D),flo,cei)]:
   box=[]
   for j in range(3):
    lo=max(F(0),max((upper[b][j]+ph[b][j])*p[j] for b in active)-r);hi=min(p[j],min((lower[b][j]+ph[b][j])*p[j] for b in active)+r);box.extend([lf(lo*G),hf(hi*G)])
   if all(box[2*j]<box[2*j+1] for j in range(3)):target.append(tuple(box))
 return inner,outer,sum(w)

def union_integer(boxes):
 cuts=[sorted({b[2*j+k] for b in boxes for k in range(2)}) for j in range(3)];lookup=[{x:i for i,x in enumerate(c)} for c in cuts];dy=np.diff(cuts[1]);dz=np.diff(cuts[2]);vol=0
 for x0,x1 in zip(cuts[0],cuts[0][1:]):
  diff=np.zeros((len(cuts[1]),len(cuts[2])),dtype=np.int32)
  for b in boxes:
   if b[0]<=x0 and x1<=b[1]:
    y0,y1=lookup[1][b[2]],lookup[1][b[3]];z0,z1=lookup[2][b[4]],lookup[2][b[5]];diff[y0,z0]+=1;diff[y1,z0]-=1;diff[y0,z1]-=1;diff[y1,z1]+=1
  occupied=diff.cumsum(axis=0).cumsum(axis=1)[:-1,:-1]>0;area=0
  for y in range(len(dy)):area+=int(dy[y])*int(dz[occupied[y]].sum())
  vol+=int(x1-x0)*area
 return vol,[len(c) for c in cuts]

def main():
 start=time.monotonic();fp=ROOT/'round004_certificate_input.json';inp=json.loads(fp.read_text());inner,outer,mass=build(inp);vi,ci=union_integer(inner);vo,co=union_integer(outer);low=F(vi,G**3)/mass;high=F(vo,G**3)/mass
 original=json.loads((ROOT.parents[1]/'baseline/original_certificates.json').read_text())[0];target=F(original['lower_bound_interval'][1]);assert low>target
 rec={'status':'independent exact integer-grid enclosure of saved finite union','input_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'root_denominator':D,'grid_denominator':G,'ratio_interval':[str(low),str(high)],'decimal_display_only':[float(low),float(high)],'lower_above_original_R3_upper':low>target,'rigorous_improvement_lower':str(low-target),'inner_boxes':len(inner),'outer_boxes':len(outer),'inner_cuts':ci,'outer_cuts':co,'elapsed_seconds':time.monotonic()-start}
 (ROOT/'round004_grid_result.json').write_text(json.dumps(rec,indent=2)+'\n');print(rec)
if __name__=='__main__':main()
