"""Exact lexicographic first jets of a union of calibrated rectangles.
Coordinates: L-units; zeta=t*L*Dr/8; output derivative wrt t is
(area coefficient in zeta)*Dr/(8L). Fixed list only, not full ideal halo.
"""
import json,math,time,os,hashlib,argparse
from fractions import Fraction as F
from pathlib import Path
S=Path(__file__).resolve().parent;P=S.parents[1];parser=argparse.ArgumentParser();parser.add_argument('--records-dir',type=Path,default=P/'round003/mechanism/generated');parser.add_argument('--output-dir',type=Path,default=S/'generated');args=parser.parse_args();D=args.output_dir.resolve();D.mkdir(parents=True,exist_ok=True);Z=json.loads((args.records_dir.resolve()/'decomposition.json').read_text());cfg=json.loads((P/'round002/verification/candidate002.json').read_text());L=int(Z['L']);Dr=int(Z['Dr']);dx=int(Z['Dx']);dm=int(Z['Dm']);px=F(cfg['Px']);py=F(cfg['Py']);pxi=int(px*dx);pyi=int(py*dx);pts=[tuple(map(F,p)) for p in cfg['pts']];sx=L//dx;sr=L//Dr

def j(c,s=0):return (int(c),F(s))
def subtract(a,b):return a[0]-b[0],a[1]-b[1]
def count(a,b,p,period):
 return max(0,(b-p)//period - (-((p-a)//period)) + 1)  # floor((b-p)/period)-ceil((a-p)/period)+1
# Explicit form above: ceil((a-p)/period) = -floor((p-a)/period).
def copies(box,idx):
 xa,xb,ya,yb=box;x,y,_=pts[idx-1];x=int(x*dx);y=int(y*dx);return count(xa,xb,x,pxi)*count(ya,yb,y,pyi)

def union(rects):
 if not rects:return F(0),F(0)
 Y=sorted({q for r in rects for q in (r[2],r[3])}); yi={v:i for i,v in enumerate(Y)};n=len(Y)-1;cover=[0]*(4*n+8);size=[(0,F(0)) for _ in cover]
 def update(node,l,r,a,b,d):
  if b<=l or r<=a:return
  if a<=l and r<=b:cover[node]+=d
  else:
   mid=(l+r)//2;update(2*node,l,mid,a,b,d);update(2*node+1,mid,r,a,b,d)
  if cover[node]>0:size[node]=subtract(Y[r],Y[l])
  elif r-l==1:size[node]=(0,F(0))
  else:
   x,y=size[2*node],size[2*node+1];size[node]=(x[0]+y[0],x[1]+y[1])
 ev=sorted((x,delta,yi[y0],yi[y1]) for x,delta,y0,y1 in [(r[0],1,r[2],r[3]) for r in rects]+[(r[1],-1,r[2],r[3]) for r in rects]);last=ev[0][0];area0=0;area1=F(0)
 for x,d,a,b in ev:
  w=subtract(x,last);h=size[1];area0+=w[0]*h[0];area1+=w[0]*h[1]+w[1]*h[0];last=x
  if n>0:update(1,0,n,a,b,d)
 return F(area0),area1

def tests():
 # Adjacent rectangles initially cover [0,2]x[0,1]; both exposed seam edges move apart.
 seam=[(j(0),j(1,-1),j(0),j(1)),(j(1,1),j(2),j(0),j(1))]; assert union(seam)==(2,-2)
 # Initially zero-width rectangle opens to order zeta, contributing exactly 2 zeta.
 line=[(j(0,-1),j(0,1),j(0),j(1))];assert union(line)==(0,2)
 # Two singleton squares of mass 4 each decreasing at rate one: area 8-2t+O(t²).
 actual=[(j(-1,F(1,8)),j(1,-F(1,8)),j(-1,F(1,8)),j(1,-F(1,8))),(j(1,F(1,8)),j(3,-F(1,8)),j(-1,F(1,8)),j(1,-F(1,8)))]
 assert union(actual)==(8,-2)
 print('selftests seam-gap, zero-width line, two mass-4 singleton receivers PASS',flush=True)
if __name__=='__main__':
 start=time.time();tests();rects=[];gaps=[];roots=[];stats={'records':0,'active':0,'signed_owner_boxes':0,'declipped':0}
 for row in Z['records']:
  stats['records']+=1;box=tuple(map(int,row[:4]));xa,xb,ya,yb=box;m=F(int(row[5]),dm)-F(copies(box,51),5000);assert m>0
  k=math.isqrt(m.numerator*Dr**2//(4*m.denominator))-1
  if k<=0:continue
  w=copies(box,27)-copies(box,34);gap=m-4*F(k,Dr)**2;assert gap>0;gaps.append(gap)
  if w<0:roots.append(4*F(k,Dr)**2/(-w))
  vel=F(w,k);r=k*sr
  x0=max(j(xb*sx-r,-vel),j(0));x1=min(j(xa*sx+r,vel),j(int(px*L)));y0=max(j(yb*sx-r,-vel),j(0));y1=min(j(ya*sx+r,vel),j(int(py*L)))
  if x0<x1 and y0<y1:
   rects.append((x0,x1,y0,y1));stats['active']+=1;stats['signed_owner_boxes']+=bool(w)
 print('pid',os.getpid(),'stats',stats,'start union',time.time()-start,flush=True);a0,a1=union(rects);area=a0/(L*L);dA=a1*Dr/(8*L);expected=json.loads((P/'round004/weights/candidate_05.check2.json').read_text());ratio=F(int(expected['ratio_num']),int(expected['ratio_den']));mcell=sum(p[2] for p in pts)-F(1,5000);assert area==ratio*mcell
 out=dict(pid=os.getpid(),stats=stats,area0=str(area),base_area_matches_R004a=True,score_dA_dt=str(dA),positive=dA>0,score_decimal=float(dA),ratio_derivative=str(dA/mcell),ratio_derivative_decimal=float(dA/mcell),elapsed=time.time()-start,minimum_strict_mass_gap=str(min(gaps)),positive_radicand_until=str(min(roots)),scope='signed 34->27 derivative of strict calibrated finite subfamily, not ideal full halo; reconstructed from R3 positive records')
 (D/'jet.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
