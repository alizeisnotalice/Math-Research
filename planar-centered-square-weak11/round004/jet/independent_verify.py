"""Independent compression-grid coefficient check; no author's segment tree.
Global periodic mass counting reconstructed directly from the R004a input.
"""
from pathlib import Path
from fractions import Fraction as F
import json,math,time,os,hashlib,argparse
import numpy as np
S=Path(__file__).resolve().parent;P=S.parents[1];parser=argparse.ArgumentParser();parser.add_argument('--records-dir',type=Path,default=P/'round003/mechanism/generated');parser.add_argument('--output-dir',type=Path,default=S/'generated');args=parser.parse_args();D=args.output_dir.resolve();D.mkdir(parents=True,exist_ok=True);start=time.time()
z=json.loads((args.records_dir.resolve()/'decomposition.json').read_text())
c=json.loads((P/'round004/weights/candidate_05.json').read_text())
dx,dm,L,dr=(int(z[k]) for k in ('Dx','Dm','L','Dr'));sx=L//dx;sr=L//dr
px,py=(F(c[k]) for k in ('Px','Py'));pxi,pyi=int(px*dx),int(py*dx)
pts=[tuple(map(F,p)) for p in c['pts']];pts=[(int(x*dx),int(y*dx),int(m*dm)) for x,y,m in pts]
def number(a,b,x,p):return max(0,(b-x)//p+((x-a)//p)+1)
R=[];gaps=[];bounds=[];clip_ties=0
for row in z['records']:
 a,b,d,e=map(int,row[:4]);counts=[number(a,b,x,pxi)*number(d,e,y,pyi) for x,y,m in pts]
 mass=sum(n*m for n,(x,y,m) in zip(counts,pts));assert mass>0
 k=math.isqrt(mass*dr*dr//(4*dm))-1;assert k>0
 w=counts[26]-counts[33];v=F(w,k);rho=k*sr
 raw=((b*sx-rho,-v),(a*sx+rho,v),(e*sx-rho,-v),(d*sx+rho,v))
 clip=((0,F(0)),(int(px*L),F(0)),(0,F(0)),(int(py*L),F(0)))
 clip_ties+=sum(q[0]==h[0] for q,h in zip(raw,clip))
 rr=(max(raw[0],clip[0]),min(raw[1],clip[1]),max(raw[2],clip[2]),min(raw[3],clip[3]))
 if rr[0]<rr[1] and rr[2]<rr[3]:R.append(rr)
 gap=F(mass,dm)-4*F(k,dr)**2;assert gap>0;gaps.append(gap)
 if w<0:bounds.append(4*F(k,dr)**2/(-w))
X=sorted({v for r in R for v in r[:2]});Y=sorted({v for r in R for v in r[2:]});xi={v:i for i,v in enumerate(X)};yi={v:i for i,v in enumerate(Y)}
inds=[(xi[r[0]],xi[r[1]],yi[r[2]],yi[r[3]]) for r in R];area0=0;area1=F(0)
for s in range(0,len(X)-1,128):
 end=min(s+128,len(X)-1);diff=np.zeros((end-s+1,len(Y)),dtype=np.int32)
 for a,b,c,d in inds:
  if b<=s or a>=end:continue
  aa,bb=max(a,s)-s,min(b,end)-s
  diff[aa,c]+=1;diff[aa,d]-=1;diff[bb,c]-=1;diff[bb,d]+=1
 covered=np.cumsum(np.cumsum(diff,axis=0),axis=1)[:end-s,:len(Y)-1]>0
 for i,cover in enumerate(covered):
  transitions=np.diff(np.concatenate(([0],cover.view(np.int8),[0])))
  starts=np.flatnonzero(transitions==1);ends=np.flatnonzero(transitions==-1)
  h0=sum(Y[b][0]-Y[a][0] for a,b in zip(starts,ends));h1=sum((Y[b][1]-Y[a][1] for a,b in zip(starts,ends)),F(0))
  w0=X[s+i+1][0]-X[s+i][0];w1=X[s+i+1][1]-X[s+i][1]
  area0+=w0*h0;area1+=w0*h1+w1*h0
area=F(area0,L*L);score=area1*dr/(8*L)
expected=json.loads((S/'jet.json').read_text())
assert area==F(expected['area0']) and score==F(expected['score_dA_dt']) and score>0
assert min(gaps)==F(expected['minimum_strict_mass_gap']) and min(bounds)==F(expected['positive_radicand_until'])
out={'pid':os.getpid(),'seconds':time.time()-start,'independent_method':'infinite periodic mass counts + coordinate compression 2D coverage, no source segment tree','records':len(R),'x_jet_lines':len(X),'y_jet_lines':len(Y),'constant_coordinate_clip_ties':clip_ties,'area0':str(area),'score':str(score),'score_decimal':float(score),'minimum_strict_gap':str(min(gaps)),'positive_radicand_until':str(min(bounds)),'source_jet_exact_match':True,'source_jet_sha256':hashlib.sha256((S/'jet.py').read_bytes()).hexdigest(),'no_finite_area_batch_rerun':True}
(D/'independent_verify.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
