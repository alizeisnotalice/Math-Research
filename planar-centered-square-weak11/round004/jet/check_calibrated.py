"""Independent strict rational inner evaluation at one very small t; sanity check, not a derivative proof."""
import sys,json,math,os,time,importlib.util
from fractions import Fraction as F
from pathlib import Path
D=Path(__file__).resolve().parent;P=D.parents[1];sys.path.insert(0,str(D));import jet
s=importlib.util.spec_from_file_location('independent',P/'checks/lower/lower_check.py');ind=importlib.util.module_from_spec(s);s.loader.exec_module(ind)
z=json.loads((jet.D/'jet.json').read_text());t=F(1,10**9);D2=2**80;grid=math.lcm(jet.dx,D2);a=grid//jet.dx;b=grid//D2;rect=[];start=time.time()
for row in jet.Z['records']:
 box=tuple(map(int,row[:4]));xa,xb,ya,yb=box;M=F(int(row[5]),jet.dm)-F(jet.copies(box,51),5000);k=math.isqrt(M.numerator*jet.Dr**2//(4*M.denominator))-1;w=jet.copies(box,27)-jet.copies(box,34);v=F(k,jet.Dr)**2+t*w/4;assert v>0;rho=math.isqrt(v.numerator*D2**2//v.denominator)*b
 x0=max(xb*a-rho,0);x1=min(xa*a+rho,int(jet.px*grid));y0=max(yb*a-rho,0);y1=min(ya*a+rho,int(jet.py*grid))
 if x0<x1 and y0<y1:rect.append((x0,x1,y0,y1))
A=F(ind.union_area(rect),grid**2);slope=(A-F(z['area0']))/t;exact=F(z['score_dA_dt']);assert slope>0 and abs(slope-exact)<F(1,1000)
print('PID',os.getpid(),'t',t,'rectangles',len(rect),'difference quotient',float(slope),'jet',float(exact),'absolute difference < 1/1000',abs(slope-exact)<F(1,1000),'wall',time.time()-start)
