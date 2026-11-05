import json,sys,math,time,os
from pathlib import Path
from fractions import Fraction as F
D=Path(__file__).resolve().parent;P=D.parents[1];sys.path.insert(0,str(D));import jet
start=time.time();n=0
for row in jet.Z['records']:
 xa,xb,ya,yb=map(int,row[:4]);box=(xa,xb,ya,yb)
 for idx in (27,34,51):
  x,y,_=jet.pts[idx-1];x=int(x*jet.dx);y=int(y*jet.dx)
  exact=sum(xa<=x+i*jet.pxi<=xb and ya<=y+j*jet.pyi<=yb for i in range(-4,5) for j in range(-4,5))
  assert exact==jet.copies(box,idx);n+=1
print('PID',os.getpid(),'period-copy counts verified',n,'wall',time.time()-start)
# No scan outside [-4,4] can intersect the selected boxes.
assert all(abs(int(row[0])) < 4*jet.pxi and abs(int(row[1])) < 4*jet.pxi and abs(int(row[2])) < 4*jet.pyi and abs(int(row[3])) < 4*jet.pyi for row in jet.Z['records'])
