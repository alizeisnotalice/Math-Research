import math,json
from pathlib import Path
from fractions import Fraction as F
D=Path(__file__).resolve().parent; c=json.loads((D/'candidate002.json').read_text()); z=json.loads((D/'candidate002.check2.json').read_text()); mass=sum(F(p[2]) for p in c['pts']); ratio=F(int(z['ratio_num']),int(z['ratio_den']));px,py,reach=map(F,(c['Px'],c['Py'],c['reach']));Dr=2**40
nb=(int(reach/px)+2)*(int(reach/py)+2);b=nb*mass*Dr**2/4;k=math.isqrt(b.numerator//b.denominator)+1; r=F(k,Dr);q=math.ceil(r/min(px,py))+1
old=F(482167717377842991882283179477574801392641890014291050323,283174368139841919615520431313000857449005056000000000000); threshold=F(17027531463,10**10); print('nb',nb,'rho upper',r,'q',q,'epsilon',F(1,2*Dr),'mass',mass);assert nb==30 and q==4
for N in (10**8,10**12):
 finite=ratio*F(N-2*q,N)**2;print('N',N,'finite ratio',float(finite),'exceeds old',finite>old,'exceeds 1.7027531463',finite>threshold);assert finite>old
 if N==10**12: assert finite>threshold
