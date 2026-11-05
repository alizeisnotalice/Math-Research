from fractions import Fraction as F
from pathlib import Path
import json,math
D=Path(__file__).resolve().parent;c=json.loads((D/'transfer_34_27.json').read_text());j=json.loads((D/'transfer_34_27.check2.json').read_text());r=F(int(j['ratio_num']),int(j['ratio_den']));old=F(4629767636980314704394281556601327346619228980165481483827,2718989330218002428308996140604808231510448537600000000000);mass=sum(F(p[2]) for p in c['pts']);px,py,reach=map(F,(c['Px'],c['Py'],c['reach']));Dr=2**40; nb=(int(reach/px)+2)*(int(reach/py)+2);T=nb*mass*Dr**2/4;rho=F(math.isqrt(T.numerator//T.denominator)+1,Dr);q=math.ceil(rho/min(px,py))+1; N=10**12; finite=r*F(N-2*q,N)**2;target=F(17027535096,10**10)
print('mass',mass,'nb',nb,'rho',rho,'q',q,'eps',F(1,2*Dr),'N',N);print('new > R2',r>old,'finite > floor',finite>target,'finite=',float(finite),'exact improvement',r-old);assert nb==30 and q==4 and r>old and finite>target
