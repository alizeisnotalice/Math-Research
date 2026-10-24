"""Independent rational enclosure of public algebraic cost, without Arb/search code."""
from fractions import Fraction as Q
from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent
p=root/'upstream-v4/cert/cert_h32_a4740_R3320_3p59123.json'
d=json.loads(p.read_text());a=Q(d['alpha']);R=Q(d['R']);h=Q(d['h']);c=Q(d['c'])
s=Q(0)
for key,value in d['a'].items():
 i,j=map(int,key.split(','));mult=1 if i==j==0 else 4 if j==0 or i==j else 8
 s+=mult*Q(value)
e=2-a;den=10**60;n=e.denominator
q=R**e.numerator * den**n
N=q.numerator//q.denominator
lo,hi=0,1<<((N.bit_length()+n-1)//n)
while lo+1<hi:
 mid=(lo+hi)//2
 if mid**n<=N:lo=mid
 else:hi=mid
assert Q(lo,den)**n<=R**e.numerator<Q(lo+1,den)**n
factor=c*a/(2-a);offset=h*h*s/2
L=factor*Q(lo,den)+offset;U=factor*Q(lo+1,den)+offset
assert L<U<Q(359123,100000)
out={'certificate_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'alpha':str(a),'R':str(R),'sum_coefficients':str(s),'exact_cost':f'({factor})*({R})^({e})+({offset})','lower':str(L),'upper':str(U),'width':str(U-L),'threshold_gap_lower':str(Q(359123,100000)-U),'method':'Integer 40th root with exact Fraction endpoint-power inequalities; no Arb and no checker import.'}
(root/'exact_cost.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
