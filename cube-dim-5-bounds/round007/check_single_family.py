"""Independent rational reconstruction of the degree-77 family certificate.
Uses saved exact mixed moments; does not import or run the producer engine.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib,time
p=Path(__file__).resolve().parents[1];start=time.monotonic();f=p/'astra_upper003/verification.json';d=json.loads(f.read_text());mom=d['exact_mixed_moments'];cs=d['constants'];src=d['whole_linear_single_correlation_family_exclusion']
A=Q(cs['A_upper']);B=Q(cs['B_lower']);c=Q(cs['c_upper']);tau=Q(15800608,26277615)
assert A**11>Q(17,11)**11/tau**6 and 0<B and B**11<Q(6,11)**11/tau**17

def atan(q,n=40):
 L=sum((Q((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(n)),Q(0));return L,L+Q(1,(2*n+1)*q**(2*n+1))
a,b=atan(5);e,h=atan(239);piL=16*a-4*h;piH=16*b-4*e
assert piL<piH and c*c*2**15*piL**5>9
Ed=[Q(v) for v in mom['E_d_r_power']];Ed2=[Q(v) for v in mom['E_d_squared_r_power']];B11=[Q(v) for v in mom['E_B11_r_power']];B33=[Q(v) for v in mom['E_B33_r_power']];B12=[Q(v) for v in mom['E_B12_r_power']]
U=[comb(7,k)*(A*Ed[k]-B*Ed2[k]) for k in range(8)]
V=[comb(11,k)*(2*B11[k]+3*B33[k]+2*B12[k]) for k in range(12)]
assert all(U[k]==0 for k in range(1,len(U),2));assert all(V[k]<=0 for k in range(1,len(V),2))

def mul(a,b):
 r=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):r[i+j]+=x*y
 return r

def power(a,n):
 r=[Q(1)]
 for _ in range(n):r=mul(r,a)
 return r
r1=power(V,7);r2=power(U,11);g=Q(src['target_ceiling']);N=max(len(r1),len(r2))-1
R=[g**4*(r1[k] if k<len(r1) else 0)-c**4*(r2[k] if k<len(r2) else 0) for k in range(N+1)]
assert N==77 and R==list(map(Q,src['residual_power_coefficients']))
bern=[sum((R[j]*Q(comb(k,j),comb(N,j)) for j in range(k+1)),Q(0)) for k in range(N+1)]
assert all(x>0 for x in bern)
mini=min(bern);assert mini==Q(src['Bernstein_partition'][0]['minimum_Bernstein_coefficient'])
out={'status':'PASS_entire_single_correlation_linear_family_ceiling','parameter_domain':'a>=0, -1<=eta<=1','degree':N,'strict_positive_Bernstein_count':len(bern),'minimum_Bernstein_coefficient':str(mini),'minimum_indices':[i for i,v in enumerate(bern) if v==mini],'A_upper_B_lower_verified_by_integer_powers':True,'c_upper_verified_with_exact_Machin_pi':True,'negative_eta_covered_by_even_U_and_odd_V_sign':True,'scope':'full polynomial certificate reconstructed; exact mixed moments are locked author inputs, not independently regenerated','input_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start}
(p/'round007/check_single_family.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['strict_positive_Bernstein_count'],out['elapsed_seconds'])
