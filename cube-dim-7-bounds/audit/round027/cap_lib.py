from pathlib import Path
import importlib.util,math
spec=importlib.util.spec_from_file_location('rational_lib',Path(__file__).resolve().parents[2]/'audit/round024/rational_lib.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
Q=base.Q
def coefficients(t):
 t=Q(t);q=1+1/(2*t);h=(1+1/t)/2
 terms=[Q(math.comb(6,j),j)*q**(6-j)*Q(-1,2)**j*(1-t**(-j))for j in range(1,7)]
 return 1-h**7-Q(7,2)*sum(terms,Q()),Q(7)-Q(7,2)*q**6,terms
def cap(t):
 t=Q(t);c,k,terms=coefficients(t);l,u=base.log(t)
 B=(c+k*l,c+k*u)if k>=0 else(c+k*u,c+k*l)
 assert B[0]>0 and B[1]<1+7*l
 F=base.sumroot((1+1/t)**7,B[0],B[1])
 return B,F,c,k,terms
