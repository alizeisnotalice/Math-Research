import sympy as s,json,time
from pathlib import Path
p=Path(__file__).resolve().parent
a=s.symbols('a',positive=True);L=1+7*s.log(a);N=a-7*s.log(a)
A=(1+2/a)**7;B=(1+1/a)**7*L
H=s.sqrt(A)+s.sqrt(B)
R=s.Rational(1,2)*((a+1)/(a+2))**s.Rational(5,2)*N/s.sqrt(L)
P=7*(1+2/a)**s.Rational(5,2)/a**2
t=time.process_time()
c1=s.simplify(s.diff(H,a)-P*(R-1))
c2=s.simplify(s.diff(N/s.sqrt(L),a)*L**s.Rational(3,2)-(7*s.log(a)-s.Rational(5,2)-7/a-s.Rational(49,2)*s.log(a)/a))
assert c1==c2==0
print(json.dumps({'status':'PASS two symbolic derivative identities','identity1':str(c1),'identity2':str(c2),'scope':'algebra only; interval domain and global monotonicity separately proved','sympy':s.__version__,'cpu_seconds':time.process_time()-t},indent=2))
