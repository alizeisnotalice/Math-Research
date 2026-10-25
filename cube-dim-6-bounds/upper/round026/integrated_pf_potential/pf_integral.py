"""Fixed x0 only. Formal recurrence checks and certified PF640 integration."""
import hashlib,json,math,sys
from fractions import Fraction as Q
from pathlib import Path
from flint import arb,fmpq,ctx
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T')
PF=ROOT/'Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round023/partial_fraction_kernel'
F=ROOT/'Math-Research-cube6-literature/cube-dim-6-bounds/audit/round023/partial_fraction_review'
T=Q(1,25)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def B(q):
 q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def endpoints(v):
 l,e=v.lower().man_exp();u,f=v.upper().man_exp();e,f=int(e),int(f);g=min(e,f)
 return dict(lower_mantissa=str(int(l)*2**(e-g)),upper_mantissa=str(int(u)*2**(f-g)),common_binary_exponent=str(g))
def bounds(d):
 f=Q(2)**int(d['common_binary_exponent']);return int(d['lower_mantissa'])*f,int(d['upper_mantissa'])*f
def ball(d):
 l,u=bounds(d);return arb(B((l+u)/2),B((u-l)/2))
def check_inputs():
 r=json.loads((P/'registration.json').read_text())
 assert r['T']=='1/25' and r['point_numerators']==[50,50,50,45,27,0] and r['denominator']==50
 for p,h in r['source_sha256'].items():assert sha(Path(p))==h,('source changed',p)
 dep=ROOT/'Math-Research-cube6-audit/cube-dim-6-bounds/upper/round010/arb_engine/dependencies.json'
 for p,h in json.loads(dep.read_text())['runtime_library_sha256'].items():assert sha(Path(p))==h,('runtime changed',p)
 br=json.loads((PF/'results.json').read_text());fr=json.loads((F/'results.json').read_text())
 assert br['Gamma_G']==fr['Gamma_G'] and fr['Gamma_matches_B023_exactly'] and fr['all_negative_powers_cancel']
 assert br['zero_even_channels_exactly_cancelled'] and br['nonzero_real_principal_slots']==417
 return r,br
def table():
 d=json.loads((PF/'partial_tables.json').read_text());out=[]
 assert d['normalization']=='p_T=(2*pi)^-6 R(T^2)'
 for z in d['zero_pole_channels']:
  if z['order']%2==0:assert Q(z['A'])==0
 for row in d['real_z_principal']:
  s=Q(row['b_squared']);assert s>0
  n,dn=math.isqrt(s.numerator),math.isqrt(s.denominator);assert n*n==s.numerator and dn*dn==s.denominator
  b=Q(n,dn);assert b*b==s and b>0
  j=row['order'];assert isinstance(j,int) and 1<=j<=6
  out.append((s,b,j,Q(row['coefficient'])))
 assert len(out)==417 and len({s for s,b,j,d in out})==359 and max(j for s,b,j,d in out)==3
 return out
def formal_checks():
 import sympy as S
 t,b=S.symbols('t b',positive=True);I=S.atan(t/b)/b;rows=[]
 for j in range(1,5):
  if j>1:I=t/(2*b*b*(j-1)*(t*t+b*b)**(j-1))+S.Rational(2*j-3,2*(j-1))*I/(b*b)
  residual=S.factor(S.diff(I,t)-(t*t+b*b)**(-j));assert residual==0
  rows.append(dict(case=len(rows)+1,kind='formal derivative of recursive I_j',j=j,residual=str(residual)))
  limit=S.limit(I,t,0,dir='+');assert limit==0
  rows.append(dict(case=len(rows)+1,kind='formal zero-time limit',j=j,limit=str(limit)))
 assert len(rows)==8
 return dict(status='PASS',positive_time_values=0,formal_scalar_cases=8,sympy=S.__version__,cases=rows)
def evaluate():
 ctx.prec=640;_,br=check_inputs();rows=table();t=B(T);integral=arb(0);praw=fmpq(0);gamma=Q(0);caches={}
 for s,b,j,d in rows:
  if s not in caches:
   sb=B(s);ib=(t/B(b)).atan()/B(b);arr=[arb(0),ib]
   for k in range(2,4):arr.append(t/(2*sb*(k-1)*(t*t+sb)**(k-1))+(2*k-3)*arr[-1]/(2*sb*(k-1)))
   caches[s]=arr
  integral+=B(d)*caches[s][j]
  den=(T*T+s)**j;praw+=fmpq(d.numerator,d.denominator)/fmpq(den.numerator,den.denominator)
  gamma+=d*Q(math.comb(2*j-2,j-1),4**(j-1))/b**(2*j-1)
 assert gamma==Q(br['Gamma_G']) # saved-table algebra, no new tuple/kernel Green pass
 scale=64*arb.pi()**6;K=integral/scale;p=arb(praw)/scale;G=B(Q(br['Gamma_G']))/(128*arb.pi()**5)
 kd,pd=endpoints(K),endpoints(p);assert bounds(kd)[0]>0 and bounds(pd)[0]>0
 return dict(status='PASS_PF640_FIXED_POINT',precision_bits=640,T='1/25',principal_entries=417,positive_square_roots_exact=359,max_order=3,K_T=kd,p_T=pd,G_from_frozen_Gamma=endpoints(G),density_raw_rational=str(praw),density_normalization='p_T=raw_rational/(64*pi^6)',Gamma_saved_table_equality=True,spatial_gradient=False,box_oracle=False,new_C6_bound=False)
