"""One new double-Hermite candidate, exact rational root/coefficient/global checks."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,sqrt,floor
import hashlib,json,time,os,sys,importlib.util
from sympy import integer_nthroot
P=Path(__file__).resolve().parents[1];HERE=Path(__file__).resolve().parent;start=time.monotonic()
def load(name,f):
 s=importlib.util.spec_from_file_location(name,f);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def add(a,b):
 c=[F(0)]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c
def mul(a,b):
 c=[F(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def scale(a,x):return [x*y for y in a]
def evaluate(a,x):return sum((q*x**i for i,q in enumerate(a)),F(0))
def deriv(a):return [i*a[i] for i in range(1,len(a))]
def root_bounds(r):
 D=10**18;k=int(integer_nthroot((D**11*r.denominator**6)//r.numerator**6,11)[0]);lo,hi=F(k,D),F(k+1,D)
 assert lo>0 and lo**11*r**6<=1<hi**11*r**6
 return lo,hi
cornerfile=P/'astra_upper006/verification.json';corner=json.loads(cornerfile.read_text());oldfile=P/'astra_upper004/corner_checkpoint.json';old=json.loads(oldfile.read_text());eta=F(27,250)
dm={k:list(map(F,corner['moments'][str(k)])) for k in range(1,5)};bm=list(map(F,old['moments']['E_b_corner_W_power']))
mu=[sum((comb(7,n)*eta**n*dm[k][n] for n in range(8)),F(0)) for k in range(1,5)]
variance=mu[0]*mu[2]-mu[1]**2;assert variance>0
A=(mu[0]*mu[3]-mu[1]*mu[2])/variance;B=(mu[1]*mu[3]-mu[2]**2)/variance;disc=A*A-4*B;assert disc>0
rf=[(float(A)-sqrt(float(disc)))/2,(float(A)+sqrt(float(disc)))/2];nodes=[F(round(r*10000),10000) for r in rf];r1,r2=nodes;assert 0<r1<r2<1
bases=[]
for r,other in [(r1,r2),(r2,r1)]:
 L=[-other/(r-other),1/(r-other)];L2=mul(L,L);diff=[-r,F(1)]
 value=add(mul(add([F(1)],scale(diff,-2/(r-other))),L2),scale(mul(diff,L2),-F(6,11)/r))
 assert evaluate(value,r)==1 and evaluate(deriv(value),r)==-F(6,11)/r and evaluate(value,other)==0 and evaluate(deriv(value),other)==0
 bases.append(value)
intervals=[root_bounds(r) for r in nodes];ell=[F(0)]
for i in range(4):ell.append(sum((base[i]*(pair[0] if base[i]>=0 else pair[1]) for base,pair in zip(bases,intervals)),F(0)))
U=sum((ell[k]*mu[k-1] for k in range(1,5)),F(0));V=sum((comb(11,n)*eta**n*bm[n] for n in range(12)),F(0));assert U>0 and V>0
D=10000;an=int(integer_nthroot((U.numerator*V.denominator*D**4)//(U.denominator*V.numerator),4)[0]);a=F(an,D);assert a>0 and a**4<=U/V<(a+F(1,D))**4
spacefile=P/'astra_upper007/space_moments.json';families={name:[{tuple(map(int,k.split(','))):F(v) for k,v in row.items()} for row in rows] for name,rows in json.loads(spacefile.read_text())['families'].items()}
power={}
for n in range(8):
 factor=F(11,4)*a**7*comb(7,n)*eta**n
 for k in range(1,5):
  for key,v in families['d'+str(k)][n].items():power[key]=power.get(key,F(0))+factor*ell[k]*v
for n in range(12):
 factor=-F(7,4)*a**11*comb(11,n)*eta**n
 for key,v in families['b'][n].items():power[key]=power.get(key,F(0))+factor*v
power={k:v for k,v in power.items() if v};helperfile=P/'astra_upper004/verify_global.py';helper=load('locked_global_bernstein',helperfile)
assert helper.corner_value(power)==F(11,4)*a**7*U-F(7,4)*a**11*V
face,degrees,tensor,index=helper.bernstein_face(power);assert degrees==(6,)*4 and len(tensor)==2401
constfile=P/'astra_upper003/verification.json';cL=F(json.loads(constfile.read_text())['constants']['c_lower']);minimum=tensor[index];g=cL*minimum
newref=P/'astra_upper009/verification.json';gref=max(F(x['global_g_exact_rational']) for x in json.loads(newref.read_text())['candidates'])
simple=F(g.numerator*10**8//g.denominator,10**8)
out={'status':'new_complete_auxiliary_g_candidate_pending_independent_audit' if g>gref else 'single_candidate_not_improved','eta':str(eta),'a':str(a),'gauss_selection_diagnostic':{'mu':[str(x) for x in mu],'A':str(A),'B':str(B),'float_nodes':rf},'nodes':[str(x) for x in nodes],'h_node_intervals':[[str(x) for x in pair] for pair in intervals],'Hermite_h_basis_coefficients':[[str(x) for x in base] for base in bases],'ell_d_coefficients':[str(x) for x in ell],'U':str(U),'V':str(V),'symmetric_power_coefficients':{','.join(map(str,k)):str(v) for k,v in sorted(power.items())},'Bernstein_count':len(tensor),'minimum_index':index,'minimum':str(minimum),'exact_g':str(g),'simple_g':str(simple),'new_minus_ref_exact':str(g-gref),'g_display_only':float(g),'positive_global_lower':minimum>0,'strictly_improves_qualified009':g>gref,'Bernstein_sha256':hashlib.sha256(json.dumps({','.join(map(str,k)):str(v) for k,v in sorted(tensor.items())},sort_keys=True).encode()).hexdigest(),'C5_bound_changed':False,'source_sha256':{str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [cornerfile,oldfile,spacefile,helperfile,constfile,newref,Path(__file__),HERE/'hermite_linear_registration.md']},'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-start,'python':sys.executable}}
(HERE/'hermite_linear.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['gauss_selection_diagnostic','symmetric_power_coefficients','source_sha256']},indent=2))
