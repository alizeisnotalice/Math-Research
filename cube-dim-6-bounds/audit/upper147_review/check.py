from pathlib import Path
from fractions import Fraction as F
import os,json,hashlib,ast,time
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
P=Path(__file__).resolve().parent
S=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round147/pair_pushforward')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
start=time.monotonic();cpu=time.process_time()
manifest=read(S/'manifest.json');assert sha(S/'manifest.json')=='2500cf19538335300f9c4f2db74c914620acba279af8e182a1f43c96eb9de69f'
for item in manifest['files']:
 p=S/item['path'];assert sha(p)==item['sha256'] and p.stat().st_size==item['bytes']
packet=read(S/'next_cases.json');reg=read(S/'next_registration.json');rel=read(S/'next_release_requirements.json')
assert sha(S/'generate_next.py')==reg['run_sha256']==rel['run_sha256']
assert sha(S/'next_cases.json')==reg['case_packet_sha256']==rel['case_packet_sha256']
assert sha(S/'PROOF.md')==reg['proof_sha256']==rel['proof_sha256']
assert sha(S/'next_registration.json')==rel['registration_sha256']
for p,h in packet['source_hashes'].items():assert sha(p)==h
cones=read(packet['cone_table_path'])['cones'];details=[]
for row,m,d in zip(packet['cases'],[3,6],[3,2]):
 c=cones[row['cone_id']];i,j=row['body_pair_zero_based'];cols=list(zip(c['nodes_by_body'][i],c['nodes_by_body'][j]));assert [list(v) for v in cols]==row['all_six_columns'];v=[tuple(map(F,x)) for x in cols];assert all(x>=0 for col in v for x in col);active=[col for col in v if any(col)]
 dirs={tuple(x/sum(col) for x in col) for col in active};assert len(active)==m and len(dirs)==d and row['alpha']==c['alpha'];assert any(a[0]*b[1]-a[1]*b[0] for a in active for b in active)
 details.append({'case_id':row['case_id'],'active_columns':m,'directions':d,'rank':2,'source_exact':True})
for n in ['generate_next.py','prototype.py']:ast.parse((S/n).read_text())
unrun=['next_parent_release.json','next_started.json','next_result.json','next_exit_receipt.json'];assert not any((S/n).exists() for n in unrun)
import sympy as sp
assert sp.__version__=='1.14.0'
a,b,x,s,t=sp.symbols('a b x s t',positive=True)
M=sp.Matrix([[sp.Rational(1,60),sp.Rational(1,30)],[sp.Rational(1,60),0]])
inv=M.inv()*sp.Matrix([a,b]);assert inv==sp.Matrix([60*b,30*(a-b)])
density=sp.exp(-sum(inv))/abs(M.det());assert sp.simplify(density-1800*sp.exp(-30*a-30*b))==0
# a=b+x makes the triangle the positive quadrant, Jacobian 1.
f=sp.simplify(density.subs(a,b+x))
def integral(g):return sp.integrate(sp.integrate(sp.expand(g.subs(a,b+x))*f,(x,0,sp.oo)),(b,0,sp.oo))
mass=integral(sp.S(1));EA=integral(a);EB=integral(b);EAB=integral(a*b);cov=sp.simplify(EAB-EA*EB)
assert (mass,EA,EB,EAB,cov)==(1,sp.Rational(1,20),sp.Rational(1,60),sp.Rational(1,900),sp.Rational(1,3600))
LT=integral(sp.exp(-s*a-t*b));target=1/((1+s/30)*(1+(s+t)/60));assert sp.cancel(LT-target)==0
r=sp.symbols('r');eta=sp.symbols('eta',nonzero=True)
for n in range(5):
 Q=sum((-1)**k*sp.diff(r**n,r,k)/eta**(k+1) for k in range(n+1));assert sp.simplify(sp.diff(Q,r)+eta*Q-r**n)==0
 assert sp.diff(r**(n+1)/sp.Integer(n+1),r)==r**n
# Pure algebraic endpoint-order lemma, independent of source density generation.
L1,L2,V1,V2,Z1,Z2=sp.symbols('L1 L2 V1 V2 Z1 Z2')
# Equality normal N=L/L(v)-K/K(v) satisfies N.v=0 identically.
K1,K2=sp.symbols('K1 K2');den1=L1*V1+L2*V2;den2=K1*V1+K2*V2
assert sp.cancel((L1/den1-K1/den2)*V1+(L2/den1-K2/den2)*V2)==0
for item in manifest['files']:assert sha(S/item['path'])==item['sha256']
assert not any((S/n).exists() for n in unrun)
out={'status':'PASS_INDEPENDENT_PROOF_SOURCE_PREFLIGHT_AND_ONE_M2_EXACT_CHECK','registration_sha256':sha(P/'registration.json'),'checker_sha256':sha(__file__),'source_manifest_sha256':sha(S/'manifest.json'),'source_runner_sha256':sha(S/'generate_next.py'),'source_proof_sha256':sha(S/'PROOF.md'),'source_cases':details,'sympy_version':sp.__version__,'m2':{'density':'1800*exp(-30*A-30*B), A>B>0','mass':str(mass),'mean_A':str(EA),'mean_B':str(EB),'mixed_moment':str(EAB),'covariance':str(cov),'Laplace':str(sp.factor(LT))},'generic_antiderivative_degrees':list(range(5)),'endpoint_order_normal_annihilates_new_ray':True,'m3_m6_densities_generated':0,'k_h_or_special_function_values':0,'old_6D_grid_calls':0,'native_calls':0,'new_bound_claimed':False,'nice':os.getpriority(os.PRIO_PROCESS,0),'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
(P/'actual_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
