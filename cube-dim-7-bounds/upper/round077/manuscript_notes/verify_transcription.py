"""New manuscript transcription check; no search and no producer program execution."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json,platform,re
HERE=Path(__file__).resolve().parent
PAPER=HERE.parents[2]/'paper/main.tex'
L=F(8192); s=F(3,2); qmin=F(2,3); qmax=F(11569,12800)
sigma=F(25,256); chi=F(63,1280); x=(1-sigma)*L
N=L**64/factorial(64); h=1/(768*L); delta=F(1,2**19)
rlo=1+delta;rhi=s*qmax-delta
eB=eN=eT=F(1,2**19);eJ=eR=eS=eD=F(1,2**20);ea=F(1,65536)
Ret=F(9717,10000); mass=F(3043,10000); gamma=F(291,1000)
P=1-eB-F(6,768); Ulo=P*(1-eN)/(1+eJ);Uhi=(1+eN)/(P*(1-eJ))
M=2*(1+eR)*Uhi;b=2*Ret*(1-eR)*Ulo-eS-eT;Gamma=mass*b/M

def H(v):
 z=(v-1)/(v+1)
 lo=sum((2*z**(2*j+1)/F(2*j+1) for j in range(64)),F(0))
 return lo,lo+2*z**129/(129*(1-z*z))
def S(v):
 return sum((v**j/F(factorial(j)) for j in range(65)),F(0))
ql,qh=H(F(6,5)); mu2lo=1-5*qh-F(5,2)*qh**2;mu2hi=1-5*ql-F(5,2)*ql**2
ML=F(3,4)*L**2/(qmin**2*N)+F(1,2)*L/(qmin*N)
VL=F(3,2)*L**2/(qmin**2*N)
eta=(1-qmax)/(2*qmax);cS=(eta-chi)*(1-sigma);cT=sigma/2-chi*(1-sigma)
mu=F(2211,50000)/(1-eJ);p=F(183,1000);vv=F(531,100000);u=(chi-mu)*x-p
G={
'fixed_parameter_geometry':2*h<delta and h<=F(1,2**20) and qmax/(x*N*N)<=h/s and 1/(x*N*N)<F(1,2) and N*qmin/x>4*rhi and N/x>rhi and x/N<1 and x*N*N>1 and 2<=x<=L,
'individual_J_and_n':N>=9000 and x/N<=F(1,1000) and x/(3*qmin*N)<=eJ<=F(1,1024),
'baseline_log_and_counts':ML<=eB/2 and 4*VL/eB**2<=ea/8 and 4/(eN**2*N)<=ea/8 and 2*h+1/N<=3*eN/4,
'target_cap_probability':1+x/(2*N)-1/(2*rhi)<=F(3,4) and F(3,4)**74<=ea/8,
'first_profile_direction_pin_probability':F(1,2)-F(5,2)*ql<F(221,5000) and N>=9000 and F(324000000)/N<=ea/8,
'second_profile_finite_h_mean':0<mu2lo<mu2hi and (1-1/N)*(mu2lo-h*(mu2hi+F(1,25)))>F(528,100000) and mu2hi+1/(250*N*qmin)<F(5295,1000000),
'profile_transport_and_probability':h<=F(1,100000) and h/5<F(1,1000000) and F(5760000)/N<=ea/8,
'face_variance_and_pin':F(5301,1000000)/(1-eJ)<vv and qh<p,
'Cantelli_gap_and_retention':chi>mu and u>0 and Ret<u*u/(vv*x+u*u),
'radial_mass_and_deficit':1/(1-(h/rlo)**2)<=1+eR and (1-ea)*H(rhi/rlo)[0]>mass and 0<P<1 and M>2 and b>0 and 0<gamma<(1-eD)*Gamma<Gamma<1,
'distinct_exponents':eta==F(1231,23138) and cS==F(13627383,3790929920) and cT==F(1447,327680) and cS!=cT and cS*L>4 and cT*L>4 and eta*x/N<=F(1,1025),
'complete_prefix_payments':F(1025,1024)*H(qmax/qmin)[1]*rhi*768*L/S(cS*L)<=eS and rhi*768*L*91/S(cT*L)<=eT,
'positivity_price_and_scalar':(L/2)**2/2>=L and L/4>gamma and H(F(256,231))[1]<F(257,2500) and sigma*L/2>4 and 1/N+1/S(sigma*L)+182/S(sigma*L/2)<F(1,5000) and gamma-F(257,2500)-F(47,250)==F(1,5000),
'complete_function_comparisons':F(172531397335235056462877787662549380239,1643616601687301779182642153860479582208)-F(7,100)+F(4,2**20)<F(47,250) and gamma-sigma+F(4,2**30)<F(23,100),
}
assert all(G.values()),G
text=PAPER.read_text();labels=re.findall(r'\\label\{([^}]+)\}',text);refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text)
assert len(labels)==len(set(labels)) and set(refs)<=set(labels)
assert not re.search(r'round0|SHA256|[0-9a-f]{64}',text)
receipt={'status':'PASS','scope':'New manuscript transcription only: exact grouped base-formula checks and complete-function comparisons; no producer rerun. Analytic all-tail and operator transfer remain a written-proof task.','python':platform.python_version(),'groups':G,'group_count':len(G),'Gamma':str(Gamma),'tex_sha256':hashlib.sha256(PAPER.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'labels':len(labels),'references':len(refs),'search':False,'MATLAB':False,'production':False}
(HERE/'transcription_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'PASS','groups':len(G),'labels':len(labels),'references':len(refs),'source_sha256':receipt['tex_sha256']}))
