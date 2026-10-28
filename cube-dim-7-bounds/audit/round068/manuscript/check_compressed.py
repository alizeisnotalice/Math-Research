"""R068 new compressed threshold checks, preservation and source/reference locks."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json,re,os,sys

ROOT=Path('/Users/zhengzhihao/Desktop/T');OUT=Path(__file__).parent
writer=ROOT/'cube-cover-manuscript/cube-dim-7-bounds'
p=writer/'paper/main.tex';body=p.read_bytes();text=body.decode()
assert hashlib.sha256(body).hexdigest()=='8135adfb51470228ed18371c6c5f1a22caac0b35523b0cd49c50b13cfb9b386b'
row_new=r'其 $\lambda$ 行列积分均至多 $D_j>0$，行积分在 $A_j$ 上 $\lambda$ 几乎处处还至少为 $r_j\ge0$。'
row_old=r'其 $\lambda$ 行列积分均至多 $D_j>0$，行积分还至少为 $r_j\ge0$。'
tau_new=r'若其 $\mu$ 行积分在 $\mu$ 几乎处处的 $A_j$ 上至多为同一个 $\tau\ge0$，则'
tau_old=r'若其 $\mu$ 行积分在 $\mu$ 几乎处处的 $A_j$ 上至多为同一个 $\tau$，则'
assert text.count(row_new)==text.count(tau_new)==1
stage=text.replace(row_new,row_old,1)
assert hashlib.sha256(stage.encode()).hexdigest()=='4c2bed1c4e227d7caeea511c0739a2457a01ca947cf35853f91a3ee9654bc0f4'
initial=stage.replace(tau_new,tau_old,1)
assert hashlib.sha256(initial.encode()).hexdigest()=='1b9e8be8a90316d44fc10e878218859bfb86dedba90beaa448c6bf27458238fb'
base=(writer/'upper/round068/manuscript_notes/private/main_R060.tex').read_text()
def block(s,start,end=None):
    a=s.index(start);return s[a:] if end is None else s[a:s.index(end,a)]
core=r'\section{预备知识}';refined=r'\subsection{固定参数、完整来源带与捕获核}'
oldwide=r'\subsection{扩大来源带后的另一组固定参数}'
newprofile=r'\subsection{乘法匹配、径向权与相对局部域}'
low=r'\subsection{较小起效门槛与可同时使用的常数族}'
assert block(base,core,refined)==block(text,core,refined)
assert block(base,refined,oldwide)==block(text,refined,newprofile)
assert block(base,low)==block(text,low)
labels=re.findall(r'\\label\{([^}]+)\}',text);refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text)
cites=[item for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}',text) for item in group.split(',')]
bib=re.findall(r'\\bibitem\{([^}]+)\}',text)
assert not(set(refs)-set(labels)) and len(labels)==len(set(labels)) and not(set(cites)-set(bib))

# A new weaker n lower bound than both author N=L^8/8! and prior soft L^8/2^16.
L=F(2)**30;N=L**8/2**18;sqrtL=F(2)**15
e=F(1,2**20);ea=F(1,65536);eret=F(1,4096);hc=F(1,2**80)
qm=F(2,3);qM=F(917,1000);sig=F(83,1000);chi=F(9,200);x=(1-sig)*L
h=hc/L;delta=F(1,2**30);rl=1+delta;rh=F(2751,2000)-delta
eta=(1-qM)/(2*qM);Csrc=(qM-qm)/qm;cS=F(917,10**7);nu=F(47,200000)
mean=(F(3,4)*L**2/qm**2+F(1,2)*L/qm)/N
var=F(3,2)*L**2/(qm**2*N)
bad_log=var/(e-mean)**2;bad_count=4/(e*e*N);bad_cap=F(3,4)**64
bad_profile=F(200000)**2/(N*qm*qm)
S11=factorial(11)*Csrc*rh/(hc*cS**11*L**10)
T11=factorial(11)*rh/(hc*nu**11*L**9*sqrtL)
R3=2/L**2+factorial(3)/(sig*L)**3+2*sqrtL*factorial(3)/(sig*L/2)**3
checks={
 'weaker_n_bound':factorial(8)<2**18 and N>=400 and N/L>16,
 'extended_local_geometry':h<min(delta/2,F(1,1000)) and qm*N/L>4*rh and 1/(x*N*N)<F(1,2) and qM/(x*N*N)<2*h/3,
 'individual_count_profile_bridges':L/(3*qm*N)<e and 2*h+1/N<3*e/4 and h<F(1,1600000),
 'baseline_pin_mean':mean<e/2 and N*qm>=200000,
 'all_five_events':max(bad_log,bad_count,bad_cap,bad_profile)<ea/8 and bad_log+2*bad_count+bad_cap+bad_profile<ea,
 'profile_probability_identity':F(200000)**2/(qm*qm)==90000000000,
 'per_face_retention':F(5000)**2/(150*x)<eret and x>=1000,
 'radial_and_cap_domain':1/(1-(h/rl)**2)<1+e and 1/rh-L/N>F(1,2),
 'finite_n_cap_rate':eta/400<(eta-chi)/2 and eta-chi-eta/400>F(1,10000) and Csrc*rh<1,
 'different_global_cap_tail':S11<e,
 'different_global_tau_tail':T11<e,
 'printed_Taylor10_coarse_margins':F(factorial(10),2**50)<e and F(factorial(10),2**44)<e,
 'printed_base_arguments':cS*L>2**16 and nu*L>2**17 and rh*L*sqrtL/hc<2**126,
 'B_positive_and_new_root_tail':L/4-4/L-F(317,1000)>0 and R3<F(17,50000),
 'price_slack':F(317,1000)-F(4333,50000)-F(23,100)==F(17,50000),
 'positive_old_B_price':F(61,200)-F(91,1000)+F(4,2**20)<F(23,100),
}
assert all(checks.values()),checks
paths={
 'final_manuscript':p,
 'tau_delta':writer/'upper/round068/manuscript_notes/directsum_scope_delta.patch',
 'row_delta':writer/'upper/round068/manuscript_notes/row_domain_delta.patch',
 'accepted_R060_snapshot':writer/'upper/round068/manuscript_notes/private/main_R060.tex',
 'R066':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/PROOF.md',
 'R066_audit':ROOT/'cube-r030-uniform-audit/cube-dim-7-bounds/audit/round066/target_profile/review.md',
 'R069':ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round069/profile_threshold/THRESHOLD_PROOF.md',
 'R069_audit':ROOT/'cube-r030-soft/cube-dim-7-bounds/audit/round069/profile_threshold/REVIEW.md'}
hashes={key:hashlib.sha256(path.read_bytes()).hexdigest() for key,path in paths.items()}
assert hashes['R069']=='d1e0420ba6ce509d2a1e9d5e0c6c01fd8aedd62c52adc3bc187b584fd646da10'
assert hashes['R069_audit']=='ebe763640e003dd87b168217be6d84a6f0bcd0dbee20bb5b38bd99483874f43a'
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,'hashes':hashes,'sources':{k:str(v) for k,v in paths.items()},
 'labels':len(labels),'references':len(refs),'undefined_refs':sorted(set(refs)-set(labels)),'duplicate_labels':[],
 'preserved_blocks':['original cover/cap/actual measure','complete R056 proof','R034/R037/bibliography tail'],
 'recovered_stages':['4c2bed1c4e227d7caeea511c0739a2457a01ca947cf35853f91a3ee9654bc0f4','1b9e8be8a90316d44fc10e878218859bfb86dedba90beaa448c6bf27458238fb'],
 'new_exact_bounds':{k:str(v) for k,v in {'target_union':bad_log+2*bad_count+bad_cap+bad_profile,'S11':S11,'T11':T11,'R3':R3}.items()},
 'scope':'New weaker-n/direct-union/Taylor11/Taylor3 compression plus references and byte/source locks. No author23/41, soft46, old60/128, log baseline, production, MATLAB, Git or compilation.'}
(OUT/'sourcehash.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'new compressed exact checks; labels/refs',len(labels),len(refs),'pid',os.getpid())
