"""R060 source lock/reference check plus a new compressed six-term certificate."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,re,os,sys

ROOT=Path('/Users/zhengzhihao/Desktop/T');OUT=Path(__file__).parent
p=ROOT/'cube-cover-manuscript/cube-dim-7-bounds/paper/main.tex'
s=p.read_text();body=p.read_bytes()
assert hashlib.sha256(body).hexdigest()=='a5ffdfac3871d4ecd3d7a8e27ecb484019d62f309ce2d7a6e434754e8e17e445'
new_line=r' L\ge64,\ n\ge75,\ nd\ge1,\ t\ge4,\ tq_{\min}\ge4r_{\rm hi},\\'
old_line=r' L\ge64,\ n\ge75,\ nd\ge1,\ t\ge4,\ tq_{\min}\ge4r_{\rm hi},\quad'
assert s.count(new_line)==1
recovered=s.replace(new_line,old_line,1).encode()
assert hashlib.sha256(recovered).hexdigest()=='93b9af22489aa2d8fc55d0ab3d01373c89175aa60efd88996a35d60268215aa4'
labels=re.findall(r'\\label\{([^}]+)\}',s);refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',s)
cites=[x for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}',s) for x in group.split(',')]
bib=re.findall(r'\\bibitem\{([^}]+)\}',s)
assert not(set(refs)-set(labels)) and len(labels)==len(set(labels)) and not(set(cites)-set(bib))

# This does not rerun the source's 64-term certificate or writer's 60 checks.
d=F(391,1905);ab=F(283,2032);bb=F(63,64);e=F(1,2**20);hc=F(1,2**48)
qm=1-bb*d;qM=1-ab*d;R0=qm+d;rl=R0+e;rh=qM+d-e
def six_term_log(u):
    lo=2*sum((u**(2*j+1)/F(2*j+1) for j in range(6)),F(0))
    return lo,lo+2*u**13/(13*(1-u*u))
loA,_=six_term_log((rh-rl)/(rh+rl));_,hiQ=six_term_log(d/(4*R0-d))
P=1-e-1280*hc;ulo=P*(1-e)/(1+e);uhi=(1+e)/(P*(1-e))
M=2*uhi/qm;b=(1-F(1,4096))*2*ulo/qM-2*e
m0=F(159477,1000000);Gamma=m0*b/M;gamma0=F(1637,12500)
old_gamma=F(172531397335235056462877787662549380239,1643616601687301779182642153860479582208)
old1=F(9775559871,281956214374400);old2=F(735063024285513,21076775872652902400)
checks={
 'six_term_mass':(1-F(1,65536))*loA>m0,
 'six_term_pin':hiQ<F(26979,250000),
 'fixed_wide_budget':M>2 and b>0 and 0<gamma0<(1-e)*Gamma<Gamma<1,
 'wide_source_width':(qM-qm)/qm>d,
 'wide_net_slack':gamma0-F(3,97)-F(1,10)==F(39,1212500),
 'old_functions_root_cost_ordering':max(old1,old2)<F(1,1000) and F(1,1000)+F(4,2**20)<F(1,40),
 'positive_B_buffer_cost_ordering':old_gamma-F(7,100)+F(4,2**20)<F(1,10),
}
assert all(checks.values()),checks
paths={
 'final_manuscript':p,
 'formatting_delta':ROOT/'cube-cover-manuscript/cube-dim-7-bounds/upper/round060/manuscript_notes/formatting_delta.patch',
 'R056':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md',
 'R056_audit':ROOT/'cube-r030-uniform-audit/cube-dim-7-bounds/audit/round056/refined_inverse/review.md',
 'R059':ROOT/'cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band/PROOF.md',
 'R059_audit':ROOT/'cube-r030-soft/cube-dim-7-bounds/audit/round059/wide_source_band/review.md',
 'local_J':ROOT/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md'}
hashes={k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in paths.items()}
assert hashes['R056']=='6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a'
assert hashes['R059']=='fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f'
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,'hashes':hashes,
 'sources':{k:str(v) for k,v in paths.items()},'recovered_initial_sha256':hashlib.sha256(recovered).hexdigest(),
 'labels':len(labels),'references':len(refs),'undefined_references':sorted(set(refs)-set(labels)),
 'duplicate_labels':sorted({x for x in labels if labels.count(x)>1}),'undefined_citations':sorted(set(cites)-set(bib)),
 'scope':'New six-term positive-log compression and complete-cost/source checks only; no writer 60 checks, author 64-log baseline, old finite gates, production, MATLAB or compilation.'}
(OUT/'sourcehash.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'new exact certificates; references',len(labels),len(refs),'pid',os.getpid())
