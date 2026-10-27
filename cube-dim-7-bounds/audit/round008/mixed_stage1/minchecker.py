#!/usr/bin/env python3
"""Read-only stdout-only StageI mixed Mellin audit. No panel generation/StageII."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import argparse,hashlib,json,time,platform
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent;CUBE=HERE.parents[2]
p=argparse.ArgumentParser();p.add_argument('--numeric-dir',type=Path,default=CUBE/'upper/round008/mixed_stage1');p.add_argument('--design',type=Path,default=CUBE/'upper/round007/theory/mixed_moment_design.md');p.add_argument('--r004-dir',type=Path,default=CUBE/'upper/round004/pair_mellin');a=p.parse_args()
NUM=a.numeric_dir.resolve();DESIGN=a.design.resolve();OLD=a.r004_dir.resolve()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=[x for x in sorted(NUM.iterdir())if x.is_file()]+[DESIGN,OLD/'certificate.json',OLD/'manifest.json']
def snap():return {str(x):sha(x)for x in files}
start=time.perf_counter();before=snap()
for n,h in json.loads((NUM/'manifest.json').read_text())['files'].items():assert sha(NUM/n)==h,n
raw=(NUM/'r004-certificate.json').read_bytes();old=json.loads(raw)
assert sha(NUM/'r004-certificate.json')==json.loads((NUM/'r004-manifest.json').read_text())['files']['certificate.json']
assert raw==(OLD/'certificate.json').read_bytes()
assert sha(OLD/'certificate.json')==json.loads((OLD/'manifest.json').read_text())['files']['certificate.json']
C=json.loads((NUM/'certificate.json').read_text());assert C['reused_q']==old['q']
def dy(d):return Q(int(d[0]))*Q(2)**int(d[1])
def box(d):return(dy(d['lower']),dy(d['upper']))
def contains(d,iv):assert dy(d['lower'])<=iv[0]<=iv[1]<=dy(d['upper'])
def A(x):x=x if isinstance(x,Q)else Q(str(x));return arb(fmpq(x.numerator,x.denominator))
def ep(x):return(Q(str(x.lower().fmpq())),Q(str(x.upper().fmpq())))
def txt(iv):return[str(x)for x in iv]
P=[json.loads(x)for x in (NUM/'panels.jsonl').read_text().splitlines()]
assert len(P)==256 and {(x['group'],x['panel'])for x in P}=={(j,k)for j in range(4)for k in range(64)} and all(x['status']=='completed'for x in P)
negative_panels=0;negative_terms=0;nonreal_constants=0;term_sum_contained=0
for row in P:
    l,u=box(row['integral']);assert l<=u
    if u<0:negative_panels+=1
    terms=row['signed_integrated_even_terms'];assert len(terms)==65
    tl=Q(0);tu=Q(0)
    for t in terms:
        v,w=box(t);assert v<=w;tl+=v;tu+=w
        if w<0:negative_terms+=1
    assert max(l,tl)<=min(u,tu)
    term_sum_contained+=int(l<=tl<=tu<=u)
    il,iu=box(row['constant_imag']);assert il<=iu
    nonreal_constants+=int(not il<=0<=iu)
ctx.prec=768
e=Q(1,7);b=Q(9,8);ar=b/2;rho=Q(27,64);h=Q(1,4);eta=h/rho
assert ar-rho>max(Q(0),7*e-Q(7,2)) and 64*2*h==32 and eta==Q(16,27)
# Independent radial/simplex parameter check for all selected counts at all envelope nodes.
for j in range(8):
    for v in (ar-rho,ar,ar+rho):
        assert Q(1,2)+v/7-e>0 and Q(1,2)+v/7>0 and 4+v-j*e>0 and 4+8*v/7-j*e>0
assert [comb(7,j)for j in range(4)]==[1,7,21,35] and sum(comb(7,j)for j in range(8))==128
# Positive real weighted Mellin formula via log-Gamma, independently arranged.
def F(v,j):return (-(A(Q(4,7))*v+A(2*j*e))*A(2).log()+(4+v-A(j*e)).lgamma()-(4+A(Q(8,7))*v-A(j*e)).lgamma()+j*(A(Q(1,2))+v/7-A(e)).lgamma()+(8-j)*(A(Q(1,2))+v/7).lgamma()-8*A(Q(1,2)).lgamma()).exp()
a0=A(ar);rr=A(rho);norm=arb.pi()*A(b).gamma();reassemblies=[]
for j in range(4):
    d=C['mixed_expectations'][j];assert d['j']==j
    L=Q(0);U=Q(0)
    for row in P:
        if row['group']==j:l,u=box(row['integral']);L+=l;U+=u
    contains(d['polynomial_lower'],(L,L));contains(d['polynomial_upper'],(U,U))
    cross1=(a0-rr).gamma()*F(a0-rr,j)*(a0+rr).gamma()*F(a0+rr,7-j)
    cross2=(a0+rr).gamma()*F(a0+rr,j)*(a0-rr).gamma()*F(a0-rr,7-j)
    contains(d['cross1'],ep(cross1));contains(d['cross2'],ep(cross2))
    cs=box(d['cstar_upper']);assert cs[0]==cs[1] and cs[0]>=max(ep(cross1)[1],ep(cross2)[1])
    E=32*A(cs[0])*A(eta)**129/(1-A(eta));tail=F(a0,j)*F(a0,7-j)*(2*(a0+16).lgamma()).exp()/(norm*31*A(32)**31)
    contains(d['raw_error'],ep(E));contains(d['normalised_tail_abs'],ep(tail))
    lo=((A(L)-E)/norm-tail).lower();hi=((A(U)+E)/norm+tail).upper()
    ilo=max(Q(0),Q(str(lo.fmpq())));ihi=Q(str(hi.fmpq()));assert 0<ilo<ihi
    declared_lo=box(d['expectation_lower']);declared_hi=box(d['expectation_upper'])
    assert declared_lo[0]==declared_lo[1]<=ilo<ihi<=declared_hi[0]==declared_hi[1]
    assert d['lower_clipped_to_zero'] is False
    reassemblies.append({'j':j,'exact_signed_polynomial_sum':txt((L,U)),'cross1':txt(ep(cross1)),'cross2':txt(ep(cross2)),'cstar_declared_upper':str(cs[0]),'raw_error':txt(ep(E)),'signed_tail_absolute_bound':txt(ep(tail)),'mixed_expectation_interval':txt((ilo,ihi))})
# Exact declared expectation endpoints, and separate monotone LH endpoints.
weights=[1,7,21,35]
Slo=2*sum(weights[j]*dy(C['mixed_expectations'][j]['expectation_lower']['lower'])**7 for j in range(4))
Shi=2*sum(weights[j]*dy(C['mixed_expectations'][j]['expectation_upper']['upper'])**7 for j in range(4))
qlo,qhi=box(old['q']);assert 0<qlo<qhi
LHlo=(8*(A(Slo).log()/7).exp()/A(qhi)).lower();LHhi=(8*(A(Shi).log()/7).exp()/A(qlo)).upper()
CL,CU=dy(C['LH_lower']['lower']),dy(C['LH_upper']['upper'])
assert CL<=Q(str(LHlo.fmpq()))<Q(str(LHhi.fmpq()))<=CU<1
assert box(C['LH_lower'])[0]==box(C['LH_lower'])[1] and box(C['LH_upper'])[0]==box(C['LH_upper'])[1]
assert box(C['mu_lower'])==(CL-1,CL-1) and C['mu_positive'] is False and C['status']=='mu_not_positive'
assert C['stage_II_authorized'] is False
# This is arithmetic validation of an out-of-domain *formal* expression, not a lower bound.
formal=A(dy(old['M_pair']['lower']))/A(CU).sqrt();stored_formal=box(C['optimistic_surrogate_M_lower'])
assert stored_formal[0]==stored_formal[1]<=ep(formal)[0]
target=Q(31,500);assert Q(C['target'])==target
# Feasible lower certificate for mu can only be0 here. Its certified correction is0.
assert CL-1<0 and CU-1<0
pair_budget_lower=dy(old['M_pair']['lower']);assert pair_budget_lower>target
assert snap()==before
out={'audit_status':'ACCEPTED StageI failed-positivity certificate; no StageII or new I','environment':{'python':platform.python_version(),'python_flint':flint.__version__,'bits':768},'new_panels_recomputed':0,'unique_signed_panels':256,'negative_panel_intervals':negative_panels,'negative_integrated_even_term_intervals':negative_terms,'nonzero_imaginary_constant_intervals':nonreal_constants,'stored_panel_contains_exact_sum_of_term_intervals_count':term_sum_contained,'R004_certificate_verbatim_reuse':True,'four_mixed_scalar_reassemblies':reassemblies,'LH_independent_interval':[str(LHlo.fmpq()),str(LHhi.fmpq())],'declared_LH_interval':txt((CL,CU)),'exact_LH_upper_below_one':True,'positive_mu_gate':False,'feasible_mu_lower_from_this_surrogate':'0','StageII_authorized':False,'formal_out_of_domain_expression_only':txt(ep(formal)),'formal_field_is_not_valid_optimistic_lower':True,'certificate_capability_only_pair_budget_lower':str(pair_budget_lower),'material_target':str(target),'materiality_interpretation':'fixed surrogate cannot clear mu>0; formal negative-mu expression has no lower-bound semantics; mu=0 gives no certifiable correction; no actual M lower or C7 exclusion','source_hashes':before,'source_before_after_equal':True,'author_checker_executed':False,'new_MC':0,'new_I_budget':False,'global_g_proved':False,'new_C7_bound':False,'elapsed_seconds':time.perf_counter()-start}
print(json.dumps(out,indent=2))
