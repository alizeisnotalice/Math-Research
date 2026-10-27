#!/usr/bin/env python3
"""Read-only R004 audit: ledger+three fixed scalar assemblies, no panel rerun."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time,platform
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent
NUM=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube7-mellin-r004/cube-dim-7-bounds/upper/round004/pair_mellin')
SOL=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube7-strategy-r004/cube-dim-7-bounds/upper/round004/strategy')
ROOT=HERE.parents[3]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
FILES=[SOL/n for n in ('pair_proof.md','pair_registration.md','theory.md','manifest.sha256')]+[p for p in sorted(NUM.iterdir())if p.is_file()]
def snapshot():return {str(p):sha(p)for p in FILES}
start=time.perf_counter();before=snapshot()
for name,h in json.loads((NUM/'manifest.json').read_text())['files'].items():assert sha(NUM/name)==h,name
# Lock all earlier sealed files without executing old checkers or science.
old={str(p.relative_to(ROOT)):sha(p)for p in sorted((ROOT/'cube-dim-7-bounds').glob('**/round00[123]/**/*'))if p.is_file()}
panels=[json.loads(line)for line in (NUM/'panels.jsonl').read_text().splitlines()]
assert len(panels)==512 and {p['panel']for p in panels}==set(range(512)) and all(p['status']=='completed'for p in panels)
def dy(p):return Q(int(p[0]))*Q(2)**int(p[1])
L=Q(0);U=Q(0)
for p in panels:
    lo,hi=dy(p['integral']['lower']),dy(p['integral']['upper'])
    assert 0<lo<=hi
    L+=lo;U+=hi
cert=json.loads((NUM/'certificate.json').read_text())
assert dy(cert['polynomial_lower']['lower'])<=L<=dy(cert['polynomial_lower']['upper'])
assert dy(cert['polynomial_upper']['lower'])<=U<=dy(cert['polynomial_upper']['upper'])
ctx.prec=768
def A(q):q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def ep(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def m(v):return (-A(Q(4,7))*v*A(2).log()+(4+v).lgamma()-(4+A(Q(8,7))*v).lgamma()+8*((A(Q(1,2))+v/7).lgamma()-A(Q(1,2)).lgamma())).exp()
a=A(Q(1,16));r=A(Q(3,64));s=A(Q(1,8));norm=arb.pi()*s.gamma()
# Scalar assembly1: entire complex-strip bound and analytical recurrence tail.
star=((a-r).lgamma()+(a+r).lgamma()).exp()*m(a-r)*m(a+r)
err=32*star*A(Q(2,3))**129*3
tail=m(a)**2*(a+16).lgamma().__mul__(2).exp()/(norm*31*A(32)**31)
# Scalar assembly2: use exact ledger sum, not any delivered sum or float value.
qlo=((A(L)-err)/norm).lower();qhi=((A(U)+err)/norm+tail).upper()
assert qlo>0 and qhi>qlo
C=1/(2*arb.pi()**A(Q(7,2)));H=2*arb.pi()**A(Q(7,2))/A(Q(7,2)).gamma()
Mlo=qlo**4/2;Mhi=qhi**4/2
Ilo=Mlo**A(Q(7,6))*C**A(Q(7,6))*H/896
Ihi=Mhi**A(Q(7,6))*C**A(Q(7,6))*H/896
final_intervals={'q':[str(qlo.fmpq()),str(qhi.fmpq())], 'M_pair':[ep(Mlo)[0],ep(Mhi)[1]],'I_pair':[ep(Ilo)[0],ep(Ihi)[1]]}
for key,x in [('M_star',star),('quadrature_absolute_error',err),('q_tail_upper_bound',tail)]:
    final_intervals[key]=ep(x)
for key,v in final_intervals.items():
    ownlo,ownhi=map(Q,v)
    assert dy(cert[key]['lower'])<=ownlo<ownhi<=dy(cert[key]['upper']),key
RI=Q(cert['rational_I_upper_numerator'],cert['rational_I_upper_denominator'])
assert RI==Q(8179950070308,10**18)
# Scalar assembly3: fixed exact R002 expression comparison only.
logMH=-A(Q(25,14))*A(2).log()+8*(A(Q(65,16)).lgamma()-A(Q(57,14)).lgamma())+64*(A(Q(57,112)).lgamma()-A(Q(1,2)).lgamma())
MH=logMH.exp();IH=(A(Q(7,6))*logMH).exp()*C**A(Q(7,6))*H/896
assert Q(final_intervals['M_pair'][1])<Q(ep(MH)[0])
assert Q(final_intervals['I_pair'][1])<RI<Q(ep(IH)[0])
assert snapshot()==before
for p,h in old.items():assert sha(ROOT/p)==h
out={'status':'PASS independent theory-interface audit; no repeated panel computation','environment':{'python':platform.python_version(),'python_flint':flint.__version__,'bits':768},'scalar_assemblies':3,'ledger_unique_completed_panels':512,'exact_polynomial_sum_interval':[str(L),str(U)],'scalar_enclosures':dict(final_intervals,R002_M_H=ep(MH),R002_I_H=ep(IH)),'independent768_contained_in_Terra_certificate':True,'rational_complete_I_upper':str(RI),'strict_comparison':'I768upper < 8179950070308/10^18 < R002_I_H768lower; Mpair768upper < R002_M_H768lower','rational_upper_minus_I768upper_exact':str(RI-Q(final_intervals['I_pair'][1])),'source_hashes':before,'source_before_after_equal':True,'earlier_sealed_files_count':len(old),'earlier_sealed_files_before_after_equal':True,'frozen_full_theory_matches_prerun':before[str(SOL/'theory.md')]=='1c6f392fe79d54ab9d749be8914dda7ba48b7e3b51e613ac985c99d98f4c906c','final_compact_proof_hash':before[str(SOL/'pair_proof.md')],'Terra_postrun_compact_proof_hash':'b5124749ddce68051b2a88b45e64e29ba2112de18281618e8083e71216f78867','postrun_compact_proof_has_later_analytical_appendix':True,'panels_recomputed':0,'author_checker_executed':False,'new_MC':0,'new_Green_nodes':0,'global_g_proved':False,'new_C7_bound':False,'elapsed_seconds':time.perf_counter()-start}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ('status','environment','scalar_assemblies','ledger_unique_completed_panels','independent768_contained_in_Terra_certificate','strict_comparison','elapsed_seconds')},indent=2))
