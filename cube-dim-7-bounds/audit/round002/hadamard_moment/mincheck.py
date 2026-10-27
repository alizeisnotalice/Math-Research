#!/usr/bin/env python3
"""Independent dual-parity/Bareiss and768-bit fixed moment audit. No MC."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import hashlib,json,time,platform
import sympy as sp
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube7-moment-r002/cube-dim-7-bounds/upper/round002/hadamard_moment')
ROOT=HERE.parents[3]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def snap():return {p.name:sha(p)for p in sorted(AUTHOR.iterdir())if p.is_file()}
def sealed():return {str(p.relative_to(ROOT)):sha(p)for p in sorted((ROOT/'cube-dim-7-bounds').glob('**/round001/**/*'))if p.is_file()}
start=time.perf_counter(); before=snap(); old=sealed()
manifest=json.loads((AUTHOR/'manifest.json').read_text())
for n,v in manifest['files'].items(): assert before[n]==v['sha256'] and (AUTHOR/n).stat().st_size==v['bytes']
author=json.loads((AUTHOR/'results.json').read_text())
receipt=json.loads((AUTHOR/'receipt.json').read_text())
for n,k in [('mincheck.py','check_sha256'),('results.json','results_sha256'),('run.stdout','stdout_sha256'),('run.stderr','stderr_sha256')]:assert before[n]==receipt[k]
provenance=json.loads((AUTHOR/'sources.json').read_text())
for n,h in provenance['sha256'].items():assert sha(Path(provenance['root'])/n)==h
# Fibres of three independent dual parity checks, not author character/coset algorithm.
groups=defaultdict(list)
for mask in range(64):
    b=[(mask>>i)&1 for i in range(6)]
    syndrome=(b[0]^b[1]^b[2],b[0]^b[3]^b[4],b[1]^b[3]^b[5])
    groups[syndrome].append(mask)
cosets=sorted(groups.values(),key=lambda x:x[0])
assert len(groups)==8 and all(len(g)==8 for g in cosets)
assert sorted(sum(cosets,[]))==list(range(64)) and cosets==author['coset_masks']
# Exact fraction-free determinant; no author implementation or determinant routine called.
def bareiss(rows):
    a=[list(x)for x in rows]; n=len(a); previous=1; sign=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n)if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        p=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                val=a[i][j]*p-a[i][k]*a[k][j]
                assert val%previous==0
                a[i][j]=val//previous
            a[i][k]=0
        previous=p
    return sign*a[-1][-1]
blocks=[];determinants=[]
for masks in cosets:
    block=[[1]*8]+[[1-2*((x>>i)&1)for x in masks]for i in range(6)]
    assert [[sum(block[i][k]*block[j][k]for k in range(8))for j in range(7)]for i in range(7)]==[[8*int(i==j)for j in range(7)]for i in range(7)]
    minors=[bareiss([[row[j]for j in range(8)if j!=omit]for row in block])for omit in range(8)]
    assert all(abs(d)==512 for d in minors)
    blocks.append(block); determinants.append(minors)
assert blocks==author['blocks'] and determinants==author['signed_7_minors']
assert len({tuple(b[i][j]for i in range(7))for b in blocks for j in range(8)})==64
assert 8*512**2==8**7 and F(512**2,4**7)==16
# Exact exponent and boundary coefficient identities.
a=sp.Rational(1,112);q=sp.Rational(7,6)
assert (4+7*a,4+8*a,sp.Rational(1,2)+a)==(sp.Rational(65,16),sp.Rational(57,14),sp.Rational(57,112))
assert -sp.Rational(1,2)/8/7==-a
assert -sp.Rational(3,2)-32*a==-sp.Rational(25,14)
C=1/(2*sp.pi**sp.Rational(7,2));H=2*sp.pi**sp.Rational(7,2)/sp.gamma(sp.Rational(7,2))
assert sp.simplify(H-16*sp.pi**3/15)==0
assert sp.simplify(C**q*H/896-1/(105*2**sp.Rational(25,6)*sp.pi**sp.Rational(13,12)))==0
ctx.prec=768
def rat(n,d=1):return arb(fmpq(n,d))
# Collapsed log-Gamma expression differs from author's block-moment product path.
logMH=-rat(25,14)*rat(2).log()+8*(rat(65,16).lgamma()-rat(57,14).lgamma())+64*(rat(57,112).lgamma()-rat(1,2).lgamma())
MH=logMH.exp()
MG=(64*(rat(65,128).lgamma()-rat(1,2).lgamma())-rat(4).log()).exp()
pi=arb.pi();ca=1/(2*pi**rat(7,2));ha=2*pi**rat(7,2)/rat(7,2).gamma()
IH=(rat(7,6)*logMH).exp()*ca**rat(7,6)*ha/896
IG=MG**rat(7,6)*ca**rat(7,6)*ha/896
def endpoints(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
ints={'M_group':endpoints(MH),'M_Gamma':endpoints(MG),'I_group':endpoints(IH),'I_Gamma':endpoints(IG)}
contain={}
for k,v in ints.items():
    lo,hi=map(F,v);source=author['intervals'][k]
    alo,ahi=F(source['lower_rational']),F(source['upper_rational'])
    assert 0<alo<lo<hi<ahi
    contain[k]=True
R=F(178444228577,20000000000000000)
assert R==F(author['rigorous_rational_I_upper'])
assert F(ints['M_group'][1])<F(ints['M_Gamma'][0])
assert F(ints['I_group'][1])<F(author['rigorous_dyadic_I_upper'])<R<F(ints['I_Gamma'][0])
assert snap()==before and sealed()==old
out={'status':'PASS independent analytical and finite Hadamard-moment audit; no new C7 bound','environment':{'python':platform.python_version(),'sympy':sp.__version__,'python_flint':flint.__version__,'arb_bits':768},'partition_generation':'dual parity syndrome fibres t1+t2+t3,t1+t4+t5,t2+t4+t6 overF2','coset_masks':cosets,'all_8_Gram8I_verified':True,'all_64_directions_exactly_once':True,'signed_64_minors_Bareiss':determinants,'all_squared_minor_coefficients':512**2,'block_Y_determinant_coefficient':'16','group_negative_moment_exponent':'1/112','evaluation_path':'collapsed log-Gamma M_H; original full-boundary C^q H/896 coefficient','expression_outward_rational_intervals_768':ints,'author512_strictly_contains_independent768_intervals':contain,'I768_rational_upper':ints['I_group'][1],'author512_dyadic_I_upper':author['rigorous_dyadic_I_upper'],'author_coarse_rational_I_upper':str(R),'strict_comparison':'I_H768upper < author512dyadicupper < R_I < I_Gamma768lower; M_H768upper < M_Gamma768lower','R_I_minus_I768upper_exact':str(R-F(ints['I_group'][1])),'source_hashes':before,'sealed_round001_hashes':old,'source_hashes_before_after_equal':True,'sealed_round001_before_after_equal':True,'author_checker_imported_or_run':False,'new_MC':0,'old_baseline_runs':0,'new_Green_nodes':0,'global_g_proved':False,'new_C7_bound':False,'endpoint_tightening':'arithmetic sharpening of the same fixed analytic expression only','elapsed_seconds':time.perf_counter()-start}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ('status','environment','all_8_Gram8I_verified','all_64_directions_exactly_once','author512_strictly_contains_independent768_intervals','strict_comparison','elapsed_seconds')},indent=2))
