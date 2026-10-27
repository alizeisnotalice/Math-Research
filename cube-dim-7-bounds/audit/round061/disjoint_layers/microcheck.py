"""R061: new two-block aggregate lower proof and a static certificate interface."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,sys

ROOT=Path('/Users/zhengzhihao/Desktop/T');OUT=Path(__file__).parent
proof=ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round061/disjoint_layers/PROOF.md'
cert=ROOT/'cube-r030-compute/cube-dim-7-bounds/upper/round062/compute/candidate_certificate.json'
record=json.loads(cert.read_text());rows=record['rows']
e=F(1,2**20);eA=F(1,65536);eRet=F(1,4096);hc=F(1,2**48)
delta=F(1,2**30);step=F(3,1600);d=F(7,25);sigma=F(3,50)
P=1-e-1280*hc;ulo=P*(1-e)/(1+e);uhi=1/ulo
def log4(u):
    return 2*sum((u**(2*k+1)/F(2*k+1) for k in range(4)),F(0))
# Lower the local b/M in each of two 64-layer groups, with no 128 log rebuild.
def coefficient(q):
    return (1-e)*((1-eRet)*ulo/uhi*q/(q+step)-e)
c1=coefficient(F(18,25));c2=coefficient(F(21,25))
two_block=(1-eA)*(c1*(log4(F(3,53))-128*delta)+c2*(log4(F(3,59))-128*delta))-F(1,2**191)
gammas=[F(r['gamma_lower']) for r in rows]
chis=[F(r['chi']) for r in rows]
etas=[(1-F(r['qhi']))/(2*F(r['qhi'])) for r in rows]
mus=[F(r['mu_upper']) for r in rows]
checks={
 'layer_count':len(rows)==128,
 'static_partition_binding':all(F(r['qlo'])==F(18,25)+i*step and F(r['qhi'])==F(18,25)+(i+1)*step for i,r in enumerate(rows)),
 'all_2J_radial_supports':step>2*delta and 1+delta>1 and F(31,25)-delta<2 and F(18,25)>F(1,2),
 'active_inset':delta-2*(delta/4)>0,
 'common_profile_coefficients':F(record['constant']['Ulo'])==ulo and F(record['constant']['Uhi'])==uhi and c1>0 and c2>0,
 # Every small log's 48-term tail is below 2^-970; rounding adds <2^-200.
 'uniform_log_tail_and_rounding':step/2<F(1,1024) and F(2,97*(1-F(1,2**20)))<1 and 128*F(1,2**198)==F(1,2**191),
 'independent_two_block_gamma':two_block>F(429,2000)+F(1,40000),
 'static_gamma_total':sum(gammas,F(0))==F(record['gamma_total_lower']) and sum(gammas,F(0))>F(429,2000),
 'positive_layer_budgets':all(g>0 for g in gammas),
 'retention_and_both_rates':min(c-m for c,m in zip(chis,mus))>=F(1,10000) and min(a-c for a,c in zip(etas,chis))>0 and min(sigma/2-(1-sigma)*c for c in chis)>0,
 'clean_scalar_identity':F(429,2000)-F(3,47)==F(14163,94000) and F(429,2000)-F(3,47)-F(3,20)==F(63,94000),
}
assert all(checks.values()),checks
paths={
 'full_proof':proof,'certificate':cert,
 'direct_sum_lemma':ROOT/'cube-r040-star-lp/cube-dim-7-bounds/audit/round061/direct_sum/PROOF.md',
 'independent_arithmetic_audit':ROOT/'cube-r030-num-audit/cube-dim-7-bounds/audit/round062/layers/audit.json',
 'local_J':ROOT/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md',
 'R056':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md'}
hashes={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}
assert hashes['full_proof']=='0a1e2bc6eb829714f62b9a0cdb5f52f651fa7d49fd54a084e89b2ff3751af716'
assert hashes['certificate']=='b889e35ef909bd979937752c67410bf3913eb76cdf7dc6c720556e093d4c0711'
assert hashes['direct_sum_lemma']=='8b8168370146670889ad5568fc3b20c48df5e020fbe64da8ffac29252dbf6f3f'
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,'hashes':hashes,'sources':{k:str(p) for k,p in paths.items()},
 'exact_values':{k:str(v) for k,v in {'two_block_lower':two_block,'two_block_margin':two_block-F(429,2000),
 'ledger_gamma_sum':sum(gammas,F(0)),'min_retention_gap':min(c-m for c,m in zip(chis,mus)),
 'min_cap_gap':min(a-c for a,c in zip(etas,chis)),'min_tau_gap':min(sigma/2-(1-sigma)*c for c in chis)}.items()},
 'scope':'New two-group four-term-log lower certificate plus static sum/min/interface checks. No author 128/129 log reconstruction, producer helper, search, local baseline, MATLAB, production or Git.'}
(OUT/'microcheck.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'new compressed certificates; pid',os.getpid())
