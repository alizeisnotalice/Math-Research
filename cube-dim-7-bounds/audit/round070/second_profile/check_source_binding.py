"""New R070 frozen-source scope/finite-guard extension check.

Reads the earlier independent H2 receipt, without rerunning its symbolic
checks or log sums. This separate run binds the source available afterward.
"""
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,sys

HERE=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/upper/round070/second_profile')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
expected={'PROOF.md':'4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8',
 'verify.py':'0f28087a2600ac84d7dfd0f1cf8b35439097be4251e0e8ab8e7ff2e5c08857ee',
 'verification.json':'8d1e1d722b2a773bda7a42db4be259ea1c0d8b8dbf1b8448c5016273742be97c'}
checks={}
def ok(k,p):
 assert bool(p),k
 checks[k]=True
for k,v in expected.items():ok('source_identity_'+k,sha(SRC/k)==v)
priorpath=HERE/'verification.json'
prior=json.loads(priorpath.read_text())
ok('independent_prior_run_passed',prior['status']=='PASS' and prior['check_count']==29)
free_lo=F(prior['values']['free_mean_lower']['rational'])
free_hi=F(prior['values']['free_mean_upper']['rational'])
hc=F(1,2**20); qlo=F(2,3); nlo=2000
full_lo=(1-F(1,nlo))*free_lo
full_hi=free_hi+F(1,250*qlo*nlo)
ok('actual_n2000_full_target_mean_lower',full_lo>F(528,100000))
ok('actual_n2000_full_target_mean_upper',full_hi<F(5295,1000000))
ok('two_sided_event_distance_at_n2000',full_lo-F(527,100000)>F(1,200000) and F(53,10000)-full_hi>F(1,200000))
ok('n2000_target_Chebyshev_constant',F(1,125)**2/(qlo**2*F(1,200000)**2)==5760000)
ok('strong_matched_interval_lower',F(527,100000)-hc/5>F(5269,10**6))
ok('strong_matched_interval_upper',F(53,10000)+hc/5<F(5301,10**6))
ok('first_profile_wide_h_transport',hc<F(1,100000))

eta=F(1,1024); c1=F(3,50); c2=F(1,125)
l2,u2=F(5269,10**6),F(5301,10**6)
varlo=(l2-c2/(nlo*qlo))/(1+eta)-c1**2*F(1,1000)/(qlo**2*(1-eta)**2)
varhi=u2/(1-eta)
ok('J_normalization_from_x_over_n',F(1,1000)/(3*qlo)<eta)
ok('centered_source_variance_lower',varlo>F(13,2500))
ok('centered_source_variance_upper',varhi<F(531,100000))
ok('source_variance_lower_fraction',varlo==F(9776676208,1862312890625))
ok('source_variance_upper_fraction',varhi==F(912,171875))

L=4096; e=F(1,2**20); delta=F(1,2**19)
hwide=F(1,256*L); hnarrow=F(1,512*L)
ok('wide_matching_h_guard',hwide==hc==delta/2 and 0<hnarrow<hwide)
ok('local_n_and_density_from_Taylor2',F(L**2,2)>=2000 and F(2,L)<F(1,1000) and F(L,3)>6)
ok('local_relative_guard',F(8,L**5)<2*hnarrow/3 and L**4>6144)
# Deliberately weaker independent lower Taylor surrogates than the author's:
# exp L >= L^3/8 and >= L^4/32, since 3!<8 and 4!<32.
ok('weaker_Taylor_lower_bounds',6<8 and 24<32)
ok('uniform_J_with_Taylor3_over8',F(4,L**2)<e)
bad=F(5760000*32,L**4)
ok('new_target_failure_with_Taylor4_over32',bad<e)
out={'status':'PASS LOCAL frozen-source extension','utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),'python':sys.executable,
 'source_hashes':{str(SRC/k):sha(SRC/k) for k in expected},'code_sha256':sha(Path(__file__)),
 'prior_independent_receipt_sha256':sha(priorpath),'check_count':len(checks),'checks':checks,
 'values':{k:{'rational':str(v),'display_only':float(v)} for k,v in {'full_target_mean_lower_n2000':full_lo,'full_target_mean_upper_n2000':full_hi,'variance_lower':varlo,'variance_upper':varhi,'profile_bad_Taylor4_over32':bad}.items()},
 'scope':'Read/bind final local source after independent formula run; n2000/variance/wide-h finite guards only. No author code or old independent symbolic/log checks rerun; no Cn conclusion.'}
(HERE/'source_binding_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'checks':len(checks),'values':{k:v['display_only'] for k,v in out['values'].items()}},indent=2))
