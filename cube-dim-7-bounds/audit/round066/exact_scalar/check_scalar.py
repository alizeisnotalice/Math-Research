"""One new seven-term log-price enclosure; not an operator certificate."""
from fractions import Fraction as F
from pathlib import Path
import datetime,hashlib,json,os,sys
root=Path(__file__).resolve().parent
source=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/PROOF.md')
expected='4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1'
hashof=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert hashof(source)==expected
u=F(83,917); gamma=F(317,1000); sigma=F(83,1000)
S6=sum(((-1)**(k+1)*u**k/k for k in range(1,7)),F())
S7=S6+u**7/7
price_upper=F(4333,50000); c=F(23,100)
checks={
 'log1p_domain':0<u<1,
 'exact_price_argument':1+u==1/(1-sigma)==F(1000,917),
 'seven_term_upper_exceeds_six_term_lower':S7>S6,
 'S7_below_chosen_rational_upper':S7<price_upper,
 'net_slack_identity':gamma-price_upper-c==F(17,50000)>0,
 'sigma_positive_decay':0<sigma<1,
 'exact_bound_coefficient_lower':gamma-price_upper==F(11517,50000)>c,
 'display_coefficient_improves_old_display':c>F(9,40),
}
assert all(checks.values()),[name for name,ok in checks.items() if not ok]
receipt={
 'status':'PASS: conditional scalar corollary for unchanged Phi66; full operator ACCEPT required',
 'actual_pid':os.getpid(),'actual_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,
 'method':'seven-term odd alternating log1p upper, six-term even lower; integral remainder signs proved in SCALAR_PROOF.md',
 'parameters':{'gamma':str(gamma),'sigma':str(sigma),'u':str(u),'price_upper':str(price_upper),'chosen_c':str(c)},
 'log_lower_S6':str(S6),'log_upper_S7':str(S7),'upper_gap':str(price_upper-S7),'net_slack':str(gamma-price_upper-c),
 'checks':checks,'check_count':len(checks),
 'source':{'path':str(source),'sha256':hashof(source)},
 'own_sources':{name:hashof(root/name) for name in ('REGISTRATION.md','SCALAR_PROOF.md','check_scalar.py')},
 'scope':'same complete F(A,B-gamma*n); improved scalar display only; no new kernel/gamma/Phi; no explicit N; general L1/all radii deferred to full review',
 'no_author_verifier_replay':True,'no_old_grid':True,'no_production':True,'no_MATLAB':True,'no_Git':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),'actual_pid':os.getpid(),'S7':str(S7),'price_upper_gap':str(price_upper-S7),'net_slack':'17/50000','proof_sha256':receipt['own_sources']['SCALAR_PROOF.md']},indent=2))
