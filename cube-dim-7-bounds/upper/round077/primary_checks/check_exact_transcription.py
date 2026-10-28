from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')
src=P/'upper/round073/cantelli_frontier/verification.json'
v=json.loads(src.read_text())['key_rational_bounds']
s=F(3,2);qmax=F(11569,12800);sig=F(25,256);chi=F(63,1280)
eB=eN=eT=F(1,2**19);eJ=eRad=eS=eD=F(1,2**20)
Ret=F(9717,10000);ml=F(3043,10000);P0=1-eB-F(6,768)
ulo=P0*(1-eN)/(1+eJ);uhi=(1+eN)/(P0*(1-eJ));M=2*(1+eRad)*uhi
b=2*Ret*(1-eRad)*ulo-eS-eT;Gamma=ml*b/M
x0=F(231,256)*8192;m=F(2211,50000)/(1-eJ);p=F(183,1000);var=F(531,100000);gap=(chi-m)*x0-p
cap=((1-qmax)/(2*qmax)-chi)*(1-sig);tau=sig/2-chi*(1-sig)
checks={
 'original_tuple_not_equal_x_fraction':qmax!=1-sig,
 'computed_budget_matches_original_exact_receipt':Gamma==F(v['Gamma']),
 'deficit_margin':(1-eD)*Gamma>F(291,1000),
 'source_gap_and_retention':gap>0 and gap*gap/(var*x0+gap*gap)>Ret,
 'distinct_cap_and_tau_rates':cap==F(13627383,3790929920) and tau==F(1447,327680) and cap!=tau,
 'recorded_log_upper_and_price':F(v['price_upper'])<F(257,2500),
 'complete_root_price_slack':F(291,1000)-F(257,2500)-F(47,250)==F(1,5000),
 'old_Phi71_lower_bound_comparison':F(13,100)-F(9,100)+F(4,8192)<F(47,250),
 'Phi69_comparison_in_its_domain':F(291,1000)-sig+F(4,2**30)<F(23,100),
}
assert all(checks.values()),checks
def sh(path):return hashlib.sha256(path.read_bytes()).hexdigest()
receipt={'status':'PASS','pid':os.getpid(),'inputs':{'R073_original_exact_receipt_sha256':sh(src)},'script_sha256':sh(Path(__file__)),'checks':checks,'scope':'Independent finite Fraction recomputation of nine formula/parameter/comparison gates for manuscript integration; not a rerun of author program or proof of geometry/all-tail claims.'}
Path(__file__).with_name('check_exact_transcription.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,ensure_ascii=False))
