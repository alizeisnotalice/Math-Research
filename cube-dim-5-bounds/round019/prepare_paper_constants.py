from fractions import Fraction as Q
from pathlib import Path
import json,hashlib,os,time
p=Path(__file__).parent;r=json.loads((p/'paper_constants_registration.json').read_bytes())
assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==r['code_sha256']
s=Path(r['source']['path']);assert hashlib.sha256(s.read_bytes()).hexdigest()==r['source']['sha256']
a=json.loads(s.read_bytes());t=time.monotonic();g=Q(a['qualified_g']);I=Q(a['complete_I_upper']);N=10**20
G=(g.numerator*N)//g.denominator;H=(I.numerator*N+I.denominator-1)//I.denominator
L,U=Q(G,N),Q(H,N);C=Q(a['simple_rational_upper']);gap=C**4*L**5-3125*U**4
assert 0<L<=g and I<=U and gap>0
x={'status':'exact_outward_paper_constants_same_certified_C5_upper','pid':os.getpid(),'seconds':time.monotonic()-t,'denominator':str(N),'g_lower_numerator':str(G),'I_upper_numerator':str(H),'g_lower':str(L),'I_upper':str(U),'C5_upper':str(C),'strict_polynomial_gap':str(gap),'from_original_exact_certificate_not_display_decimals':True,'source':r['source'],'new_H_mean_log_R_root_calls':0,'paper_or_second_audit_complete_claim':False}
with (p/'paper_constants.json').open('x') as f:json.dump(x,f,indent=2);f.write('\n')
print(json.dumps({k:x[k] for k in ['status','pid','g_lower_numerator','I_upper_numerator','denominator','C5_upper']}))
