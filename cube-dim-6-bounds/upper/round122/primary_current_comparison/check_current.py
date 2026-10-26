import sys
sys.set_int_max_str_digits(0)
import json,hashlib,datetime,os
from fractions import Fraction as F
from pathlib import Path
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
r=json.loads((P/'registration.json').read_text());data={}
for z in r['sources']:
 assert sha(z['path'])==z['sha256'];data[z['round']]=json.loads(Path(z['path']).read_text())
with (P/'started.json').open('x')as f:json.dump({'PID':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':sha(__file__),'registration_sha256':sha(P/'registration.json')},f)
a=data[122]['configurations'][3];newL,newU=map(F,a['new_far_sufficient_allocation']);oldL,oldU=map(F,a['old_far_allocation']);delta=F(a['I_recovery_exact']);assert newL==oldL+delta and newU==oldU+delta and newU>=newL>0
w116=F(data[116]['configurations'][0]['I_envelope_lower_to_fifth']);c120=data[120]['configurations'][1];coarse=F(c120['J_far_lower_fifth_exact']);fine=F(c120['J_far_refined_lower_fifth'])
checks={name:{'lower_fifth_exact':str(w),'residual_to_current_upper_fifth':str(w-newU**5),'strictly_excludes_specified_expression_current_allowance':w>newU**5}for name,w in [('116_window',w116),('120_coarse',coarse),('120_sameR_refined',fine)]}
out={'status':'COMPLETED_ONE_NEW_CURRENT_BUDGET_COMPARISON','new_far_lower':str(newL),'new_far_upper':str(newU),'delta_exact':str(delta),'comparison':checks,'actual_I_or_model_excluded':False,'new_C6_bound':False,'source_sha256':sha(__file__),'registration_sha256':sha(P/'registration.json')}
with (P/'results.json').open('x')as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps({k:v['strictly_excludes_specified_expression_current_allowance']for k,v in checks.items()}))
