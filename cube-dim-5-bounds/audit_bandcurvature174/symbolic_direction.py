#!/usr/bin/env python3
"""One new symbolic positive-numerator denominator-direction diagnostic only."""
from pathlib import Path
import argparse,hashlib,json,os,time
import sympy as s
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--registration',type=Path,required=True);ap.add_argument('--sha',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    assert sha(a.registration)==a.sha and not a.output.exists();r=json.loads(a.registration.read_bytes())
    assert r['new_symbolic_case']=='positive_Kupper_over_Hlower_residual_only'
    for p,h in {**r['source_locks'],**r['code_locks']}.items():assert sha(p)==h
    t=time.monotonic();h,l=s.symbols('h l',positive=True);k,d=s.symbols('k d',nonnegative=True)
    # KU=K+d, h>=l; positive residual avoids the153 negative-tail ratio rule.
    residual=s.factor((k+d)/(16*l)-k/(16*h));expected=(d*h+k*(h-l))/(16*h*l)
    assert s.cancel(residual-expected)==0
    result={'status':'PASS_new_symbolic_denominator_direction_identity','pid':os.getpid(),'seconds':time.monotonic()-t,'residual':str(residual),'identity_difference':'0',
      'nonnegative_reason':'d>=0,k>=0,h>=l>0','registration_sha256':a.sha,'real_H_fields_read':0,'real_q_bandproducts_corners':0,'old153_weight_or_symbolic_baseline_replays':0}
    with a.output.open('x') as f:json.dump(result,f,sort_keys=True,indent=2);f.write('\n')
    print(json.dumps(result))
if __name__=='__main__':main()
