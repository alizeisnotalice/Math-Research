"""Read-only: exact saved endpoints/weights/hashes; no anchor or native replay."""
import sys
sys.dont_write_bytecode=True
import json,hashlib,math
from pathlib import Path
from fractions import Fraction as Q
from run_full_once import validate,bounds
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 m=json.loads((P/'manifest.json').read_text())
 for f,h in m['whitelist_sha256'].items():assert sha(P/f)==h,('local',f)
 for f,h in m['external_readonly_sha256'].items():assert sha(Path(f))==h,('external',f)
 v=json.loads((P/'scalar_validation.json').read_text());n=json.loads((P/'scalar_native.json').read_text())
 assert v['status']=='PASS' and v['scalar_Hessian_checks']==len(v['rows'])==len(n['scalar_Hessian_checks'])==16
 assert v['native_sha256']==sha(P/'scalar_native.json')
 for ref,native in zip(v['rows'],n['scalar_Hessian_checks']):
  assert (ref['set'],ref['i'],ref['j'])==(native['set'],native['i'],native['j'])
  assert ref['diagonal_multiplier']==(2 if ref['i']==ref['j'] else 1)
  w=list(map(Q,ref['polynomial_leading_weights_exact']));conds=[(Q(x),k) for x,k in ref['conditions']]
  for d in range(8):assert sum(q*(0 if d<k else math.comb(d,k)*x**(d-k)) for q,(x,k) in zip(w,conds))==Q(d==7)
  lo,hi=bounds(native['hessian']);a,b=bounds(ref['reference_dyadic']);assert lo<=a<=b<=hi
 rec=json.loads((P/'full_execution.json').read_text());assert rec['attempt']==1 and rec['exit_code']==0 and not rec['stderr']
 b=json.loads((P/'build_record.json').read_text());assert rec['source_sha256']==b['source_sha256']==sha(P/'band_hessian.cpp') and rec['binary_sha256']==b['binary_sha256']
 old=json.loads((P.parents[1]/'round016/band_engine/native_result.json').read_text());full=json.loads((P/'native_hessian.json').read_text())
 check=validate(full,old);assert check==json.loads((P/'full_validation.json').read_text())
 assert not check['complete_Hessian_independent_reference'] and len(full['hessian_upper_lex'])==21
 print('PASS saved-only: hashes,16 interpolation receipts/weights,old9-field cross,21 Hessian endpoints,existing Euler;0 scalar/0 native replays/0 writes')
if __name__=='__main__':main()
