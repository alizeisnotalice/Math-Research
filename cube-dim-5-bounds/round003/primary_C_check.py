from pathlib import Path
from fractions import Fraction as F
import importlib.util,time,json,hashlib
root=Path(__file__).resolve().parents[1]
audit=root/'audit/012_portable_C_verify.py'
spec=importlib.util.spec_from_file_location('independent_C',audit);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
start=time.monotonic();boxes=m.retained_C_boxes();volume,segments,nodes=m.bitmask_union(tuple(b for _,b in boxes),start+120);ratio=32*volume;lo,hi=ratio.interval()
result_path=root/'lower_round003/a2_equals_a3_results.json'
lock_path=root/'lower_round003/lock_c_a2_equals_a3.json'
b_path=root/'lower_round003/lock_b_1_32.json'
prod=json.loads(result_path.read_text());lock=json.loads(lock_path.read_text());b=json.loads(b_path.read_text())
assert [str(x) for x in ratio.coeff]==[prod['coefficients'].get(f'{i//5},{i%5}','0') for i in range(25)]
assert F(lock['exact_interval'][0])<=lo<=hi<=F(lock['exact_interval'][1])
old_source=json.loads((root/'lower/rectangular_certificate.json').read_text())
old=max(old_source['candidates'],key=lambda c:F(c['ratio_interval'][1]))['ratio_interval']
assert lo>F(b['exact_interval'][1]) and lo>F(old[1]) and lo>F(19943,10000)
out={'status':'PASS_primary_actual_independent_C_reconstruction_and_exact_comparison','retained':len(boxes),'nodes':nodes,'ratio_interval':[str(lo),str(hi)],'new_lower_minus_B_upper':str(lo-F(b['exact_interval'][1])),'new_lower_minus_locked001_upper':str(lo-F(old[1])),'new_lower_minus_19943_over_10000':str(lo-F(19943,10000)),'seconds':time.monotonic()-start,'audit_source_sha256':hashlib.sha256(audit.read_bytes()).hexdigest(),'C_lock_sha256':hashlib.sha256(lock_path.read_bytes()).hexdigest(),'C_result_sha256':hashlib.sha256(result_path.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
raw=json.dumps(out,indent=2)+'\n';Path(__file__).with_suffix('.json').write_text(raw);print('PASS: 42 boxes; exact C > B,001,19943/10000; seconds',out['seconds'])
