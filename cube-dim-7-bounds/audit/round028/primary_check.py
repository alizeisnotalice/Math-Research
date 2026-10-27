from pathlib import Path
from fractions import Fraction as Q
import importlib.util,json,time
p=Path(__file__).resolve().parents[2];spec=importlib.util.spec_from_file_location('cap_lib',p/'audit/round027/cap_lib.py');lib=importlib.util.module_from_spec(spec);spec.loader.exec_module(lib)
r=json.loads((p/'upper/round028/cap_candidate/result.json').read_text());start=time.perf_counter();raw=float.fromhex(r['raw_t_hex']);q=round(raw*2**20);assert q==r['integer_q'] and r['denominator']==2**20;t=Q(q,2**20);assert str(t)==r['exact_t'];B,F,c,k,terms=lib.cap(t);U=Q(r['candidate_grid_upper']);assert F[1]<U and lib.base.grid(F[1])==U
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
assert U<dy(r['R027_winner_interval']['lower'])
print(json.dumps(dict(status='PASS exact Fraction final-candidate certificate',raw_t_hex=r['raw_t_hex'],exact_t=str(t),Bcap=lib.base.encode(B),F=lib.base.encode(F),rational_constant=str(c),log_coefficient=str(k),strict_C7_upper=str(U),old_complete_lower=str(dy(r['R027_winner_interval']['lower'])),local_search_global_optimality_claim=False,wall_seconds=time.perf_counter()-start),indent=2))
