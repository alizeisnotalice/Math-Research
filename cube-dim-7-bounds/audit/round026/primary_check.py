from pathlib import Path
from fractions import Fraction as Q
import importlib.util,json,time
p=Path(__file__).resolve().parents[2];spec=importlib.util.spec_from_file_location('rational_lib',p/'audit/round024/rational_lib.py');lib=importlib.util.module_from_spec(spec);spec.loader.exec_module(lib)
r=json.loads((p/'upper/round026/buffer_parameter_numeric/result.json').read_text());start=time.perf_counter();lo=Q(r['root_bracket_lower']);hi=Q(r['root_bracket_upper']);t=Q(r['candidate']);assert hi-lo==Q(1,2**80) and t==(lo+hi)/2
l=lib.G26(lo);u=lib.G26(hi);assert l[1]<0<u[0];f=lib.F25(t);U=Q(r['candidate_grid_upper']);assert f[1]<U and lib.grid(f[1])==U
def dy(x):return Q(int(x[0]))*Q(2)**int(x[1])
old=json.loads((p/'upper/round025/buffer_scalar/result.json').read_text())['winner'];assert r['reused_R025_row']==old and U<dy(old['interval']['lower'])
print(json.dumps(dict(status='PASS independent exact Fraction root endpoints and candidate scalar',source_candidate=str(t),root_bracket=[str(lo),str(hi)],G_endpoints=[lib.encode(l),lib.encode(u)],F_candidate=lib.encode(f),strict_C7_upper=str(U),old_complete_lower=str(dy(old['interval']['lower'])),scope='candidate interval not exact minimum-value enclosure; root location only this upper-expression family; no old F re-evaluation',wall_seconds=time.perf_counter()-start),indent=2))
