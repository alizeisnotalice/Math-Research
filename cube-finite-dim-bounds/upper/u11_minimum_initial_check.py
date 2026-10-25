"""One initial <=2048-cell G lower diagnostic. No I/cost/tau search."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util,hashlib,json,time
from flint import arb,ctx
import u11_minimum_dd as dd
ctx.prec=128
GEOM=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-integral4/cube-finite-dim-bounds/upper/continuous_integral_n4/ci4_partition.py')
sp=importlib.util.spec_from_file_location('u11_minimum_shared_geometry',GEOM);geom=importlib.util.module_from_spec(sp);sp.loader.exec_module(geom)
start=time.time();g,polys,hp,counts=geom.build();assert len(polys)<=2048
TARGET=Q(21,25000)
def aq(q):q=Q(q);return arb(q.numerator)/q.denominator
# Explicit lower cap: one-large-jump event in a Hadamard body basis.
# Remaining four variables |C_t|<=t/2 have prob >=(4/(5*pi))^4.
# For |a|<=1, |z|<=epsilon, p_t>=256/(12500000*pi^8)t^-2 on epsilon<=t<=1.
eps=geom.EPS
cap_lower=Q(256,12500000)/Q(22,7)**8*(1/eps-1)
assert cap_lower>TARGET
stats={'certified':0,'finite_lower_but_not_certified':0,'individual_lower_pole':0};vols={k:Q(0) for k in stats};rows=[]
for ix,p in enumerate(polys):
 vertices=[(Q(1),)+v for v,a in p]
 nodes=dd.normalized_vertex_nodes(vertices)
 rr=dd.green_coefficient_interval(nodes)
 vol=geom.volume(p)
 if rr['lower_pi3'] is None:status='individual_lower_pole'
 elif aq(rr['lower_pi3'])/(arb.pi()**3)>=aq(TARGET):status='certified'
 else:status='finite_lower_but_not_certified'
 stats[status]+=1;vols[status]+=vol
 rows.append({'id':ix,'status':status,'vertices':[[str(z) for z in x] for x in vertices],
              'volume':str(vol),'lower_pi3':str(rr['lower_pi3']) if rr['lower_pi3'] is not None else None,
              'lower_pole_count':rr['lower_pole_count'],'upper_pole_count':rr['upper_pole_count'],
              'first_lower_pole':rr['first_lower_pole']})
 if ix%50==0:print('initial cell',ix,'/',len(polys),stats,flush=True)
assert sum(vols.values())==Q(1,6)-eps**3/36
out={'batch':'11K','stage':'initial affine mixed-sign DD, no representation revision yet','target':str(TARGET),
     'Arb_precision_bits':ctx.prec,'initial_cell_count':len(polys),'max_budget':2048,'counts':stats,
     'volumes':{k:str(v) for k,v in vols.items()},'cap_epsilon':str(eps),'cap_lower_exact':str(cap_lower),
     'cap_lower_display':float(cap_lower),'cap_g_certified':True,'global_g_certified':stats['certified']==len(polys),
     'shared_geometry_path':str(GEOM),'shared_geometry_sha256':hashlib.sha256(GEOM.read_bytes()).hexdigest(),
     'stable_DD_sha256':hashlib.sha256(Path(dd.__file__).read_bytes()).hexdigest(),
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cells':rows,'elapsed_seconds':time.time()-start}
Path(__file__).with_name('u11_minimum_initial_results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print('INITIAL DONE',stats,'cap',float(cap_lower),'seconds',out['elapsed_seconds'],flush=True)
