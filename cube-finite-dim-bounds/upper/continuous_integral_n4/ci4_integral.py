"""Continuous per-region integral range diagnosis, not a target certificate.
Shared sign-DD evaluator is upper_theory's implementation, not duplicated here.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util,hashlib,json,argparse,time
from flint import arb,ctx
import ci4_partition as partition
ctx.prec=192
DD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-upper/cube-finite-dim-bounds/upper/u11_minimum_dd.py')
DD_SHA='2ed3070118111b93b7f5091e8d66145d42a49ac6297ba4f4974e148cf6743797'
def a(q):
 q=Q(q);return arb(q.numerator)/q.denominator

def interval(x):return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display':str(x)}
def p43(x):return arb(0) if x==0 else x.root(3)**4

def run(save):
 start=time.time();assert hashlib.sha256(DD.read_bytes()).hexdigest()==DD_SHA
 assert hashlib.sha256(DD.with_name('u11_green.py').read_bytes()).hexdigest()==partition.SOURCE_SHA
 sp=importlib.util.spec_from_file_location('shared_dd',DD);dd=importlib.util.module_from_spec(sp);sp.loader.exec_module(dd)
 g,ps,hp,counts=partition.build();pi=arb.pi();D=1/(48*pi**4);B=arb(5632)*32/pi**3
 finite_u=arb(0);finite_l=arb(0);fallback=arb(0);count_u=0;count_l=0;volume_u=Q(0);volume_l=Q(0);rows=[];upper_failed=[]
 for i,p in enumerate(ps):
  vol=partition.volume(p);vertices=[(Q(1),)+v for v,active in p]
  nodes=dd.normalized_vertex_nodes(vertices);out=dd.green_coefficient_interval(nodes)
  # Independent universal radial upper range, valid on the entire polytope.
  # sum corner deficits <= sqrt(3)*r, so r >= min(sum deficits)/sqrt(3).
  deficit=min(Q(3)-sum(v[1:]) for v in vertices);assert deficit>=partition.EPS
  radial=D+B*arb(3).sqrt()/a(deficit)
  if out['lower_pi3'] is not None:
   clo=max(Q(0),out['lower_pi3']);finite_l+=a(vol)*p43(a(clo)/pi**3);count_l+=1;volume_l+=vol
  if out['upper_pi3'] is not None:
   chi=out['upper_pi3'];assert chi>=0
   finite_u+=a(vol)*p43(a(chi)/pi**3);count_u+=1;volume_u+=vol
  else:
   upper_failed.append(i)
  # Every region has a rigorous radial fallback; this is usually uselessly large.
  fallback+=a(vol)*p43(radial)
  if len(rows)<12:rows.append({'region':i,'volume':str(vol),'lower_pi3':None if out['lower_pi3'] is None else str(out['lower_pi3']),'upper_pi3':None if out['upper_pi3'] is None else str(out['upper_pi3']),'upper_poles':out['upper_pole_count'],'lower_poles':out['lower_pole_count'],'first_upper_pole':out['first_upper_pole']})
 eps=a(partition.EPS)
 def E(r):return arb(2).root(3)*4*pi*(p43(D)*r**3/3+arb(3)/5*p43(B)*r.root(3)**5)
 cap=E(eps)/8
 out={'status':'target not certified; continuous range diagnostic','precision_bits':ctx.prec,'shared_evaluator':str(DD),'shared_evaluator_sha256':DD_SHA,'shared_source_sha256':partition.SOURCE_SHA,'regions':len(ps),'finite_upper_regions':count_u,'finite_lower_regions':count_l,'finite_upper_ordered_volume':str(volume_u),'finite_lower_ordered_volume':str(volume_l),'finite_upper_region_I4_upper_expression':interval(6*finite_u),'certified_I4_lower_expression':interval(6*finite_l),'radial_fallback_all_regions_I4_upper_expression':interval(6*fallback+cap),'crude_global_I4_upper_expression':interval(E(arb(3).sqrt())/8),'omitted_corner_I4_budget':interval(cap),'unresolved_upper_regions':upper_failed,'first_region_diagnostics':rows,'seconds':time.time()-start}
 if save:Path(save).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ['first_region_diagnostics','unresolved_upper_regions']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--save');run(p.parse_args().save)
