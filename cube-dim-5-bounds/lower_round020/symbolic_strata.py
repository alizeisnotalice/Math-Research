#!/usr/bin/env python3
"""Three exact symbolic capture-area arrangements, no radius or a1 grid."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import hashlib,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def save(p,d):p.write_text(json.dumps(d,separators=(',',':'))+'\n')
def rv(f,rho):return f[0]*rho+f[1]
def pv(p,rho):return p[0]*rho*rho+p[1]*rho+p[2]
def ps(p):return [str(x) for x in p]
def add(p,q):return tuple(x+y for x,y in zip(p,q))
def product_area(f,g):return (Q(f[0]*g[0],25),Q(f[0]*g[1]+f[1]*g[0],25),Q(f[1]*g[1],25))
def runtime(start):
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
class Limit(Exception):pass
def symbolic(m,atoms):
 start=time.monotonic();lo,hi=Q(m,10),Q(m+1,10);mid=Q(2*m+1,20);cells=0
 def check():
  if cells>10000 or time.monotonic()-start>30:raise Limit('30seconds_or10000cells')
 # Each tuple is (coefficient of rho, integer constant), in z=5*x coordinates.
 raws=[]
 for axis in ('z_x','z_y'):
  raws.append([(0,0),(0,25)]+[(sign*5,a[axis]) for a in atoms for sign in (-1,1)])
 full=set(raws[0])|set(raws[1]);pairs=0;roots_at_boundary=0;roots_outside=0;parallel=0
 for f,g in combinations(sorted(full),2):
  pairs+=1;check();delta=f[0]-g[0]
  if delta==0:parallel+=1;assert f[1]!=g[1];continue
  zero=Q(g[1]-f[1],delta)
  assert (zero*10).denominator==1 and not lo<zero<hi
  if zero in (lo,hi):roots_at_boundary+=1
  else:roots_outside+=1
 # The midpoint is used only for selecting the strict endpoint order and masks.
 bs=[]
 for a in atoms:
  b=[]
  for axis in ('z_x','z_y'):
   lower=max(((0,0),(-5,a[axis])),key=lambda f:rv(f,mid));upper=min(((0,25),(5,a[axis])),key=lambda f:rv(f,mid));b.extend((lower,upper))
  if rv(b[0],mid)<rv(b[1],mid) and rv(b[2],mid)<rv(b[3],mid):bs.append((a['id'],tuple(b)))
 axes=[];widths=[]
 for q in (0,2):
  funcs=sorted({(0,0),(0,25)}|{t for i,b in bs for t in b[q:q+2]},key=lambda f:rv(f,mid))
  assert len({rv(f,mid) for f in funcs})==len(funcs)
  w=[]
  for a,b in zip(funcs,funcs[1:]):
   delta=(b[0]-a[0],b[1]-a[1]);assert rv(delta,lo)>=0 and rv(delta,hi)>=0 and rv(delta,mid)>0;w.append(delta)
  assert rv(funcs[0],mid)==0 and rv(funcs[-1],mid)==25
  axes.append(funcs);widths.append(w)
 xm=[sum(1<<i for i,b in bs if rv(b[0],mid)<=rv(a,mid) and rv(z,mid)<=rv(b[1],mid)) for a,z in zip(axes[0],axes[0][1:])]
 ym=[sum(1<<i for i,b in bs if rv(b[2],mid)<=rv(a,mid) and rv(z,mid)<=rv(b[3],mid)) for a,z in zip(axes[1],axes[1][1:])]
 hist={};masks=[];counts={};trace=hashlib.sha256()
 for xi,dx in enumerate(widths[0]):
  row=[]
  for yi,dy in enumerate(widths[1]):
   cells+=1;check();mask=xm[xi]&ym[yi];count=mask.bit_count();area=product_area(dx,dy);hist[count]=add(hist.get(count,(Q(0),)*3),area);counts[count]=counts.get(count,0)+1
   row.append(format(mask,'x'));trace.update(f'{xi},{yi}:{mask:x}\n'.encode())
  masks.append(row)
 total=(Q(0),)*3
 for p in hist.values():total=add(total,p)
 assert total==(0,0,25)
 tails={}
 for n in range(1,17):
  p=(Q(0),)*3
  for count,v in hist.items():
   if count>=n:p=add(p,v)
  tails[str(n)]=ps(p)
 return {'status':'PASS_symbolic_open_stratum','m':m,'rho_open_stratum':[str(lo),str(hi)],'order_sample_rho':str(mid),'scaled_domain':[0,25,0,25],'coordinate_scaling':'z=5*x_normalized; normalized area=scaled area/25','affine_convention':'[rho_coefficient,integer_constant]','polynomial_convention':'[rho_squared,rho,constant] for normalized area','raw_endpoint_counts':[len(r) for r in raws],'unique_full_endpoint_counts':[len(set(r)) for r in raws],'endpoint_identities_removed':[len(r)-len(set(r)) for r in raws],'all_pair_linear_zero_guard':{'unique_union_functions':len(full),'nonidentical_pairs_checked':pairs,'parallel_distinct_pairs':parallel,'roots_at_stratum_boundary':roots_at_boundary,'roots_outside_closed_stratum':roots_outside,'proof':'all slopes are-5,0,5 and constants integers, so every finite zero is in(1/10)Z; none lies in the registered open stratum'},'endpoint_affines':axes,'width_affines':widths,'strict_width_guard':'each affine width is nonnegative at both closed stratum endpoints and positive at its midpoint, hence positive throughout the open stratum','active_window_count':len(bs),'cells':cells,'cell_count_histogram':{str(k):v for k,v in sorted(counts.items())},'cell_covering_atom_masks_hex':masks,'cell_mask_trace_sha256':trace.hexdigest(),'count_area_polynomials':{str(k):ps(v) for k,v in sorted(hist.items())},'all_counts_area_polynomial':ps(total),'threshold_area_polynomials':tails,'minimum_open_cell_count':min(hist),'maximum_open_cell_count':max(hist),'closed_capture_extrema_note':'inside this open stratum2rho is never in(1/5)Z; at every event line/corner one neighboring quadrant contains all touching atoms; period seam is harmless','runtime':runtime(start),'scope':'one full symbolic rho stratum only; no a1 union or parameter optimization'}
def main():
 start=time.monotonic();reg=load(HERE/'registration.json');assert sha(Path(__file__))==reg['source_sha256']
 for row in reg['inputs']:assert sha(PROJECT/row['path'])==row['sha256']
 print(json.dumps({'stage':'before_symbolic_arrangements','runtime':runtime(start),'strata':[14,15,16]}),flush=True)
 if os.environ.get('LOWER020_HANDOFF')=='1':
  deadline=time.monotonic()+20
  while not (HERE/'start.flag').exists() and time.monotonic()<deadline:time.sleep(.05)
 raw=load(PROJECT/'lower_round019/atoms.json');atoms=[]
 for a in raw:
  x=Q(5*a['x_grid'],100000);y=Q(5*a['y_grid'],100000);assert x.denominator==y.denominator==1
  assert x==4*a['m']-3*a['n'] and y==3*a['m']+4*a['n'];atoms.append({'id':a['id'],'m':a['m'],'n':a['n'],'z_x':int(x),'z_y':int(y)})
 assert len(atoms)==len({(a['m'],a['n']) for a in atoms})==225
 save(HERE/'scaled_atoms.json',atoms);done=[];out={}
 for m in (14,15,16):
  q=HERE/f'stratum_m{m}.json'
  if q.exists():r=load(q);assert r['m']==m and r['status']=='PASS_symbolic_open_stratum'
  else:
   try:r=symbolic(m,atoms)
   except Limit as e:save(HERE/'checkpoint.json',{'status':'CHECKPOINT_LIMIT','next_stratum':m,'completed':done,'reason':str(e),'runtime':runtime(start)});return
   r['scaled_atom_manifest_sha256']=sha(HERE/'scaled_atoms.json');r['source_sha256']=sha(Path(__file__));r['registration_sha256']=sha(HERE/'registration.json');save(q,r)
  done.append({'m':m,'path':q.name,'sha256':sha(q)});out[m]=r
  save(HERE/'checkpoint.json',{'status':'saved_after_one_symbolic_stratum','completed':done,'runtime':runtime(start)})
  print(json.dumps({'m':m,'cells':r['cells'],'histogram_polynomials':r['count_area_polynomials'],'sum':r['all_counts_area_polynomial'],'seconds':r['runtime']['wall_seconds']}),flush=True)
 # Evaluate saved019 endpoints only, never run another arrangement.
 prior=load(PROJECT/'lower_round019/results.json');bundle=load(PROJECT/'lower_round019/arrangement_certificates.json');checks=[]
 for row in bundle['records']:
  r=row['record'];rho=Q(r['rho']);m=(10*rho).numerator//(10*rho).denominator;assert m in out
  assert Q(m,10)<rho<Q(m+1,10);h=out[m]['count_area_polynomials']
  assert set(h)==set(r['normalized_area_by_count'])
  for count,p in h.items():assert pv(tuple(map(Q,p)),rho)==Q(r['normalized_area_by_count'][count])
  assert out[m]['minimum_open_cell_count']==r['minimum_closed_period_capture'] and out[m]['maximum_open_cell_count']==r['maximum_closed_period_capture']
  checks.append({'N':r['N'],'mode':r['mode'],'rho':str(rho),'m':m,'all_histogram_polynomials_match_saved019_exactly':True,'original_saved_record_sha256':row['original_sha256']})
 cases=[]
 for c in prior['cases']:
  N=c['N'];lo,hi=map(Q,c['rho_pair']);m=(10*lo).numerator//(10*lo).denominator;assert Q(m,10)<lo<hi<Q(m+1,10)
  p=out[m]['threshold_area_polynomials'][str(c['n'])];coef=tuple(map(Q,p));assert [str(pv(coef,rho)) for rho in (lo,hi)]==c['normalized_capture_area_interval']
  cases.append({'K':c['K'],'n':c['n'],'N':N,'m':m,'normalized_exact_area_polynomial':p,'rho_exact':'w/(1+u), w^5=N,u^5=2','physical_exact_area_formula':f'(({p[0]})*w^2+({p[1]})*w*(1+u)+({p[2]})*(1+u)^2)/4; w^5={N},u^5=2','matches019_inner_outer_area_evaluations':True,'was019_partial':[c['K'],c['n']] in prior['remaining_partial_cases'],'prior_status':c['status'],'status':'EXACT_PLANAR_AREA_FOR_THIS_CASE'})
 solved=[c for c in cases if c['was019_partial']];assert len(solved)==7 and len(checks)==28
 result={'status':'COMPLETE_three_symbolic_strata_and_seven_partial_exact_areas','strata':done,'case_count':len(cases),'previous_partial_cases_resolved':[[c['K'],c['n']] for c in solved],'remaining_partial_exact_area_cases':[],'all_case_exact_formulas':cases,'saved019_endpoint_evaluations':checks,'all_histogram_sums_identically25':True,'all_widths_strict_inside_strata':True,'all_nonidentical_endpoint_zeros_outside_open_strata':True,'atom_completeness':'all normalized rho<1.7<5;25 period representatives plus h,k=-1,0,1 include every relevant distinct atom','distinct_atoms_and_common_radius':'Cartesian product of n distinct planar indices with K distinct unrotated indices givesKn distinct atoms; common r=(K*n)^(1/5)/2. This result evaluates planar areas only, not five-dimensional witness unions.','old_exclusions_not_claimed_as_first':'019/earlier empty-case conclusions are only crosschecked; new content is whole-stratum exact topology and symbolic areas','source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),'total_symbolic_cells':sum(out[m]['cells'] for m in out),'runtime':runtime(start),'scope':'exact planar capture areas for current16 cases and three open rho strata; no sum of case areas as union, a1 optimization, fullscale optimum orC5 breakthrough'}
 save(HERE/'results.json',result);save(HERE/'checkpoint.json',{'status':'COMPLETE_three_registered_strata_only','completed':done,'result_sha256':sha(HERE/'results.json'),'runtime':runtime(start),'no_further_task_started':True});print(json.dumps({k:v for k,v in result.items() if k not in ('all_case_exact_formulas','saved019_endpoint_evaluations')},indent=2),flush=True)
if __name__=='__main__':
 assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
 main()
