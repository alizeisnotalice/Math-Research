"""Independently replay r094 finite certificates and selected high-precision identities."""
import argparse,hashlib,json,re,itertools,sys,multiprocessing as mp
from pathlib import Path
from fractions import Fraction as Q
import mpmath as m

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(s):return m.mpf(re.sub(r'`[^*]*','',str(s)).replace('*^','e'))
def tilt_job(spec):
 si,ri=spec;out=[]
 for wp in (80,160):
  m.mp.dps=wp+70;w=[m.mpf(1)/8,m.mpf(3)/8,m.mpf(1),m.mpf(5)/2];r=[[Q(1,5),Q(7,10),Q(2,5),Q(1)],[Q(1,10),Q(3,10),Q(3,5),Q(9,10)],[Q(1)]*4,[Q(1,2),Q(1,3),Q(1,4),Q(1,5)]][ri-1];r=[m.mpf(t.numerator)/t.denominator for t in r];s=[m.mpf(1)/1000,m.mpf(1)/100,m.mpf(3)/25,m.mpf(1)/2,m.mpf(11)/10,m.mpf(3)][si-1];c=[-m.log(t) for t in r];p=[t/(1+t) for t in r];a=[v/(1+s*v) for v in w];d=[(1+s*v)*z for v,z in zip(w,c)];pt=[1/(1+m.exp(z)) for z in d];q=m.fsum(v+x*z*b for v,x,z,b in zip(a,w,c,pt));psi=m.fsum(m.log1p(m.exp(-z))-m.log1p(t)-m.log1p(s*v) for z,t,v in zip(d,r,w))
  def mix(sc,pr,tilt):
   coeff=[m.fprod(x/(x-y) for j,y in enumerate(sc) if j!=i) for i,x in enumerate(sc)];terms=[]
   for mask in itertools.product((0,1),repeat=4):
    t=q-m.fsum(b*x*z for b,x,z in zip(mask,w,c))
    if t<=0:continue
    prob=m.fprod(x if b else 1-x for b,x in zip(mask,pr));v=m.fsum(k*(m.exp(-s*t)*(-m.expm1(-(1/x-s)*t))/(1-s*x) if tilt else -m.expm1(-t/x)) for x,k in zip(sc,coeff));terms.append(prob*v)
   return m.fsum(terms)
  orig=mix(w,p,False);tilt=m.exp(psi+s*q)*mix(a,pt,True);assert abs(orig-tilt)<m.power(10,-wp);out.append([m.nstr(q,65),m.nstr(orig,65)])
 return spec,out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--project-root',required=True);ap.add_argument('--output',required=True);a=ap.parse_args();root=Path(a.project_root);data=root/'rounds/r094';sys.path.insert(0,str(root/'code'));from laplace94 import job
 load=lambda p:json.loads(p.read_text());wf=load(data/'wolfram/verification94.json');py=load(data/'python/curvature94.json');tf=load(data/'python/tilted94.json');nr=load(data/'python/normal94.json');res=load(data/'resources94.json')
 assert wf['pass_all'] and wf['symbolic_residuals']==['0']*4 and wf['kernels_closed'] and wf['full_curve_count']==4225 and wf['all_full_curves_pass']
 st=wf['steiner_rational'];mean=1-Q(1,320)**2/(2*Q(1,16));assert mean==Q(st['witness_mean'])==Q(12799,12800);assert Q(st['mean_margin'])==mean-Q(9999,10000)>0;assert Q(st['full_halo_difference_lower'])==Q(1,640)-Q(40000,99980001)>0
 for z in res.values():assert z['exit_code']==0 and z['stop_reason'] is None and z['reserve_GiB']==1 and z['min_system_available_MiB']>=1024
 assert py['pass_batch'] and py['cases']==4225 and py['actual_workers']==12 and len(py['selected_rows'])==67
 byj={r['j']:r for r in wf['curves']};pyj={r['dyadic_exponent_numerator']:r for r in py['selected_rows']};assert set(byj)==set(pyj)
 maxerr=m.mpf(0)
 with mp.get_context('spawn').Pool(12) as pool:
  for r in pool.imap_unordered(job,sorted(byj)):
   j=r['dyadic_exponent_numerator'];wr=byj[j];assert wr['pass']
   for z,w,old in zip(r['rows'],wr['rows'],pyj[j]['rows']):
    assert z['wp']==w['wp']==old['wp']
    with m.workdps(90):
     for key in ('c','curvature_ratio'):
      x=val(z[key]);e=abs(x-val(w[key]))/max(abs(x),m.mpf(1));assert e<m.mpf('1e-60'),(j,key,str(e));maxerr=max(maxerr,e)
     assert val(w['root_residual'])<m.power(10,-w['wp']) and val(w['stationary_value_error_over_u'])<m.power(10,-w['wp'])
  bytilt={(r['s_index'],r['r_index']):r for r in wf['tilts']};assert len(bytilt)==24
  for spec,rows in pool.imap_unordered(tilt_job,sorted(bytilt)):
   wr=bytilt[tuple(spec)];assert wr['pass']
   for x,w in zip(rows,wr['rows']):
    with m.workdps(90):
     e=max(abs(val(x[0])-val(w['q']))/max(1,abs(val(x[0]))),abs(val(x[1])-val(w['original_cdf'])));assert e<m.mpf('1e-60');maxerr=max(maxerr,e);assert val(w['absolute_error'])<m.power(10,-w['wp'])
 assert tf['pass_all'] and tf['total_cases']==315 and tf['actual_workers']==12 and tf['exact_identity_cases']==8;assert tf['max_two_tolerance_logratio_difference']<1e-10;assert tf['max_weighted_combined_tail_bound']<3e-14
 for row in tf['selected_rows']:
  assert abs(row['a']+(row['N']-1)*row['b']-row['N'])<1e-9
  for z in row['tol_rows']:assert 0<z['F_weighted']<1 and 0<z['F_reference']<1
 assert nr['pass_batch'] and nr['cases']==30848 and nr['actual_workers']==12 and nr['negative_observed_gain_count']==30848 and nr['positive_DKW_lower_bound_count']==0 and nr['unique_seeds']==30848
 sources=['check94.py','reduce94.py','verify94.wl','laplace94.py','normal_probe94.py','tilted_probe94.py','math_guard.py'];inputs=['python/curvature94.json','python/normal94.json','python/tilted94.json','wolfram/verification94.json','resources94.json'];docs=['core94.md','audit94.md','report.md','certification94.md','sources94.md'];prior=['rounds/r011/report.md','rounds/r021/report.md','rounds/r022/report.md','rounds/r028/report.md','rounds/r092/core92.md','rounds/r093/core93.md','rounds/r093/audit93.md']
 output=dict(pass_all=True,goal_status='active',general_upper_closed=False,cQ_class='all finite-volume unconditional coordinatewise downsets with convex logarithmic image',finite_onset_certified=False,exact_rational_Steiner_pass=True,symbolic_residual_count=4,selected_curvature_replayed=67,tilt_identities_replayed=24,max_independent_relative_or_absolute_error=m.nstr(maxerr,30),normal_cases=30848,positive_mixture_cases=315,source_sha256={f:sha(root/'code'/f) for f in sources},input_sha256={f:sha(data/f) for f in inputs},document_sha256={f:sha(data/f) for f in docs},prior_dependency_sha256={f:sha(root/f) for f in prior},scope='Finite rational and symbolic certificates, independent high-precision diagnostics and records; asymptotic full-weight and physical-class results rely on the mathematical proof. No general H_N or C_N upper closure, onset dimension or new global lower bound.')
 Path(a.output).write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
if __name__=='__main__':main()
