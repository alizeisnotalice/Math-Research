"""Replay physical periodic witnesses and directed Lorenz-moment certificates."""
import argparse,hashlib,itertools,json,re,sys,math
from pathlib import Path
from fractions import Fraction as Q
import mpmath as m

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def num(s):return m.mpf(re.sub(r'`[^*]*','',str(s)).replace('*^','e'))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--project-root',required=True);ap.add_argument('--output',required=True);args=ap.parse_args();root=Path(args.project_root);data=root/'rounds/r096';sys.path.insert(0,str(root/'code'));from information96 import job;from sliver79 import logi
 load=lambda f:json.loads((data/f).read_text());info=load('information/information96.json');wf=load('information/information_wolfram96.json');grid=load('python/periodic96.json');resources=load('resources96.json')
 assert info['pass_all'] and info['cases']==56 and wf['pass_all'] and wf['symbolic_residuals']==['0']*4 and wf['kernels_closed'];by_r={x['r']:x for x in wf['rows']};assert len(by_r)==7
 best=info['best_certificate'];base=info['base_certificate'];assert Q(best['information_coefficient_lower'])>Q('0.12324');assert Q(base['information_coefficient_lower'])>Q('0.12')
 for rec in (base,best):
  replay=job((rec['r'],rec['mesh'],rec['mode']));assert replay==rec
 for rec in info['rows']:
  assert rec['all_root_upper_checks'];r=Q(rec['r']);assert Q(rec['information_coefficient_lower'])==((r-1)*Q(rec['mu_lower'])-Q(rec['log_J_upper']))/r
  with m.workdps(90):
   v=by_r[rec['r']]['rows'];assert abs(num(v[0]['J'])-num(v[1]['J']))<m.mpf('1e-60');ju=Q(rec['J_upper']);assert m.mpf(ju.numerator)/ju.denominator>num(v[1]['J'])
 for z in resources.values():assert z['exit_code']==0 and z['stop_reason'] is None and z['reserve_GiB']==1 and z['min_system_available_MiB']>=1024
 assert grid['pass_all_candidates'] and grid['actual_workers']==12 and grid['case_count']==12
 physical=[];rect_count=0;minmargin=None
 for case in grid['cases']:
  dim,mx=case['dimension'],case['cells_per_axis'];p=Q(case['p']);theta=Q(case['theta']);br=list(map(Q,case['axis_breaks']));assert 0<p<1 and 0<theta<=1 and dim>=1 and mx>=1;assert br[0]==0 and br[-1]==1 and len(br)==mx+1 and all(a<b for a,b in zip(br,br[1:]))
  coords=list(itertools.product(range(mx),repeat=dim));index={x:i for i,x in enumerate(coords)};vol=[math.prod(br[j+1]-br[j] for j in x) for x in coords];cand=case['candidate'];den=cand['source_rounding_denominator'];xs=cand['source_cell_fraction_numerators'];assert cand['certificate_pass'] and den==2**24 and len(xs)==len(coords) and all(isinstance(x,int) and 0<=x<=den for x in xs);assert cand['source_cell_coordinates_lexicographic']==[list(x) for x in coords]
  density=sum((v*Q(x,den) for v,x in zip(vol,xs)),Q(0));assert density==Q(cand['source_density_fraction'])>0;union=set()
  for rect in cand['certified_selected_rectangles']:
   axes=rect['axis_origin_length_cells'];assert len(axes)==dim;assert all(len(x)==2 and all(isinstance(y,int) and not isinstance(y,bool) for y in x) and 0<=x[0]<mx and 1<=x[1]<=mx for x in axes);rc=[tuple((int(o)+j)%mx for j in range(int(length))) for o,length in axes];ids={index[x] for x in itertools.product(*rc)};assert ids==set(rect['receiver_cell_ids']);volume=sum((vol[i] for i in ids),Q(0));mass=sum((vol[i]*Q(xs[i],den) for i in ids),Q(0));mean=mass/volume;assert mean==Q(rect['mean_fraction']);margin=mean-p;assert margin==Q(rect['mean_minus_p_fraction'])>0;minmargin=margin if minmargin is None else min(minmargin,margin);union|=ids;rect_count+=1
  assert union==set(cand['receiver_cell_ids_exact_union']);covered=sum((vol[i] for i in union),Q(0));assert covered==Q(cand['theta_actual_fraction'])>=theta;ratio=covered/density;assert ratio==Q(cand['theta_actual_over_density_fraction']);K=1000*dim;finite=ratio*Q(K-2,K)**dim;log_upper=logi(1/p)[1];orlicz_upper=(1/p)*(1+log_upper)**(dim-1);c_lower=finite/orlicz_upper;physical.append(dict(N=dim,p=str(p),theta_actual=str(covered),source_density=str(density),periodic_halo_lower=str(ratio),finite_box_periods_per_axis=K,finite_physical_indicator_lower=str(finite),Orlicz_denominator_upper=str(orlicz_upper),finite_endpoint_C_lower=str(c_lower),finite_endpoint_C_lower_float=float(c_lower),solver_status=case['status'],solver_optimality_certified=False))
 # A genuinely independent finite LP primal/dual equality and its Gram relaxation.
 p=Q(3,4);f=(Q(1,2),Q(1),Q(1,2));cost=sum(f)/3;mass=(f[0]+f[1])/3;dual=p*Q(4,3)-Q(1,3);energy=(Q(3,4)**2+Q(3,2)**2+Q(3,4)**2)/3;assert cost==dual==Q(2,3) and mass==p*Q(2,3) and energy==Q(9,8) and p*p/energy==Q(1,2)
 sources=['information96.py','information96.wl','periodic96.py','check96.py','math_guard.py','sliver79.py'];inputs=['information/information96.json','information/information_wolfram96.json','python/periodic96.json','resources96.json'];docs=['core96.md','audit96.md','report.md','certification96.md','sources96.md'];prior=['rounds/r011/report.md','rounds/r021/report.md','rounds/r092/core92.md','rounds/r095/core95.md','rounds/r095/audit95.md']
 out=dict(pass_all=True,goal_status='active',general_upper_closed=False,finite_onset_certified=False,strict_information_coefficient_lower=best['information_coefficient_lower'],information_coefficient_exceeds='0.12324',integer_moment_symbolic_residuals=wf['symbolic_residuals'],rational_Jensen_best_and_base_replayed=True,periodic_source_candidates=12,strict_rectangles_replayed=rect_count,minimum_strict_mean_margin=str(minmargin),finite_LP_primal_dual='2/3',finite_Gram_lower='1/2',physical_certificates=physical,source_sha256={x:sha(root/'code'/x) for x in sources},input_sha256={x:sha(data/x) for x in inputs},document_sha256={x:sha(data/x) for x in docs},prior_dependency_sha256={x:sha(root/x) for x in prior},scope='Strict original periodic sources and explicit finite physical/Orlicz lower certificates; directed information necessary coefficient. Solver bounds remain numerical, not continuum or finite-grid optimality certificates, and no general upper closure.')
 Path(args.output).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('physical_certificates','source_sha256','input_sha256','document_sha256','prior_dependency_sha256')},indent=2))
if __name__=='__main__':main()
