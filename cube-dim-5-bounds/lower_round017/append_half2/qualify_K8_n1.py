#!/usr/bin/env python3
"""One analytic nearest-point qualification; no cell arrangement or union."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import hashlib,importlib.util,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;PREV=HERE.parent;PROJECT=PREV.parent
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();assert all(os.environ.get(k)=='1' for k in ENV);assert os.getpriority(os.PRIO_PROCESS,0)>=10
 reg=json.loads((HERE/'registration_qualify_K8_n1.json').read_text());assert reg['case']==[8,1]
 hp=PREV/'endpoint_diagnostic.py';assert sha(hp)==reg['exact_helper_sha256'];spec=importlib.util.spec_from_file_location('round017_exact_helpers_qualification',hp);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
 old=json.loads((PREV/'endpoint_diagnostic_results.json').read_text());pending=[tuple(c) for c in old['new_planar_K_n_cases_needed_for_window']];assert len(pending)==18 and (8,1) in pending
 ul,uh=h.roots(2);ss=((1+ul)/2,(1+uh)/2);side=h.roots(8);rr=(side[0]/2,side[1]/2);rlo,rhi=rr
 assert (2*rlo)**5<8<=(2*rhi)**5 and Q(3,2)<side[0]<=side[1]<2
 assert ss[1]/2<rlo<rhi<ss[0]
 corners=list(product((Q(-1,2),Q(1,2)),repeat=2));max_error=Q(0);checks=0
 for em,en in corners:
  ex,ey=(4*em-3*en)/5,(3*em+4*en)/5
  assert abs(ex)<=Q(7,10) and abs(ey)<=Q(7,10);checks+=2;max_error=max(max_error,abs(ex),abs(ey))
 assert max_error==Q(7,10)
 cover_upper=max_error*ss[1];assert cover_upper<rlo
 # At the left parameter endpoint, verify every unrotated vertex explicitly.
 amin=Q(3,2);spaces=((amin,amin),ss,ss);box=[]
 for ab in spaces:box.extend((max(Q(0),ab[1]-rlo),min(ab[0]/2,rlo)))
 assert all(box[2*q]<box[2*q+1] for q in range(3));vertices=list(product((0,1),repeat=3));assert len(vertices)==len(set(vertices))==8
 for v in vertices:
  for q,(j,ab) in enumerate(zip(v,spaces)):
   pl,ph=j*ab[0],j*ab[1];assert ph-rlo<=box[2*q] and box[2*q+1]<=pl+rlo;checks+=2
 ratio_lo=8*ss[0]**2*(rlo-amin/2)*(rlo-ss[1]/2)**2
 ratio_hi=8*ss[1]**2*(rhi-amin/2)*(rhi-ss[0]/2)**2
 cp=PROJECT/'lower_round003/a2_equals_a3_results.json';C=json.loads(cp.read_text());assert sha(cp)==old['C_source_sha256'];assert 0<ratio_lo<=ratio_hi<Q(C['ratio_interval'][0])
 remaining=[list(c) for c in pending if c!=(8,1)]
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 out={'status':'PASS_one_analytic_planar_cover_qualification','K':8,'n':1,'exact_common_radius':'8^(1/5)/2=u^3/2','strict_inner_radius':str(rlo),'selected_mass':8,
      'nearest_point_cover_proof':'round both original rotated-lattice coordinates to integers; error lies in[-1/2,1/2]^2 and R maps its corners to supnorm<=7/10; multiply by s',
      'supnorm_cover_radius_upper':str(cover_upper),'strict_radius_margin_lower':str(rlo-cover_upper),'plane_capture_set':'entire full[0,5s]^2 period even at the strict inner radius','plane_area_exact_expression':'25*s^2=(25/4)+(25/2)*u+(25/4)*u^2','plane_area_interval':[str(25*ss[0]**2),str(25*ss[1]**2)],
      'first_axis_activation_parameter_interval':[str(t) for t in side],'single_family_ratio_expression':'8*s^2*(r-a/2)*(r-s/2)^2 for3/2<=a<2r;0 for a>=2r','single_family_maximum_at_a':['3/2'],'single_family_maximum_ratio_interval':[str(ratio_lo),str(ratio_hi)],'single_family_alone_below_C_cert':True,
      'unrotated_inner_box_at_a_3_2':[str(t) for t in box],'distinct_unrotated_vertices':8,'exact_corner_and_endpoint_inequalities_checked':checks,
      'original_pending_materialized_plane_packets':18,'analytically_resolved_capture_cases':[[8,1]],'remaining_unresolved_capture_count':len(remaining),'remaining_unresolved_capture_cases':remaining,
      'fixed_atom_ID_partition_materialized':False,'downstream_usage':'use wholeplane with its nearest-point selector proof; do not assignonefixedatomID totheentireperiod; materialize a fixed-IDpartition only if a laterpipeline requiresit',
      'source_sha256':sha(Path(__file__)),'previous_result_sha256':sha(PREV/'endpoint_diagnostic_results.json'),'C_source_sha256':sha(cp),'nodes':checks,'runtime':{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}},
      'limits':{'seconds':120,'nodes':200000},'scope':'one previouslyunevaluated planarcase only; no arrangement scan,fullunion,newparameter grid or lowerrecord'}
 assert checks<200000 and time.monotonic()-start<120
 (HERE/'qualification_K8_n1_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
