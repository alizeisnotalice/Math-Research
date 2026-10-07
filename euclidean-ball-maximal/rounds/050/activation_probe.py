#!/usr/bin/env python3
"""Exact retained-activation bridge and future-occupation probes (n=1)."""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left,bisect_right
import argparse,importlib.util,json,time,hashlib,sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
LABELS={'minimal_four_atoms','four_atom_perturbation_s471007','four_atom_perturbation_s471016','minimal_cloud_with_two_coarse_atoms','geometric_contact_clouds_s491101','geometric_contact_clouds_s491102','geometric_contact_clouds_s491108','chain_L2_alpha1','chain_L4_alpha1','chain_L8_alpha1','chain_L16_alpha1','chain_L8_alpha2'}

def module(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def outside(shear,poly,R):
 return sum((shear.area(p) for p in (shear.clip(poly,F(-1),F(1),-R),shear.clip(poly,F(1),F(-1),-R))),F(0))

def evaluate(shear,globalmod,flow,old,batch,label,atoms,alpha,params):
 began=time.monotonic();saved,radii,cells,_=flow.context(old,atoms,alpha);lo=alpha/4;hi=alpha/2;I=hi-lo;X=alpha*F(saved['eligible_volume'])
 targets,classified=shear.target_cells(old,atoms,alpha,radii,saved);globalcells=[]
 for a,b,kind,_ in classified:
  mid=(a+b)/2;trace=[sum((w for y,w in atoms if abs(y-mid)<R),F(0))/(2*R) for R in radii[1:]]
  rec,infinity=globalmod.global_records(trace,lo,hi);globalcells.append((a,b,kind,rec))
 raw_global_count=len(globalcells);merged=[]
 for a,b,kind,rec in globalcells:
  if merged and merged[-1][1]==a and merged[-1][2:]==(kind,rec):
   previous=merged[-1];merged[-1]=(previous[0],b,kind,rec)
  else:merged.append((a,b,kind,rec))
 globalcells=merged;targets=[globalcells[0][0]]+[b for a,b,kind,rec in globalcells]
 data=[]
 for a,b,J,trace in cells:
  rec,infinity=globalmod.global_records(trace,lo,hi);assert infinity==0
  mid=(a+b)/2;members={k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k]) for k in rec}
  data.append((a,b,trace,rec,members))
 raw_data=data;data=[]
 for a,b,trace,rec,members in raw_data:
  if data and data[-1][1]==a and data[-1][3:]==(rec,members):
   previous=data[-1];data[-1]=(previous[0],b,previous[2],rec,members)
  else:data.append((a,b,trace,rec,members))
 C_D=C_G=C_D_fine=T=T_outer=T_fine=T_direct=R_cross=B=original_O=F(0);target_pieces=triples=T_terms=0;CD_by_kind={s:F(0) for s in ('eligible_E','MP_le_alpha','MP_gt_2alpha','band_J1')}
 CD_l_eq_j=CD_l_gt_j=F(0)
 for xa,xb,tx,ix,mx in data:
  for za,zb,tz,iz,mz in data:
   rectangle=[(xa,za),(xb,za),(xb,zb),(xa,zb)]
   for j,(ja,jb) in ix.items():
    R=radii[j]
    for k,(ka,kb) in iz.items():
     if j>=k:continue
     bl=max(ja,ka);br=min(jb,kb)
     if bl>=br:continue
     r=radii[k];pairwidth=(br-bl)/I
     pairstrip=globalmod.strip_area(xa,xb,za,zb,R);B+=pairstrip*tz[k-1]*pairwidth/(2*R)
     for source in mx[j]&mz[k]:
      y,w=atoms[source];density=w/(4*R*r)
      R_cross+=(xb-xa)*(zb-za)*density*pairwidth
      original_O+=((xb-xa)*(zb-za)-pairstrip)*density*pairwidth
      # Since 2^(j-l)u_l=u_j, retaining D_l cancels the scale coefficient.
      # T removes F_l and sums every l in [j,k), still with same x/z beta.
      for l in range(j,k):
       Rl=radii[l];xl=max(xa,y-Rl);xr=min(xb,y+Rl);A=max(F(0),xr-xl)*(zb-za)
       if not A:continue
       T_terms+=1;amount=A*density*pairwidth;T+=amount
       outerA=A-globalmod.strip_area(xl,xr,za,zb,R);T_outer+=outerA*density*pairwidth
       finearea=max(F(0),min(xr,y+r)-max(xl,y-r))*(zb-za);T_fine+=finearea*density*pairwidth
       if l>j:assert outerA==0
       # Independent product-convolution integration after dropping F_l.
       xlength=max(F(0),min(xb,y+Rl)-max(xa,y-Rl));zlength=max(F(0),min(zb,y+r)-max(za,y-r))
       T_direct+=w*F(2)**(j-l)*xlength*zlength/(4*Rl*r)*pairwidth
      # C_D keeps the q first-record condition and integrates the same shear.
      qmin=xa+za-y;qmax=xb+zb-y;start=max(0,bisect_right(targets,qmin)-1);stop=bisect_left(targets,qmax)
      assert targets[0]<=qmin and qmax<=targets[-1]
      for a,b,kind,iq in globalcells[start:stop]:
       if b<=qmin or a>=qmax:continue
       common=[(l,max(bl,la),min(br,lb)) for l,(la,lb) in iq.items() if j<=l<k and max(bl,la)<min(br,lb)]
       if not common:continue
       sub=shear.slab(rectangle,y,a,b);A=shear.area(sub)
       if not A:continue
       target_pieces+=1
       for l,u,v in common:
        poly=shear.x_slab(sub,y,radii[l]);PA=shear.area(poly)
        if not PA:continue
        triples+=1;factor=density*(v-u)/I;amount=PA*factor;C_D+=amount;CD_by_kind[kind]+=amount
        C_G+=outside(shear,poly,R)*factor;C_D_fine+=shear.area(shear.x_slab(poly,y,r))*factor
        if l==j:CD_l_eq_j+=amount
        else:CD_l_gt_j+=amount
 assert T==T_direct and C_G<=C_D<=T and C_D_fine<=T_fine and T_outer+T_fine<=T and C_G<=T_outer
 assert C_D<=C_G+B and C_D<=R_cross<=T
 assert T_outer==original_O and T-R_cross<=F(9,8)*X
 assert sum(CD_by_kind.values(),F(0))==C_D and CD_l_eq_j+CD_l_gt_j==C_D
 # Independent first-moment convolution for the entire retained activation.
 activation=F(0)
 for a,b,trace,ix,members in raw_data:
  for j,(ja,jb) in ix.items():
   for l in range(j,saved['D']+1):
    activation+=(b-a)*trace[l-1]*F(2)**(j-l)*(jb-ja)/I
 assert activation<=F(5,2)*X
 return dict(batch=batch,label=label,params=params,N=len(atoms),D=saved['D'],alpha=str(alpha),X=str(X),C_D=str(C_D),T=str(T),C_G=str(C_G),R=str(R_cross),B=str(B),original_O=str(original_O),R_over_X=str(R_cross/X),B_over_X=str(B/X),T_over_R=str(T/R_cross) if R_cross else None,extra_T_over_X=str((T-R_cross)/X),deep_budget_check='T-R<=9X/8 passed',outer_identity_check='T_outer=original_O passed',retained_bridge_checks='C_G<=C_D<=C_G+B and C_D<=R<=T passed',C_D_over_X=str(C_D/X),T_over_X=str(T/X),C_G_over_X=str(C_G/X),T_outer=str(T_outer),T_outer_over_X=str(T_outer/X),T_source_inner_fine=str(T_fine),T_source_inner_fine_over_X=str(T_fine/X),T_inside_coarse_outside_source_fine=str(T-T_outer-T_fine),T_inside_coarse_outside_source_fine_over_X=str((T-T_outer-T_fine)/X),C_D_source_inner_fine=str(C_D_fine),C_D_source_inner_fine_over_X=str(C_D_fine/X),C_D_inside_coarse=str(C_D-C_G),C_D_inside_coarse_over_X=str((C_D-C_G)/X),target_gate_loss=str(T-C_D),target_gate_loss_over_X=str((T-C_D)/X),C_D_by_target_kind={s:str(v) for s,v in CD_by_kind.items()},C_D_same_old_label=str(CD_l_eq_j),C_D_later_old_label=str(CD_l_gt_j),activation_first_moment=str(activation),activation_first_moment_over_X=str(activation/X),independent_T_product_check='passed',observer_cells_original=len(raw_data),observer_cells_merged=len(data),global_cells_original=raw_global_count,global_cells_merged=len(globalcells),target_pieces=target_pieces,three_record_intersections=triples,T_terms=T_terms,source_weight_preserved=True,exact_error='0',elapsed_seconds=time.monotonic()-began)

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);p.add_argument('--labels',nargs='*');a=p.parse_args()
 shear=module(a.repo_root/'rounds/048/shear_probe.py','shear50');glob=module(a.repo_root/'rounds/049/global_probe.py','global50');flow=module(a.repo_root/'rounds/047/flow_probe.py','flow50');record=module(a.repo_root/'rounds/046/threshold_probe.py','record50');prefix=module(a.repo_root/'rounds/046/prefix_probe.py','prefix50');old=record.load(a.repo_root,a.output)
 receipt=json.loads((a.repo_root/'rounds/049/verification.json').read_text());reference={c['label']:c for c in receipt['global_records']['cases']}
 chosen=LABELS if a.labels is None else set(a.labels);inputs=list(flow.cases(old,prefix));inputs.append((0,'minimal_cloud_with_two_coarse_atoms',[(F(-507,256),F(3,10)),(F(0),F(13,256)),(F(1,16),F(3,64)),(F(3,32),F(9,128)),(F(515,256),F(3,10)),(F(10),F(297,1280))],F(1),dict(family='two_coarse_atoms_for_J1_target')))
 began=time.monotonic();out=dict(status='running',scope='actual n1 retained-activation bridge with shared threshold and original source',cases=[],finite_tests_prove_uniform_bound=False)
 for batch,label,atoms,alpha,params in inputs:
  if label not in chosen:continue
  c=evaluate(shear,glob,flow,old,batch,label,atoms,alpha,params)
  independent=flow.flow_probe(old,record,batch,label,atoms,alpha,params);assert F(c['B'])==F(independent['mean']['B']);assert F(c['R'])==(F(independent['mean']['Q'])-F(independent['mean']['diag']))/2;c['independent_flow_R_B_check']='passed'
  if label in reference:
   ref=reference[label];expected=F(ref['O_over_X'])*F(ref['new_paid_over_O']) if ref['new_paid_over_O'] is not None else F(0);assert F(c['C_G_over_X'])==expected;c['independent_round49_CG_crosscheck']='passed'
  out['cases'].append(c);a.output.write_text(json.dumps(out)+'\n');print(batch,label,'C_D/X',float(F(c['C_D_over_X'])),'T/X',float(F(c['T_over_X'])),'C_G/X',float(F(c['C_G_over_X'])),'sec',round(c['elapsed_seconds'],2),flush=True)
 out.update(status='passed',case_count=len(out['cases']),batch_counts=[sum(c['batch']==b for c in out['cases']) for b in range(3)],elapsed_seconds=time.monotonic()-began,sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),exact_error='0');a.output.write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['case_count'],out['batch_counts'],out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
