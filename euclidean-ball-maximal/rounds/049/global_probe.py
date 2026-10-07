#!/usr/bin/env python3
"""Exact full-space first-threshold records for source-preserving shear.
Only requested external receipts are written; original P and common beta retained.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from bisect import bisect_left,bisect_right
import argparse,importlib.util,json,time,hashlib,sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
EXTRA={'four_atom_perturbation_s471007','four_atom_perturbation_s471010','four_atom_perturbation_s471013','four_atom_perturbation_s471017','narrow_multiple_clouds_s481104','narrow_multiple_clouds_s481109','geometric_contact_clouds_s491102','geometric_contact_clouds_s491107','chain_L2_alpha1','chain_L4_alpha1','chain_L8_alpha1','chain_L8_alpha3/2','chain_L32_alpha2'}

def module(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def global_records(trace,lo,hi):
 """K=infinity is allowed: finite labels need not exhaust I."""
 running=F(0);result={}
 for l,g in enumerate(trace,1):
  a=max(lo,running);b=min(hi,g)
  if a<b:result[l]=(a,b)
  running=max(running,g)
 infinity=max(F(0),hi-max(lo,running))
 assert sum((b-a for a,b in result.values()),F(0))+infinity==hi-lo
 # Independent scalar first-exit check, including non-exit.
 cuts=sorted({lo,hi}|{g for g in trace if lo<g<hi});direct=defaultdict(F)
 for a,b in zip(cuts,cuts[1:]):
  beta=(a+b)/2;l=next((k for k,g in enumerate(trace,1) if g>beta),None);direct[l]+=b-a
 assert direct[None]==infinity and all(direct[l]==b-a for l,(a,b) in result.items())
 assert sum(direct.values())==hi-lo
 return result,infinity

def strip_area(xa,xb,za,zb,R):
 square=lambda t:max(F(0),t)**2
 value=sum((sign*(square(zb+shift)-square(za+shift))/2 for shift,sign in ((R-xa,1),(R-xb,-1),(-R-xa,-1),(-R-xb,1))),F(0))
 assert 0<=value<=(xb-xa)*(zb-za)
 return value

def evaluate(shear,flow,old,batch,label,atoms,alpha,params):
 began=time.monotonic();saved,radii,cells,_=flow.context(old,atoms,alpha);lo=alpha/4;hi=alpha/2;I=hi-lo;X=alpha*F(saved['eligible_volume'])
 targets,classified=shear.target_cells(old,atoms,alpha,radii,saved);globalcells=[]
 for a,b,kind,_ in classified:
  mid=(a+b)/2;trace=[sum((w for y,w in atoms if abs(y-mid)<R),F(0))/(2*R) for R in radii[1:]]
  rec,infinity=global_records(trace,lo,hi);globalcells.append((a,b,kind,trace,rec,infinity))
 data=[]
 for a,b,J,trace in cells:
  rec,infinity=global_records(trace,lo,hi);assert infinity==0
  mid=(a+b)/2;members={k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k]) for k in rec}
  data.append((a,b,trace,rec,members))
 # Global comparison directly integrates every global first-record interval.
 Bglobal=Beligible=F(0);comparison_pairs=0
 for a,b,kind,trace,iq,infinity in globalcells:
  if not iq:continue
  for za,zb,tz,iz,mz in data:
   for l,(la,lb) in iq.items():
    A=strip_area(a,b,za,zb,radii[l])
    if not A:continue
    for k,(ka,kb) in iz.items():
     if l>=k:continue
     W=max(F(0),min(lb,kb)-max(la,ka))
     if not W:continue
     value=A*tz[k-1]*W/(2*radii[l]*I);Bglobal+=value;comparison_pairs+=1
     if kind=='eligible_E':Beligible+=value
 # Exact shear classification against full-space target records.
 kinds=('eligible_E','MP_le_alpha','MP_gt_2alpha','band_J1');classes={k:F(0) for k in kinds};paidclass={k:F(0) for k in kinds};sameclass={k:F(0) for k in kinds}
 relations={k:F(0) for k in ('Kq_lt_j','Kq_eq_j','j_lt_Kq_lt_k','Kq_eq_k','Kq_gt_k','Kq_infinity')}
 O=paid=same=F(0);polygons=pieces=record_pairs=0;unpaid_witness=None
 for xa,xb,tx,ix,mx in data:
  for za,zb,tz,iz,mz in data:
   rectangle=[(xa,za),(xb,za),(xb,zb),(xa,zb)]
   for j,(ja,jb) in ix.items():
    R=radii[j];polys=[p for p in (shear.clip(rectangle,F(-1),F(1),-R),shear.clip(rectangle,F(1),F(-1),-R)) if shear.area(p)>0]
    assert sum((shear.area(p) for p in polys),F(0))==(xb-xa)*(zb-za)-strip_area(xa,xb,za,zb,R)
    for k,(ka,kb) in iz.items():
     if j>=k:continue
     bl=max(ja,ka);br=min(jb,kb)
     if bl>=br:continue
     r=radii[k]
     for source in mx[j]&mz[k]:
      y,w=atoms[source];density=w/(4*R*r)
      for poly in polys:
       polygons+=1;cost=shear.area(poly)*density*(br-bl)/I;O+=cost;subtotal=F(0)
       qmin=min(x+z-y for x,z in poly);qmax=max(x+z-y for x,z in poly);start=max(0,bisect_right(targets,qmin)-1);stop=bisect_left(targets,qmax)
       for a,b,kind,tq,iq,infinity in globalcells[start:stop]:
        if b<=qmin or a>=qmax:continue
        sub=shear.slab(poly,y,a,b);A=shear.area(sub)
        if not A:continue
        pieces+=1;value=A*density*(br-bl)/I;classes[kind]+=value;subtotal+=value;finite=F(0);localpaid=F(0)
        for l,(la,lb) in iq.items():
         u=max(bl,la);v=min(br,lb)
         if u>=v:continue
         record_pairs+=1;factor=density*(v-u)/I;finite+=A*factor
         relation='Kq_lt_j' if l<j else 'Kq_eq_j' if l==j else 'j_lt_Kq_lt_k' if l<k else 'Kq_eq_k' if l==k else 'Kq_gt_k';relations[relation]+=A*factor
         if j<=l<k:
          compA=shear.area(shear.x_slab(sub,y,radii[l]));amount=compA*factor;paid+=amount;paidclass[kind]+=amount;localpaid+=amount
          if l==j:same+=amount;sameclass[kind]+=amount
        nonexit=value-finite;assert nonexit>=0;relations['Kq_infinity']+=nonexit
        if value>localpaid and unpaid_witness is None:
         unpaid_witness=dict(source_y=str(y),j=j,k=k,beta_lo=str(bl),beta_hi=str(br),target_kind=kind,target_interval=[str(a),str(b)],piece_cost=str(value),piece_paid=str(localpaid),polygon=[[str(x),str(z)] for x,z in sub],target_records={str(l):[str(u),str(v)] for l,(u,v) in iq.items()},target_nonexit_width=str(infinity))
       assert subtotal==cost
 assert sum(classes.values(),F(0))==O and sum(relations.values(),F(0))==O and sum(paidclass.values(),F(0))==paid
 assert paid<=Bglobal and same<=Bglobal and Beligible<=Bglobal
 remainder=O-paid;assert remainder>=0
 escaped=O-classes['eligible_E'];caught_escape=paid-paidclass['eligible_E']
 return dict(batch=batch,label=label,params=params,N=len(atoms),D=saved['D'],alpha=str(alpha),X=str(X),O=str(O),O_over_X=str(O/X),Bglobal=str(Bglobal),Bglobal_over_X=str(Bglobal/X),Beligible=str(Beligible),Beligible_over_X=str(Beligible/X),new_paid=str(paid),new_paid_over_O=str(paid/O) if O else None,new_paid_over_X=str(paid/X),same_label=str(same),remainder=str(remainder),remainder_over_O=str(remainder/O) if O else None,original_classes={k:str(v) for k,v in classes.items()},new_paid_by_class={k:str(v) for k,v in paidclass.items()},same_label_by_class={k:str(v) for k,v in sameclass.items()},caught_fraction_by_class={k:str(paidclass[k]/classes[k]) if classes[k] else None for k in kinds},original_escape_cost=str(escaped),caught_original_escape=str(caught_escape),caught_original_escape_fraction=str(caught_escape/escaped) if escaped else None,global_record_relation={k:str(v) for k,v in relations.items()},global_cell_count=len(globalcells),global_nonexit_cells=sum(infinity>0 for a,b,kind,trace,rec,infinity in globalcells),comparison_pairs=comparison_pairs,source_polygons=polygons,target_pieces=pieces,target_record_pairs=record_pairs,unpaid_witness=unpaid_witness,source_weight_preserved=True,exact_error='0',elapsed_seconds=time.monotonic()-began)

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);p.add_argument('--labels',nargs='*');a=p.parse_args()
 shear=module(a.repo_root/'rounds/048/shear_probe.py','shear49');flow=module(a.repo_root/'rounds/047/flow_probe.py','flow49');record=module(a.repo_root/'rounds/046/threshold_probe.py','record49');prefix=module(a.repo_root/'rounds/046/prefix_probe.py','prefix49');old=record.load(a.repo_root,a.output)
 receipt=json.loads((a.repo_root/'rounds/048/verification.json').read_text());reference={c['label']:c for c in receipt['shear']['cases']} if 'shear' in receipt else {}
 # The archived round48 script owns the original nine labels and the sixth-atom input.
 chosen=set(shear.LABELS)|EXTRA if a.labels is None else set(a.labels)
 inputs=list(flow.cases(old,prefix));inputs.append((0,'minimal_cloud_with_two_coarse_atoms',[(F(-507,256),F(3,10)),(F(0),F(13,256)),(F(1,16),F(3,64)),(F(3,32),F(9,128)),(F(515,256),F(3,10)),(F(10),F(297,1280))],F(1),dict(family='two_coarse_atoms_for_J1_target')))
 began=time.monotonic();out=dict(status='running',scope='finite actual n1 full-space common-threshold source-preserving shear',cases=[],finite_tests_prove_uniform_general_bound=False)
 for batch,label,atoms,alpha,params in inputs:
  if label not in chosen:continue
  c=evaluate(shear,flow,old,batch,label,atoms,alpha,params)
  if label in reference:
   ref=reference[label];assert c['O_over_X']==ref['O_over_X'];O=F(c['O'])
   assert {k:str(F(v)/O) if O else None for k,v in c['original_classes'].items()}==ref['class_fraction']
   assert (str(F(c['new_paid_by_class']['eligible_E'])/O) if O else None)==ref['compatible_over_O']
   independent=flow.flow_probe(old,record,batch,label,atoms,alpha,params);assert c['Beligible']==independent['mean']['B'];c['round48_crosscheck']='passed';c['independent_eligible_flow_B_check']='passed'
  out['cases'].append(c);a.output.write_text(json.dumps(out)+'\n');print(batch,label,'Bglobal/X',float(F(c['Bglobal_over_X'])),'paid/O',c['new_paid_over_O'],'catch',c['caught_original_escape_fraction'],'seconds',round(c['elapsed_seconds'],2),flush=True)
 out.update(status='passed',case_count=len(out['cases']),batch_counts=[sum(c['batch']==b for c in out['cases']) for b in range(3)],elapsed_seconds=time.monotonic()-began,sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),exact_error='0');a.output.write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['case_count'],out['batch_counts'],out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
