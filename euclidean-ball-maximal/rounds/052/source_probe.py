#!/usr/bin/env python3
"""Exact actual one-dimensional outer cost split by coarse-source boundary depth."""
from pathlib import Path
from fractions import Fraction as F
import argparse,importlib.util,json,sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
LABELS={'minimal_four_atoms','four_atom_perturbation_s471007','four_atom_perturbation_s471016','geometric_contact_clouds_s491101','geometric_contact_clouds_s491102','geometric_contact_clouds_s491108','chain_L4_alpha1','chain_L8_alpha1','chain_L8_alpha2'}
def module(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def evaluate(flow,record,globalmod,old,batch,label,atoms,alpha):
 saved,radii,cells,_=flow.context(old,atoms,alpha);lo=alpha/4;hi=alpha/2;X=alpha*F(saved['eligible_volume']);data=[]
 for a,b,J,trace in cells:
  rec=record.intervals(trace,lo,hi);mid=(a+b)/2
  members={k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k]) for k in rec}
  if data and data[-1][1]==a and data[-1][2:]==(rec,members):
   previous=data[-1];data[-1]=(previous[0],b,rec,members)
  else:data.append((a,b,rec,members))
 O=F(0);rhos=[1-F(1,2**b) for b in range(1,5)];deep={rho:F(0) for rho in rhos};terms=0
 for xa,xb,ix,mx in data:
  for za,zb,iz,mz in data:
   for j,(ja,jb) in ix.items():
    R=radii[j]
    for k,(ka,kb) in iz.items():
     if k<=j:continue
     width=min(jb,kb)-max(ja,ka)
     if width<=0:continue
     r=radii[k];outside=(xb-xa)*(zb-za)-globalmod.strip_area(xa,xb,za,zb,R)
     for source in mx[j]&mz[k]:
      y,w=atoms[source];density=w*width/((hi-lo)*4*R*r)
      O+=density*outside;terms+=1
      for rho in rhos:
       a=max(xa,y-rho*R);b=min(xb,y+rho*R)
       if a<b:deep[rho]+=density*((b-a)*(zb-za)-globalmod.strip_area(a,b,za,zb,R))
 assert deep[F(1,2)]==0
 assert list(deep.values())==sorted(deep.values())
 rows=[]
 for b,rho in enumerate(rhos,1):
  A=sum((max(F(0),F(1,2)*(1-F(2)**d*(1-rho))) for d in range(1,b+1)),F(0))
  A1=max(F(0),rho-F(1,2));G=(A+A1)/2;assert G==F(b-1,4)
  bound=2*X*rho*G
  assert 0<=deep[rho]<=O and deep[rho]<=bound
  rows.append(dict(rho=str(rho),inner_cost=str(deep[rho]),inner_over_X=str(deep[rho]/X),inner_over_O=str(deep[rho]/O) if O else None,shell_cost=str(O-deep[rho]),exact_upper_over_X=str(2*rho*G),old_triangle_upper_over_X=str(2*rho*(b-1)),fine_sup_escape=str(G)))
 return dict(batch=batch,label=label,D=saved['D'],N=len(atoms),X=str(X),O=str(O),O_over_X=str(O/X),source_terms=terms,profiles=rows)

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 flow=module(a.repo_root/'rounds/047/flow_probe.py','flow52');rec=module(a.repo_root/'rounds/046/threshold_probe.py','record52');pre=module(a.repo_root/'rounds/046/prefix_probe.py','prefix52');glob=module(a.repo_root/'rounds/049/global_probe.py','global52');old=rec.load(a.repo_root,a.output)
 refs={c['label']:c for c in json.loads((a.repo_root/'rounds/049/verification.json').read_text())['global_records']['cases']};cases=[]
 for batch,label,atoms,alpha,params in flow.cases(old,pre):
  if label not in LABELS:continue
  c=evaluate(flow,rec,glob,old,batch,label,atoms,alpha);assert F(c['O_over_X'])==F(refs[label]['O_over_X']);c['independent_round49_O_check']='passed';cases.append(c)
 assert len(cases)==9
 out=dict(status='passed',batch_counts=[sum(c['batch']==b for c in cases) for b in range(3)],cases=cases,scope='Exact n=1 source/observer/common-threshold integrals; no higher-dimensional uniform conclusion from finite cases.',main_theorem_proved=False)
 a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status='passed',cases=len(cases),source_terms=sum(c['source_terms'] for c in cases))))
if __name__=='__main__':main()
