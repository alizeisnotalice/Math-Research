"""Read-only local exact verification, independent of numerical union values."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
import round011_geometry as g
R=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mul(a,b):
 z=[x*y for x in a for y in b];return min(z),max(z)
def nearest(x):
 q=x*10**12+F(1,2);return F(q.numerator//q.denominator,10**12)
def main():
 d=json.loads((R/'round011d_direction.json').read_text());h=F(d['s0']);a=list(map(F,d['period_velocity']));v=[list(map(F,row)) for row in d['phase_velocity']];delta=F(d['delta']);assert delta==g.WTS[2] and 0<h<=F(1,4);assert d['steps']==['0',str(h/2),str(-h/2),str(h),str(-h)]
 piv=[(g.P[j]*(1-abs(a[j])*h),g.P[j]*(1+abs(a[j])*h)) for j in range(3)];phiv=[[(g.PH[b][j]-abs(v[b][j])*h,g.PH[b][j]+abs(v[b][j])*h) for j in range(3)] for b in range(3)];assert all(q[0]>1 for q in piv) and all(0<=q[0]<=q[1]<=1 for row in phiv for q in row)
 for s in [-h,h]:assert g.WTS[1]+delta*s>0 and g.WTS[2]-delta*s>0 and sum([g.WTS[0],g.WTS[1]+delta*s,g.WTS[2]-delta*s])==g.W
 mindiff=None;data=[]
 for key in g.SEED['keys']:
  lo,hi,n=g.block_counts(key);act=[b for b in range(3) if n[b]];bh=[max(act,key=lambda b:hi[b][j]+g.PH[b][j]) for j in range(3)];bl=[min(act,key=lambda b:lo[b][j]+g.PH[b][j]) for j in range(3)];m=sum(n[b]*g.WTS[b] for b in act);dm=delta*(n[1]-n[2]);assert m-abs(dm)*h>0;rl,rh=g.radius_interval(m);sp=(dm/(24*rh**2),dm/(24*rl**2)) if dm>=0 else (dm/(24*rl**2),dm/(24*rh**2));box=[]
  for j in range(3):
   for inds,b,sgn in [(hi,bh[j],1),(lo,bl[j],-1)]:
    for c in act:
     if b==c:continue
     gap=sgn*(inds[b][j]+g.PH[b][j]-inds[c][j]-g.PH[c][j])-h*abs(v[b][j]-v[c][j]);assert gap>0;mindiff=gap if mindiff is None else min(mindiff,gap)
   b=bh[j];c=bl[j];lc=mul(piv[j],(hi[b][j]+phiv[b][j][0],hi[b][j]+phiv[b][j][1]));hc=mul(piv[j],(lo[c][j]+phiv[c][j][0],lo[c][j]+phiv[c][j][1]));rmin=g.radius_interval(m-abs(dm)*h)[0];rmax=g.radius_interval(m+abs(dm)*h)[1];box.append(((lc[0]-rmax,lc[1]-rmin),(hc[0]+rmin,hc[1]+rmax)))
  data.append((lo,hi,n,bh,bl,m,sp,box))
 assert str(mindiff)==d['guards']['minimum_normalized_atom_winner_margin'];assert d['guards']['all_original_raw_winner_keys_checked']==639
 contact=[]
 for cc,stored in zip(d['constraints'],d['guards']['contacts']):
  j=cc['axis'];B=data[cc['lower_key_id']];A=data[cc['upper_key_id']];b=B[3][j];c=A[4][j];assert b==cc['lower_owner'] and c==cc['upper_owner'] and b!=c
  coeff=[g.P[j]*(B[1][b][j]+g.PH[b][j]-A[0][c][j]-g.PH[c][j]),g.P[j]*(int(b==1)-int(c==1)),g.P[j]*(int(b==2)-int(c==2))];beta=[B[6][k]+A[6][k] for k in range(2)];assert list(map(str,coeff))==cc['coefficients'] and list(map(str,beta))==cc['beta_interval'];u=[a[j],v[1][j],v[2][j]];norm=sum(x*x for x in coeff);assert u==[nearest(sum(beta)*x/(2*norm)) for x in coeff];dot=sum(x*y for x,y in zip(coeff,u));assert [str(dot-beta[1]),str(dot-beta[0])]==cc['residual_interval'];assert not(beta[0]<=0<=beta[1])
  marg=min(B[7][j][0][0],piv[j][0]-B[7][j][0][1],A[7][j][1][0],piv[j][0]-A[7][j][1][1]);tang=[min(A[7][t][1][0],B[7][t][1][0])-max(A[7][t][0][1],B[7][t][0][1]) for t in range(3) if t!=j];assert marg>0 and min(tang)>0;assert str(marg)==stored['unclipped_margin_lower'] and list(map(str,tang))==stored['tangent_overlap_lower'];contact.append({'axis':j,'raw_face_guard':True,'positive_tangent_overlap':True,'residual_interval':cc['residual_interval']})
 lock=json.loads((R/'round011d_query_lock.json').read_text());assert all(sha(R/name)==digest for name,digest in lock.items());family=json.loads((R/'round011d_family.json').read_text());assert set(map(tuple,g.SEED['keys']))<=set(map(tuple,family['keys']));res=json.loads((R/'round011d_results.json').read_text());assert len(res['queries'])==4;assert [q['s'] for q in res['queries']]==d['steps'][1:];assert all(sum(map(F,q['parameters']['weights']))==g.W for q in res['queries'])
 out={'status':'exact local-direction and pre-query-lock verification PASS','input_sha256':sha(R/'round008_certificate_input.json'),'direction_sha256':sha(R/'round011d_direction.json'),'verified_keys':639,'mass_per_cell_exact_constant':str(g.W),'delta':str(delta),'all_original_raw_extrema_keep_unique_winners':True,'whole_interval_period_phase_positive_mass_guards':True,'minimum_atom_winner_normalized_margin':str(mindiff),'contact_constraints':contact,'quantized_compensation_not_exact_zero':True,'pre_query_source_lock_unchanged':True,'old639_in_full_family':True,'new_full_joint_union_parameter_queries':4,'baseline_union_queries':0,'scope':'local exact analysis only; numerical union values not certified; no entire-direction exclusion'}
 (R/'round011d_local_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='contact_constraints'},indent=2))
if __name__=='__main__':main()
