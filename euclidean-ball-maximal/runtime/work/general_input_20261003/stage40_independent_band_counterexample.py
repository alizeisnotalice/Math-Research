"""Independent exact full-mask and bridge-saturation reconstruction.
No producer imports; stop law is read from saved actual common-c certificate.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,sys
D=Path(__file__).resolve().parents[2]/'output/general_input_20261003'
data=json.loads((D/'stage40_band_search.json').read_text())
name='stage40_band_counterexample'
p=data['retained_complete_scopes'][data['simplest_optimized_positive_case']]
(D/(name+'.json')).write_text(json.dumps(p,indent=2)+'\n')
s=p['summary'];a=p['complete_actual_matched_height_audit'];stop=a['complete_actual_observer_stop']
atoms=[(F(v['location']),F(v['mass'])) for v in p['original_source']]
alpha=F(s['alpha']);aa=alpha/8;b=alpha/2;dd=s['actual_D']+1;R={j:F(1,2**j)/(2*aa) for j in range(1,dd+1)}
assert R[dd]<min(w for y,w in atoms)/(2*alpha)
eta=[(F(v['lo']),F(v['hi']),F(v['density'])) for v in stop['complete_eta_density']]
assert sum(d*(r-l) for l,r,d in eta)==1 and all(0<=d<=2*aa for l,r,d in eta)
assert sum(d*(r*r-l*l)/2 for l,r,d in eta)==sum(w*y for y,w in atoms)
def maximal(x):
 ds=sorted({abs(x-y) for y,w in atoms});assert ds[0]>0
 return max(sum(w for y,w in atoms if abs(x-y)<=r)/(2*r) for r in ds)
knots={y for y,w in atoms}|{y+sg*r for y,w in atoms for r in R.values() for sg in [-1,1]}
for l in range(len(atoms)):
 for r in range(l,len(atoms)):
  mass=sum(w for y,w in atoms[l:r+1]);mid=(atoms[l][0]+atoms[r][0])/2;half=(atoms[r][0]-atoms[l][0])/2
  for beta in [alpha,2*alpha]:
   distance=mass/(2*beta)-half
   if distance>=0:knots.update([mid-distance,mid+distance])
knots=sorted(knots);cells=[];full_volume=F(0)
for l,r in zip(knots,knots[1:]):
 x=(l+r)/2
 if not alpha<maximal(x)<=2*alpha:continue
 full_volume+=r-l
 g={j:sum(w for y,w in atoms if abs(x-y)<R[j])/(2*R[j]) for j in R}
 J=next(j for j in R if g[j]>aa);K=next(j for j in R if g[j]>b)
 if J>=2:cells.append((l,r,J,K))
volume=sum(r-l for l,r,J,K in cells);assert volume==F(a['eligible_volume'])
assert full_volume==F(a['original_full_volume'])
def h(x,k):return sum((max(F(0),min(r,x+R[k])-max(l,x-R[k]))/(2*R[k]) for l,r,J,K in cells if K==k),F(0))
ks=sorted({K for l,r,J,K in cells})
def A(x):return sum(h(x,k) for k in ks)
def u(x):return sum(d*(x*(min(x,r)-l)-(min(x,r)**2-l*l)/2) for l,r,d in eta if x>l)-sum(w*max(x-y,0) for y,w in atoms)
minima=[];pieces=0
for k in ks:
 t=alpha*R[k]**2/16
 nodes=sorted({v for l,r,J,K in cells if K==k for v in [l-R[k],l+R[k],r-R[k],r+R[k]]}|{y for y,w in atoms}|{v for l,r,d in eta for v in [l,r]})
 for l,r in zip(nodes,nodes[1:]):
  mid=(l+r)/2
  if h(mid,k)==0:continue
  d=sum(d for el,er,d in eta if el<mid<er);length=r-l;c=u(l);aa2=d/2;bb=(u(r)-c-aa2*length*length)/length
  assert u(mid)==c+bb*length/2+aa2*length*length/4
  values=[c,u(r)]
  if aa2 and 0<-bb/(2*aa2)<length:
   v=-bb/(2*aa2);values.append(c+bb*v+aa2*v*v)
  minima.append(min(values)/t);pieces+=1
assert min(minima)==F(s['minimum_u_over_t_exact']) and min(minima)>1
nodes=sorted({v for l,r,J,k in cells for v in [l-R[k],l+R[k],r-R[k],r+R[k]]}|{v for l,r,d in eta for v in [l,r]})
source=sum(w*max(F(0),A(y)-1) for y,w in atoms);returned=L=F(0)
for l,r in zip(nodes,nodes[1:]):
 vl,vr=A(l),A(r);d=sum(d for el,er,d in eta if el<(l+r)/2<er);L+=d*(r-l)*(vl+vr)/2
 sub=[l,r]
 if min(vl,vr)<1<max(vl,vr):sub.insert(1,l+(r-l)*(1-vl)/(vr-vl))
 for left,right in zip(sub,sub[1:]):
  def B(x):return max(F(0),vl+(vr-vl)*(x-l)/(r-l)-1)
  returned+=d*(right-left)*(B(left)+B(right))/2
M=sum(w*A(y) for y,w in atoms);E=source-returned;X=alpha*volume
assert E==F(s['E1_exact']) and E/X==F(s['E1_over_X_exact']) and M-L==F(s['S_exact'])
assert E>0
out=dict(status='passed',original_band_mask_independently_reconstructed=True,actual_J_K_independently_reconstructed=True,stopping_law_reexecuted=False,original_D=s['actual_D'],audit_D=dd,stop_mass_mean_cap_verified=True,bridge_support_pieces=pieces,minimum_u_over_t=str(min(minima)),M=str(M),S=str(M-L),X=str(X),E1=str(E),source_excess=str(source),return_excess=str(returned),E1_over_X=str(E/X),strict_margin_E1_minus_X_over_16=str(E-X/16),finite_D_captures_smoothed_full_E=True,source_file_sha256=hashlib.sha256((D/(name+'.json')).read_bytes()).hexdigest(),main_theorem_proved=False,independent_band_cells=[{'lo':str(l),'hi':str(r),'J':J,'K':K} for l,r,J,K in cells])
(D/'stage40_independent_band_counterexample.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='independent_band_cells'},indent=2))
