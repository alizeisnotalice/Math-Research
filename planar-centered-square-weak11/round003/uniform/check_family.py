from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json

def union_area(rects):
 xs=sorted({v for a,b,c,d in rects for v in [a,b]});area=F(0)
 for a,b in zip(xs,xs[1:]):
  mid=(a+b)/2;intervals=sorted((c,d) for l,r,c,d in rects if l<mid<r)
  if not intervals:continue
  l,r=intervals[0];length=F(0)
  for c,d in intervals[1:]:
   if c>r:length+=r-l;l,r=c,d
   else:r=max(r,d)
  length+=r-l;area+=(b-a)*length
 return area
out=[];subset_count=0
for N in [1,2,3,8,16,32,64]:
 q=F(1,16*4**N);alpha=1+q;eps=1-1/alpha**2
 a=[F(0)]+[F(2**(i-1),2**N)/alpha for i in range(1,N+1)]
 masses=[4*(a[i]**2-a[i-1]**2) for i in range(1,N+1)]
 assert all(x>0 for x in masses) and sum(masses)+eps==1
 assert 0<eps<=a[1]**2 and q<=eps<=2*q and a[-1]<F(1,2)
 prefix=F(0);atom_rects=[];function_rects=[];margins=[]
 t=[None]+[eps/(16*a[i]) for i in range(1,N+1)];h=t[N]/4
 for i in range(1,N+1):
  prefix+=masses[i-1];assert prefix==4*a[i]**2
  rho=a[i]+t[i];margin=4*a[i]**2+eps-4*rho**2
  assert margin>=F(31,64)*eps>0
  margins.append(margin/eps)
  assert h<=t[i]/4
  assert a[i]+h+t[i]/2<rho and a[i]+h<rho
  assert t[i]/2+h<rho
  atom_rects.extend([(-t[i],F(0),a[i]/2,a[i]),(-t[i],F(0),-a[i],-a[i]/2)])
  function_rects.extend([(-t[i]/2,F(0),a[i]/2,a[i]),(-t[i]/2,F(0),-a[i],-a[i]/2)])
 assert union_area(atom_rects)==N*eps/16
 assert union_area(function_rects)==N*eps/32
 assert h<F(1,128) and a[-1]+h<F(3,4)
 # Exact test of the all-subset old-receiver inequality for modest N only.
 if N<=8:
  for mask in range(1,1<<N):
   ids=[i for i in range(N) if mask>>i&1];j=max(ids)
   assert sum(masses[i] for i in ids)<=4*a[j+1]**2
   subset_count+=1
 out.append({'N':N,'epsilon':str(eps),'total_mass':'1','atom_certificate_area_over_epsilon':str(union_area(atom_rects)/eps),'L1_certificate_area_over_epsilon':str(union_area(function_rects)/eps),'minimum_witness_margin_over_epsilon':str(min(margins)),'old_subset_tested':N<=8,'all_passed':True})
 if N==8:
  sample={'N':N,'epsilon':str(eps),'atoms':[{'position':[str(a[i]),'0'],'mass':str(masses[i-1])} for i in range(1,N+1)]+[{'position':['0','0'],'mass':str(eps)}],'h':str(h),'atom_rectangles':[[str(z) for z in r] for r in atom_rects],'L1_rectangles':[[str(z) for z in r] for r in function_rects]}
  Path(__file__).with_name('sample-N8.json').write_text(json.dumps(sample,indent=2)+'\n')
print(json.dumps({'scope':'Exact input, strict witness, support and rational-certificate area checks. Does not certify all positive radii, arbitrary N, true full halo areas or a global weak constant; those statements require the analytic proof.','cases':out,'old_nonempty_subsets_checked':subset_count,'all_passed':True},indent=2))
