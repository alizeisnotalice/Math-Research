from fractions import Fraction as F
from itertools import product
import json,sys,time
T=time.time(); n1,n2,n3=3,3,2
phi=[[F(0),F(2),F(-1)],[F(1),F(-2),F(3)],[F(-1),F(1),F(0)]]
def vtx(i,j,k):
 i%=n1;j%=n2;k%=n3
 # Each displacement omits its own coordinate; all diagonal derivatives are 1.
 return (F((j-1)*k,3), F((i-1)*(1-k),5), phi[i][j], F((i+2*j)%3-1,7))
def deriv(i,j,k,x,y,z,u):
 ds=[[F(0) for _ in range(3)] for _ in range(3)]
 xyz=(x,y,z)
 for bits in product((0,1),repeat=3):
  V=vtx(i+bits[0],j+bits[1],k+bits[2]); delta=(u**3*V[0],u**3*V[1],u*V[2]+u**3*V[3])
  for d in range(3):
   weight=F(1 if bits[d] else -1)
   for q in range(3):
    if q!=d:weight*=xyz[q] if bits[q] else 1-xyz[q]
   for r in range(3):ds[r][d]+=delta[r]*weight
 for d in range(3):ds[d][d]+=1
 return ds
def fstr(x): return str(x)
def calc(u):
 wt=(F(1,6),F(4,6),F(1,6)); pos=(F(0),F(1,2),F(1))
 results=[]
 for i,j,k in product(range(n1),range(n2),range(n3)):
  integ=[F(0)]*5
  for ia,ib,ic in product(range(3),repeat=3):
   m=deriv(i,j,k,pos[ia],pos[ib],pos[ic],u); q=wt[ia]*wt[ib]*wt[ic]
   xy=m[0][1]*m[1][0]; xz=m[0][2]*m[2][0]; yz=m[1][2]*m[2][1]
   cyc=m[0][1]*m[1][2]*m[2][0]+m[0][2]*m[1][0]*m[2][1]
   for h,zv in enumerate((1-xy-xz-yz+cyc,xy,xz,yz,cyc)): integ[h]+=q*zv
  results.append(integ)
 sums=[sum(x[h] for x in results) for h in range(5)]
 diffs=[x[4] for x in results]; pairvals=[x[2] for x in results]
 return {'u':str(u),'cells':len(results),'mean_volume':fstr(sums[0]/len(results)),'mean_pair':fstr(sum(sums[1:4])/len(results)),'signed_cycle_over_u7':fstr(sums[4]/len(results)/u**7),'max_abs_cell_cycle_over_u7':fstr(max(map(abs,diffs))/u**7),'cell_pair_xz_over_u4_range':[fstr(min(pairvals)/u**4),fstr(max(pairvals)/u**4)],'max_abs_cell_identity_residual':fstr(max(abs(x[0]-(1-sum(x[1:4])+x[4])) for x in results))}
vals=[calc(F(1,2**k)) for k in (4,6,8,10)]
out={'grid':'3x3x2 periodic, trilinear perturbation, tensor Simpson exact','phi':[[str(t) for t in row] for row in phi],'diagonals':'all exactly 1','samples':vals,'seconds':round(time.time()-T,3),'scope':'oriented parameter-cell Jacobian volume; invertibility requires separately bounded u; not an actual halo or packing proof'}
print(json.dumps(out,indent=2))
