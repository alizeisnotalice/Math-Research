"""Frozen deterministic rational R037 cube families; no search or batch entry point."""
import hashlib,random
from fractions import Fraction as Q
from math import prod

STRATA=((8,4,48,6),(16,8,64,8),(32,16,80,8),(64,32,96,12),(32,8,96,12))
TEMPLATES=('A','B','C')

def seed(template,stratum,j,domain):
 raw=f'R038-JP-v1|{template}|{stratum}|{j}|{domain}'.encode();return int.from_bytes(hashlib.sha256(raw).digest(),'big')
def generate(template,si,j):
 n,a,N,levels=STRATA[si];delta=(Q(1,64),Q(1,16),Q(1,8),Q(1,4))[j%4];D=a*(Q(1,4),Q(1,2),Q(1),Q(2))[j//4];rng=random.Random(seed(template,si,j,'geometry'))
 if template=='A':
  cubes=[];slots=[]
  for ell in range(levels):
   r=Q(1,4**ell)
   for ux,uy in ((1,0),(-1,0),(0,1),(0,-1)):
    zx,zy=D*ux,D*uy;cubes.append(((zx+2*r*ux,zy+2*r*uy),r))
    for sign in (-1,1):slots.append(((zx+(2-a+2*delta)*r*ux+sign*delta*r*(-uy)/2,zy+(2-a+2*delta)*r*uy+sign*delta*r*ux/2),delta*r/a))
  rng.shuffle(slots);return cubes+slots[:N-len(cubes)],0
 if template=='B':
  cubes=[];slots=[]
  for ell in range(levels):
   r=Q(1,4**ell)
   for sign in (-1,1):
    for uy in (-1,1):cubes.append(((sign*D,2*r*uy),r))
   for tv in (Q(1,4),Q(1,2),Q(3,4)):
    for sign in (-1,1):slots.append((((1-2*tv)*D,sign*(2+delta/4)*r),r))
  rng.shuffle(slots);return cubes+slots[:N-len(cubes)],0
 for retry in range(32):
  rng=random.Random(seed(template,si,j,'retry-'+str(retry)));base=[[rng.choice((-1,1)) for _ in range(n)] for _ in range(2)];seen=set();cubes=[];ok=True
  for i in range(N):
   ell=i%levels;group=(i//levels)%2;r=Q(1,4**ell)
   signs=tuple(base[group][k]*(-1 if rng.randrange(8)==0 else 1) for k in range(n))
   if (ell,group,signs) in seen:ok=False;break
   seen.add((ell,group,signs));numer=(1,n//4,n//2,3*n//4)[j//4];den=n if j//4==0 else n
   coords=tuple((D if group and k==0 else 0)+a*r*signs[k]*(1+delta if rng.randrange(den)<numer else 1-delta) for k in range(n))
   for z,s in cubes:
    if all(abs(coords[k]-z[k])<r+s for k in range(n)):ok=False;break
   if not ok:break
   cubes.append((coords,r))
  if ok:return cubes,retry
 return None,32

def witness_graph(cubes,a):
 n=len(cubes);adj=[0]*n
 for i in range(n):
  for j in range(i):
   if all(abs(x-y)<a*(cubes[i][1]+cubes[j][1]) for x,y in zip(cubes[i][0],cubes[j][0])):
    adj[i]|=1<<j;adj[j]|=1<<i
 return adj

def volumes(n,a,cubes):return [(2*a*r)**n for _,r in cubes]

def geometry_fingerprint(cubes,n):
 """Order-independent exact rational cube-multiset fingerprint; duplicate cubes retained."""
 rows=sorted((tuple(str(v) for v in z)+(('0',)*(n-len(z))),str(r)) for z,r in cubes)
 payload=repr((n,rows)).encode('utf-8')
 return hashlib.sha256(payload).hexdigest()
