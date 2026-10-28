"""Explicit binary cubic-extension subspaces and two-segment Walsh diagnostics."""
from pathlib import Path
import numpy as np,json,math,hashlib
OUT=Path(__file__).resolve().parent
reg=OUT/'registration.json'
if reg.exists():raise RuntimeError('Refusing repeated registered run')
plan={'b':[2,3,4,5],'primitive_polynomial_candidates':{'6':0x43,'9':0x211,'12':0x1053,'15':0x8003},
 'Walsh_B':[4,8,16,32],'Walsh_theta_denominator_factor':4,
 'scope':'finite exact field counts and floating exact-grid Walsh integrals; no C_n bound',
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
reg.write_text(json.dumps(plan,indent=2));fields=[]
for b in plan['b']:
 B=1<<b;d=3*b;q=1<<d;K=q-1;poly=plan['primitive_polynomial_candidates'][str(d)]
 ex=[];z=1
 for j in range(K):
  ex.append(z);z<<=1
  if z&q:z^=poly
 assert z==1 and len(set(ex))==K
 ex=np.array(ex,dtype=np.int64);log=np.full(q,-1,dtype=np.int64);log[ex]=np.arange(K)
 sub=np.r_[0,ex[np.arange(B-1)*(K//(B-1))]]
 assert set(int(x)^int(y)for x in sub for y in sub)==set(map(int,sub))
 W=np.array([0 if u==0 else int(u)^int(ex[(2*log[u]+1)%K])for u in sub],dtype=np.int64)
 assert len(set(W))==B and set(int(x)^int(y)for x in W for y in W)==set(map(int,W))
 w=W[W>0];logs=log[w]
 products=[int(ex[(logs[i]+logs[j])%K])for i in range(len(w))for j in range(i,len(w))]
 assert len(set(products))==len(products)
 D=np.unique((logs[:,None]-logs[None,:])%K);D=D[D>0]
 assert len(D)==(B-1)*(B-2)
 flags=np.zeros(K,dtype=bool);flags[D]=True
 overlaps=np.bincount(((D[:,None]-D[None,:])%K).ravel(),minlength=K)
 primitive=np.array([s for s in range(1,K)if math.gcd(s,K)==1],dtype=np.int64)
 short=primitive[(~flags[primitive])&(~flags[2*primitive%K])]
 assert len(short)>0
 best=int(short[np.argmin(overlaps[short])]);minoverlap=int(overlaps[best])
 assert minoverlap==len(set(map(int,D)) & set(map(int,(D+best)%K)))
 fields.append({'B':B,'q':q,'K':K,'polynomial_hex':hex(poly),'W':W.tolist(),
   'Sidon_product_count':len(products),'distinct_nontrivial_ratio_count':len(D),
   'primitive_count':len(primitive),'primitive_short_lag_good_count':len(short),
   'minimum_adjacent_defect_overlap_among_short_good_primitive':minoverlap,
   'minimizing_primitive_exponent':best})
walsh=[]
for B in plan['Walsh_B']:
 N=4*B;y=np.arange(B);a=np.arange(1,B)
 chi=np.array([[1-2*((int(aa)&int(yy)).bit_count()%2)for yy in y]for aa in a])
 cumulative=np.c_[np.zeros(B-1),np.cumsum(chi,axis=1)]/B
 r=(np.arange(N)+.5)/N
 index=np.searchsorted((y+.5)/B,r,side='right');F=cumulative[:,index]
 mask=~np.eye(B-1,dtype=bool);maxsum=-1;arg=None
 for k in range(N+1):
  W1=F[:,:N-k]@F[:,k:].T/N if k<N else np.zeros((B-1,B-1))
  W2=F[:,N-k:]@F[:,:k].T/N if k else np.zeros((B-1,B-1))
  A1=np.where(mask,abs(W1),-1);A2=np.where(mask,abs(W2),-1)
  val=float(A1.max()+A2.max())
  if val>maxsum:
   p1=np.unravel_index(A1.argmax(),A1.shape);p2=np.unravel_index(A2.argmax(),A2.shape)
   maxsum=val;arg={'theta_numerator':k,'theta_denominator':N,
    'pair1':[int(a[p1[0]]),int(a[p1[1]])],'pair2':[int(a[p2[0]]),int(a[p2[1]])]}
 walsh.append({'B':B,'max_sum_absolute_disjoint_segment_cross_terms':maxsum,'argmax':arg})
res={'numpy':np.__version__,'registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),
     'fields':fields,'Walsh':walsh}
(OUT/'results.json').write_text(json.dumps(res,indent=2))
for row in fields:print('field',row['B'],'short-good',row['primitive_short_lag_good_count'],'min adjacent overlap',row['minimum_adjacent_defect_overlap_among_short_good_primitive'])
for row in walsh:print('Walsh',row['B'],row['max_sum_absolute_disjoint_segment_cross_terms'],row['argmax'])
