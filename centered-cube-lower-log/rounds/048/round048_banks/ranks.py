from pathlib import Path
from fractions import Fraction as F
import json,time,itertools,collections
HERE=Path(__file__).resolve().parent

def mul(a,b,d,p):
 z=0
 while b:
  if b&1:z^=a
  b>>=1;a<<=1
  if a&(1<<d):a^=p
 return z

def power(a,e,d,p):
 z=1
 while e:
  if e&1:z=mul(z,a,d,p)
  a=mul(a,a,d,p);e>>=1
 return z

def rank(rows):
 piv={}
 for x in rows:
  while x:
   j=x.bit_length()-1
   if j in piv:x^=piv[j]
   else:piv[j]=x;break
 return len(piv)

def frac(x):return x-(x.numerator//x.denominator)

def main():
 m=json.loads((HERE/'MANIFEST.json').read_text());plan=json.loads((HERE/'FROZEN_PLAN.json').read_text());d=m['d'];p=m['modulus'];K=m['K'];banks=m['banks'];out=[];start=time.perf_counter()
 for ls in plan['scales']:
  L=F(ls);offset=K*L/2;cuts=sorted({F(0),F(1),frac(offset),frac(-offset)})
  for left,right in zip(cuts,cuts[1:]):
   theta=(left+right)/2;shifts=[int((theta+sgn*offset)//1)%K for sgn in [1,-1]];factors=[power(m['generator'],s,d,p) for s in shifts]
   rows=[[[mul(x,g,d,p) for x in b['basis']] for g in factors] for b in banks]
   diagonal=[rank(r[0]+r[1]) for r in rows];pair=[]
   for a,b in itertools.combinations(range(4),2):
    pair.append({'banks':[a,b],'four_endpoint_joint_rank':rank(rows[a][0]+rows[a][1]+rows[b][0]+rows[b][1]),'cross_endpoint_pair_ranks':[[rank(rows[a][i]+rows[b][j]) for j in range(2)] for i in range(2)]})
   out.append({'L':str(L),'theta_interval':[str(left),str(right)],'modK_endpoint_offsets':shifts,'same_bank_two_endpoint_ranks':diagonal,'pair_banks':pair})
 hist=collections.Counter(t['four_endpoint_joint_rank'] for r in out for t in r['pair_banks']);cross=collections.Counter(x for r in out for t in r['pair_banks'] for row in t['cross_endpoint_pair_ranks'] for x in row)
 result={'rank_cases':out,'four_endpoint_joint_rank_histogram':dict(hist),'cross_endpoint_pair_rank_histogram':dict(cross),'seconds':time.perf_counter()-start,'meaning':'Exact F2 ranks by polynomial arithmetic and bit elimination. Phase boundaries have Haar measure zero; the continuous clip model is defined there too.'}
 (HERE/'RANKS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'intervals':len(out),'rank44_histogram':dict(hist),'cross_pair_histogram':dict(cross),'seconds':result['seconds']}))
if __name__=='__main__':main()
