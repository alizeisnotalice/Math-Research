"""Exact finite safe-packet check via trace-row rank, without discrete logs or K table."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
from verify_gf import cmul,cpower
H=Path(__file__).resolve().parent

def rank(rows):
 p={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k in p:x^=p[k]
   else:p[k]=x;break
 return len(p)
def main():
 target=H/'SAFE_PACKET.json';assert not target.exists();path=H.parent/'round046_gf/MANIFEST.json';m=json.loads(path.read_text());K=m['K'];p=F(99,100);n=10**32;v=F(9,29)*p*p*n/(K*K);basis=m['basis'];poly=m['modulus'];g=m['generator'];N=64;A=F(1,64);T=1;gap=1;P0=2;cache={};nodes=[]
 for i in range(16):
  packet=[]
  for a in range(N+1):
   L=F(10,9)+F(2*i*P0,1)/(K*v)+F(2*a,1)*A/(K*v);z=K*(L-1);floor=z.numerator//z.denominator;ks=[k for k in range(floor-2,floor+3)if abs(z-k)<=2]
   assert F(10,9)<=L<=F(11,9)
   for k in ks:
    kk=k%K
    if kk not in cache:
     x=cpower(g,kk,poly);r=rank(basis+[cmul(x,b,poly)for b in basis]);cache[kk]=dict(k=kk,trace_rank=r,in_defect_set=kk!=0 and r<2*m['b'])
    assert not cache[kk]['in_defect_set']
   packet.append(dict(a=a,L=str(L),nearby_integer_defect_candidates=[k%K for k in ks]))
  nodes.append(dict(i=i,all_nodes_safe=True,first_L=packet[0]['L'],last_L=packet[-1]['L'],nodes_checked=len(packet)))
 out=dict(status='PASS',n=str(n),p=str(p),w=64,v=str(v),T=T,A=str(A),N=N,gap=gap,P0=P0,finite_packets=16,nodes_checked=16*(N+1),exact_rank_checks=list(cache.values()),packets=nodes,manifest_SHA256=hashlib.sha256(path.read_bytes()).hexdigest(),script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Only16 finite packets checked, not the whole n-dependent scale collection, tail events or weak-type constant.')
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],nodes=out['nodes_checked'],distinct_rank_checks=len(cache),ranks=[r['trace_rank']for r in cache.values()]),indent=2))
if __name__=='__main__':main()
