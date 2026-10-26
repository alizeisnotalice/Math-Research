#!/usr/bin/env python3
"""Portable exact recursive-DD and independent band/box readback; no new H/log."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,permutations
from collections import Counter
from math import factorial,prod,ceil,floor
import gzip,json,hashlib,time,os,resource,sys
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rq(x):
 x=F(str(x));return fmpq(x.numerator,x.denominator)
def ep(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def orbit(t):return factorial(4)//prod(factorial(v) for v in Counter(t).values())
def coarse(x,scale,up):return F(ceil(x*scale) if up else floor(x*scale),scale)
def recursive_weights(nodes,cones):
 w={}
 for cone in cones:
  ks=sorted(nodes[i] for i in cone['vertices']);tab=[{k*k:fmpq(1,2)} for k in ks]
  for order in range(1,5):
   nexttab=[]
   for i in range(5-order):
    keys=set(tab[i])|set(tab[i+1]);den=ks[i+order]-ks[i]
    nexttab.append({k:(tab[i+1].get(k,fmpq(0))-tab[i].get(k,fmpq(0)))/den for k in keys})
   tab=nexttab
  for k,v in tab[0].items():w[k]=w.get(k,fmpq(0))-cone['det_numerator']*v
 return {k:v for k,v in w.items() if v}
def main():
 began=time.monotonic();ctx.prec=384;reg=json.loads((H/'registration.json').read_text());receipt=json.loads((H/'run_receipt.json').read_text());assert sha(H/'engine.py')==receipt['engine_sha256'];assert sha(H/'registration.json')==receipt['registration_sha256'];assert sha(H/'mechanism.md')==reg['mechanism_sha256'];assert sha(H/'geometry017.json.gz')==reg['sources']['table017']['sha256']
 table=json.loads(gzip.decompress((H/'geometry017.json.gz').read_bytes()));seen=set();summaries=[]
 # N^{-1} tail columns, from independently multiplying the exact source matrix.
 A=[[F(a,6) for a in row] for row in [[1,1,1,1],[-1,-1,-1,2],[-1,-1,2,-1],[-1,2,-1,-1],[2,-1,-1,-1]]]
 Q=[[sum(A[k][i]*A[k][j] for k in range(5)) for j in range(4)] for i in range(4)]
 assert all(Q[i][j]==F(i==j,4)-F(1,36) for i in range(4) for j in range(4))
 for row in reg['cells']:
  r=json.loads((H/(row['id']+'.json')).read_text());assert r['fingerprint']==receipt['fingerprint'];cell=tuple(row['indices']);assert cell not in seen;seen.add(cell);assert cell==tuple(sorted(cell,reverse=True));assert r['weight']==orbit(cell)
  c=list(map(F,r['anchor']));D=r['D'];nodes=[sum(int(x*D)*z for x,z in zip(c,v)) for v in table['vertices']];assert len(set(nodes))==280
  w={int(k):rq(v) for k,v in r['weights'].items()};assert recursive_weights(nodes,table['simplices'])==w;assert len(w)==140 and sum(w.values(),fmpq(0))==0
  assert sum((v*k for k,v in w.items()),fmpq(0))==0
  assert sum((v*k*k for k,v in w.items()),fmpq(0))==-fmpq(sum(s['det_numerator'] for s in table['simplices']),2)
  lo=[F(i,64)-x for i,x in zip(cell,c[1:])];hi=[F(i+1,64)-x for i,x in zip(cell,c[1:])];mean=[(l+h)/2 for l,h in zip(lo,hi)];second=[(l*l+l*h+h*h)/3 for l,h in zip(lo,hi)]
  eta=sum(Q[i][j]*(second[i] if i==j else mean[i]*mean[j]) for i in range(4) for j in range(4));B2=max(sum(Q[i][j]*u[i]*u[j] for i in range(4) for j in range(4)) for u in product(*zip(lo,hi)));h=max(map(abs,lo+hi));lam=max(map(abs,mean));assert F(r['geometry']['eta'])==eta and F(r['geometry']['B2'])==B2 and F(r['geometry']['h_inf'])==h
  rows=r['endpoints'];assert [x['t'] for x in rows]==reg['time_endpoints'] and len(rows)==73
  L=arb(0);U=arb(0);M=arb(0);minmass=None
  for i,x in enumerate(rows):
   a=F(x['t']);nxt=rows[i+1]['H'] if i+1<len(rows) else {'lower':'0','upper':'0'};ml=F(x['H']['lower'])-F(nxt['upper']);mu=F(x['H']['upper'])-F(nxt['lower']);assert 0<ml<=mu;minmass=ml if minmass is None else min(minmass,ml)
   # Independent asinh branch by Chebyshev cosh identity.
   z=arb(rq(1+B2/(10*a*a)));P=16*z**5-20*z**3+5*z;old=P+(P*P-1).sqrt();ex=arb(rq(h/(2*a))).exp()**2;RU=min(ep(old,True),ep(ex,True));C=min(RU,1+lam/a+F(21,8)*eta*RU/(a*a))
   L+=arb(rq(ml/RU));U+=arb(rq(mu*RU));M+=arb(rq(mu*C))
  l=ep(L);u=ep(U,True);m=min(u,ep(M,True));assert 0<l<=m<=u
  q=arb(5)/4;v=arb(rq(F(1,64)**4));lp=arb(rq(l))**q;up=arb(rq(u))**q;chord=v*(lp+(up-lp)*arb(rq((m-l)/(u-l))));rel=chord/(v*arb(rq(rows[0]['H']['lower']))**q)
  b=r['bounds'];cert={'H_lower':str(coarse(F(b['H_lower']),10**16,False)),'H_upper':str(coarse(F(b['H_upper']),10**16,True)),'mean_upper':str(coarse(F(b['mean_upper']),10**16,True)),'chord_upper':str(coarse(F(b['chord_integral_upper']),10**22,True)),'relative_baseline_upper':str(coarse(F(b['relative_center_baseline_upper']),10**8,True))}
  assert F(cert['H_lower'])<=l and u<=F(cert['H_upper']) and m<=F(cert['mean_upper']) and ep(chord,True)<=F(cert['chord_upper']) and ep(rel,True)<=F(cert['relative_baseline_upper'])
  # Exact synchronized-child orbit identity, without a global mesh.
  children={tuple(sorted((2*i+d for i,d in zip(cell,bits)),reverse=True)) for bits in product((0,1),repeat=4)};assert sum(orbit(ch) for ch in children)==16*orbit(cell);assert all(tuple(sorted((i//2 for i in ch),reverse=True))==cell for ch in children)
  summaries.append({'id':r['id'],'class':r['class'],'indices':cell,'H_center_lower':rows[0]['H']['lower'],'certificates':cert,'ratio_display':float(F(b['relative_center_baseline_upper'])),'zero_crossing_ray_count':r['zero_crossing_ray_count'],'taylor_band_count':b['taylor_band_count'],'exp_band_count':b['exp_band_count'],'runtime':r['runtime'],'source_sha256':sha(H/(r['id']+'.json')),'all73_bands_positive':True})
 # Small abstract exact complete-layer checks supplement general combinatorial proof.
 from itertools import combinations_with_replacement
 for m0 in (1,2,3,4,8):assert sum(orbit(c) for c in combinations_with_replacement(range(m0),4))==m0**4
 result={'status':'pass','seconds':time.monotonic()-began,'pid':os.getpid(),'precision_bits':384,'recursive_stencils_checked':32,'positive_bands':32*73,'new_H_or_log_evaluations':0,'new_jet_count':0,'portable_input_table_sha256':sha(H/'geometry017.json.gz'),'readback_source_sha256':sha(__file__),'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'cells':summaries}
 (H/'readback.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='cells'},indent=2))
if __name__=='__main__':main()
