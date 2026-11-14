"""Exact solver-free certificates for a genuine interval catalogue family."""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
import json
HERE=Path(__file__).resolve().parent
def verify(rec):
 M,a,b=map(F,(rec['M'],rec['a'],rec['b']));assert M/4<=b<M/2 and 0<a<M/2-b
 true=M/2;l=a+b;u=M+2*a;pr=[M/u,a/u,a/u]
 lengths=[a,b,M-2*b,b,a];inc=[[0,1,1,1,0],[1,1,0,0,0],[0,0,0,1,1]]
 for mask in range(1,8):
  chosen=[j for j in range(3) if mask>>j&1]
  union=sum((m for z,m in enumerate(lengths) if any(inc[j][z] for j in chosen)),F(0))
  assert sum((pr[j] for j in chosen),F(0))<=union/u
  # Each side alone costs l/2; disjoint sides cost l.
  need=true if 0 in chosen else len(chosen)*l/2
  wl=next(z for z in rec['subcatalogues'] if z['rows']==[j+1 for j in chosen])
  assert F(wl['source'])==need and F(wl['gain'])==union/need
  if mask!=7:assert union/need<u/true
 peak=(b+2*a)/(l*u);resp=1/(2*peak)
 assert resp<=2*b and resp<=true
 # Uniformly fill the two inner M/4 pieces to get the upper certificate.
 w=2*a/u;source_probability=F(1,2)+(M/(4*l)-F(1,2))*w
 upper=true/(2*source_probability)
 assert upper==resp==F(rec['responsible_testing'])
 assert true==F(rec['source'])
 deficit=1-resp/true
 assert deficit==a*(M-2*a-2*b)/(M*(b+2*a))==F(rec['relative_loss'])
 assert 2-4*deficit>0 and (2-4*deficit)**2>=3
 return {'M':str(M),'a':str(a),'b':str(b),'true_half_volume':str(true),'responsible_half_volume':str(resp),'relative_loss':str(deficit),'all_seven_subsets_and_Hall_constraints':True}
if __name__=='__main__':
 d=json.loads((HERE/'wolfram.json').read_text())
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(verify,d['results']))
 assert d['maximum_identity']==0
 out={'status':'verified','workers':12,'cases':len(rows),'method':'Fraction only, no optimizer; checks source, responsible density upper/lower, all subsets, Hall inequalities, rational bound below symbolic family maximum','results':rows}
 (HERE/'replay.json').write_text(json.dumps(out,indent=2)+'\n');print({'status':'verified','cases':len(rows),'workers':12})
