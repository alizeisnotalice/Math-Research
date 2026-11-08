import json,sys
from pathlib import Path
from fractions import Fraction as Q
from columns78 import master
from allocate77 import phi_interval
def make(extended=False):
 v=Q(16,25);rad=Q(1,10**7);sig=1+Q(1,10**8)
 pts=[(Q(0),Q(2,1000)),(Q(1,1000),Q(21,20000)),(Q(2,1000),Q(0))]
 hs=[(Q(21,20),(1-Q(2077,10**6))/Q(21,20)),(Q(1),1-Q(2025,10**6)),((1-Q(2075,10**6))/Q(21,20),Q(21,20))]
 if extended:hs.append((Q(9,10),Q(9,10)))
 vols=[[(a+x+rad)*(b+y+rad) for x,y in pts] for a,b in hs]
 mins=[[(a+x-rad)*(b+y-rad) for x,y in pts] for a,b in hs]
 eligible=[[i for i,V in enumerate(row) if sig*V<=1] for row in vols]
 assert eligible[:3]==[[1,2],[0,2],[0,1]]
 assert all(mins[j][j]>1 for j in range(3))
 if extended:assert eligible[-1]==[0,1,2]
 cols=[]
 for j,ids in enumerate(eligible):
  em=sum(1<<i for i in ids)
  for mask in range(1,8):
   if mask&em!=mask:continue
   mass=sig*max(vols[j][i] for i in ids if mask>>i&1);g=mass/v
   assert 1<g<=Q(25,16)
   lo,hi=phi_interval(g)
   cols.append(dict(anchor_index=j,mask=mask,mass=str(mass),density=str(g),lower=str(v*lo),upper=str(v*hi)))
 c=dict(name='triangle_plus_corner' if extended else 'triangle_three_anchors',receiver_cells=[list(map(str,p)) for p in pts],source_volume=str(v),source_density='25/16',sigma=str(sig),receiver_radius=str(rad),anchor_sides=[[str(z) for z in h] for h in hs],maximum_volumes=[[str(z) for z in row] for row in vols],minimum_volumes=[[str(z) for z in row] for row in mins],eligible=eligible,columns=cols)
 c['cover']=master(c,cols)
 return c
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
 rows=[make(False),make(True)]
 (out/'triangle78.json').write_text(json.dumps(dict(cases=rows),indent=2)+'\n')
 for c in rows:
  d=c['cover'];print(c['name'],'columns',len(c['columns']),'IP/FP certified',d['certified_integer_over_fractional_lower'],d['certified_integer_over_fractional_upper'],'weights',d['weights'])
 proofs=[]
 for c in rows:
  cols=c['columns'];d=c['cover']
  if len(c['anchor_sides'])==3:
   pairs={z['mask']:z for z in cols if z['mask'].bit_count()==2}
   lam=[]
   for i in range(3):
    pos=[z for mask,z in pairs.items() if mask>>i&1];neg=next(z for mask,z in pairs.items() if not(mask>>i&1))
    l=(sum((Q(z['lower']) for z in pos),Q(0))-Q(neg['upper']))/2
    u=(sum((Q(z['upper']) for z in pos),Q(0))-Q(neg['lower']))/2
    assert l>=0 and all(u<=Q(z['lower']) for z in cols if z['mask']==1<<i)
    lam.append([str(l),str(u)])
   fl=sum((Q(z['lower']) for z in pairs.values()),Q(0))/2;fu=sum((Q(z['upper']) for z in pairs.values()),Q(0))/2
   proofs.append(dict(name=c['name'],analytic_fractional_lower=str(fl),analytic_fractional_upper=str(fu),supporting_price_intervals=lam,
     integer_over_fractional_lower=str(Q(d['integer_lower'])/fu),integer_over_fractional_upper=str(Q(d['integer_upper'])/fl)))
  else:
   whole=next(z for z in cols if z['anchor_index']==3 and z['mask']==7)
   pivot=next(z for z in cols if z['anchor_index']==3 and z['mask']==2)
   assert whole['mass']==pivot['mass']
   assert all(Q(z['mass'])>=Q(whole['mass']) for z in cols if z['mask']&2)
   proofs.append(dict(name=c['name'],pivot_cell=1,exact_optimum_lower=whole['lower'],exact_optimum_upper=whole['upper'],integer_equals_fractional_exactly=True))
 (out/'analytic78.json').write_text(json.dumps(dict(cases=proofs),indent=2)+'\n')
 print('PASS analytic half-pair dual and exact natural-corner repair')
