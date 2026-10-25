"""Deterministic exact contact preparation; never integrates union volume."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,itertools
import round011_geometry as g
R=Path(__file__).resolve().parent
TARGET=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds/audit/022_results.json')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(name,v): (R/name).write_text(json.dumps(v,indent=2)+'\n')
def ceil(x):return -((-x.numerator)//x.denominator)
def rd(x): # nearest, ties upward
 return F((x*10**12+F(1,2)).numerator//(x*10**12+F(1,2)).denominator,10**12)
def ivadd(a,b):return a[0]+b[0],a[1]+b[1]
def ivmul(a,b):
 q=[x*y for x in a for y in b];return min(q),max(q)
def speed(m,d):
 l,h=g.radius_interval(m);x=g.WTS[2]*d/24
 return (x/h**2,x/l**2) if x>=0 else (x/l**2,x/h**2)
def main():
 assert sha(R/'round008_certificate_input.json')=='8fa317bdb9d2b4209eb18f8fb1bb0fb6bc95abd6861728113c02abec3ab236bf'
 tj=json.loads(TARGET.read_text());dump('round011d_target.json',{'path':str(TARGET),'sha256':sha(TARGET),'input_sha256':tj['input_sha256'],'ratio_interval':tj['ratio_interval']})
 rows=[]
 for i,key in enumerate(g.SEED['keys']):
  box=g.receiver_box(key)
  if box is None:continue
  lo,hi,n=g.block_counts(key);act=[b for b in range(3) if n[b]];m=sum(n[b]*g.WTS[b] for b in act);lc=[];hc=[];bh=[];bl=[]
  for j in range(3):
   bh.append(max(act,key=lambda b:g.P[j]*(hi[b][j]+g.PH[b][j])));bl.append(min(act,key=lambda b:g.P[j]*(lo[b][j]+g.PH[b][j])))
   lc.append(g.P[j]*(hi[bh[-1]][j]+g.PH[bh[-1]][j]));hc.append(g.P[j]*(lo[bl[-1]][j]+g.PH[bl[-1]][j]))
  bb=[[g.bounds(e) for e in ax] for ax in box]
  rows.append({'id':i,'key':key,'lo':lo,'hi':hi,'n':n,'m':m,'lc':lc,'hc':hc,'bh':bh,'bl':bl,'bounds':bb,'speed':speed(m,n[1]-n[2])})
 constraints=[];a=[F(0)]*3;v=[[F(0)]*3 for _ in range(3)]
 for j in range(3):
  best=None
  for B in rows:
   lb=g.bounds(g.expr(B['lc'][j],-1,B['m']))
   if not(g.P[j]/1000<=lb[0] and lb[1]<=999*g.P[j]/1000):continue
   for A in rows:
    if B['id']==A['id'] or B['bh'][j]==A['bl'][j]:continue
    ua=g.bounds(g.expr(A['hc'][j],1,A['m']))
    if not(g.P[j]/1000<=ua[0] and ua[1]<=999*g.P[j]/1000):continue
    beta=ivadd(B['speed'],A['speed'])
    if beta[0]<=0<=beta[1]:continue
    if any(min(B['bounds'][t][1][0],A['bounds'][t][1][0])-max(B['bounds'][t][0][1],A['bounds'][t][0][1])<g.P[t]/10000 for t in range(3) if t!=j):continue
    gap=(lb[0]-ua[1],lb[1]-ua[0]);score=max(abs(x) for x in gap)/g.P[j];order=(score,B['id'],A['id'])
    if best is None or order<best[0]:best=(order,B,A,beta,gap)
  assert best is not None,j
  order,B,A,beta,gap=best;b=B['bh'][j];c=A['bl'][j];co=[B['lc'][j]-A['hc'][j],g.P[j]*(int(b==1)-int(c==1)),g.P[j]*(int(b==2)-int(c==2))];norm=sum(x*x for x in co);mid=sum(beta)/2;u=[rd(mid*x/norm) for x in co];a[j],v[1][j],v[2][j]=u;dot=sum(x*y for x,y in zip(co,u));res=(dot-beta[1],dot-beta[0]);constraints.append({'axis':j,'lower_key_id':B['id'],'upper_key_id':A['id'],'lower_owner':b,'upper_owner':c,'coefficients':list(map(str,co)),'beta_interval':list(map(str,beta)),'gap_interval':list(map(str,gap)),'normalized_gap_upper':str(order[0]),'velocity':list(map(str,u)),'residual_interval':list(map(str,res))})
 def path(s):return ([g.P[j]*(1+a[j]*s) for j in range(3)],[[g.PH[b][j]+v[b][j]*s for j in range(3)] for b in range(3)],[g.WTS[0],g.WTS[1]+g.WTS[2]*s,g.WTS[2]*(1-s)])
 def guard(h):
  ends=[path(-h),path(h)];piv=[(min(e[0][j] for e in ends),max(e[0][j] for e in ends)) for j in range(3)];phiv=[[(min(e[1][b][j] for e in ends),max(e[1][b][j] for e in ends)) for j in range(3)] for b in range(3)]
  if min(x[0] for x in piv)<=1 or any(x<0 or x>1 for e in ends for b in e[1] for x in b):return False,'period/phase guard',None
  minwinner=None
  # All 639 keys, including 84 currently empty receivers: do not silently
  # drop their raw atomic extrema from the branch guard.
  for kk in g.SEED['keys']:
   klo,khi,kn=g.block_counts(kk);act=[b for b in range(3) if kn[b]]
   for j in range(3):
    for inds,sgn in [(khi,1),(klo,-1)]:
     b=max(act,key=lambda b:sgn*(inds[b][j]+g.PH[b][j]))
     for c in act:
      if b==c:continue
      diff=sgn*(inds[b][j]+g.PH[b][j]-inds[c][j]-g.PH[c][j]);margin=diff-h*abs(v[b][j]-v[c][j]);minwinner=margin if minwinner is None else min(minwinner,margin)
      if margin<=0:return False,'all639 atom winner guard',None
  for rr in rows:
   act=[b for b in range(3) if rr['n'][b]]
   for j in range(3):
    for side,idx in [('hi','bh'),('lo','bl')]:
     b=rr[idx][j]
     for c in act:
      if b==c:continue
      diff=rr[side][b][j]+g.PH[b][j]-rr[side][c][j]-g.PH[c][j];margin=abs(diff)-h*abs(v[b][j]-v[c][j]);minwinner=margin if minwinner is None else min(minwinner,margin)
      if margin<=0:return False,'atom winner guard',None
  byid={r['id']:r for r in rows}
  def interval_box(rr):
   ds=g.WTS[2]*(rr['n'][1]-rr['n'][2]);ml=rr['m']-abs(ds)*h;mh=rr['m']+abs(ds)*h;rl=g.radius_interval(ml)[0];rh=g.radius_interval(mh)[1];out=[]
   for j in range(3):
    b=rr['bh'][j];c=rr['bl'][j];lc=ivmul(piv[j],(rr['hi'][b][j]+phiv[b][j][0],rr['hi'][b][j]+phiv[b][j][1]));hc=ivmul(piv[j],(rr['lo'][c][j]+phiv[c][j][0],rr['lo'][c][j]+phiv[c][j][1]));out.append(((lc[0]-rh,lc[1]-rl),(hc[0]+rl,hc[1]+rh)))
   return out
  contactguards=[]
  for cc in constraints:
   j=cc['axis'];B=interval_box(byid[cc['lower_key_id']]);A=interval_box(byid[cc['upper_key_id']]);marg=min(B[j][0][0],piv[j][0]-B[j][0][1],A[j][1][0],piv[j][0]-A[j][1][1]);tang=[min(A[t][1][0],B[t][1][0])-max(A[t][0][1],B[t][0][1]) for t in range(3) if t!=j]
   if marg<=0 or min(tang)<=0:return False,'contact clipping/tangent guard',None
   contactguards.append({'axis':j,'unclipped_margin_lower':str(marg),'tangent_overlap_lower':list(map(str,tang))})
  return True,'all exact interval guards PASS',{'all_original_raw_winner_keys_checked':len(g.SEED['keys']),'minimum_normalized_atom_winner_margin':str(minwinner),'contacts':contactguards,'period_intervals':[[str(t) for t in x] for x in piv],'phase_intervals':[[[str(t) for t in x] for x in row] for row in phiv]}
 trials=[]
 for k in range(9):
  h=F(1,4*2**k);ok,msg,gd=guard(h);trials.append({'s0':str(h),'pass':ok,'reason':msg})
  if ok:break
 out={'status':msg,'input_sha256':sha(R/'round008_certificate_input.json'),'target_sha256':sha(TARGET),'delta':str(g.WTS[2]),'period_velocity':list(map(str,a)),'phase_velocity':[[str(t) for t in row] for row in v],'constraints':constraints,'guard_trials':trials,'guards':gd,'s0':str(h),'steps':['0',str(h/2),str(-h/2),str(h),str(-h)],'baseline_union_reused_without_query':True,'rank':3,'rank_reason':'three nonzero rows have disjoint axis supports','quantization_denominator':10**12,'weights_constant':str(g.W),'symmetric_positive_weight_factor_lower':str(1-h)}
 dump('round011d_direction.json',out)
 if not ok:return
 # Canonicalize only known existing patterns, then conservatively include every
 # periodic placement which can intersect the receiver cell over the path.
 patterns=set()
 raw=json.loads((R/'round011_new_event_keys.json').read_text());raw=raw.get('keys',raw) if isinstance(raw,dict) else raw
 for kk in list(g.SEED['keys'])+raw:
  lo,hi,n=g.block_counts(kk);anchor=next(b for b in range(3) if n[b]);q=lo[anchor];nl=[];nh=[]
  for b in range(3):
   nl+=([lo[b][j]-q[j] for j in range(3)] if n[b] else [0,0,0]);nh+=([hi[b][j]-q[j] for j in range(3)] if n[b] else [-1,-1,-1])
  patterns.add(tuple(nl+nh))
 ends=[path(-h),path(h)];pmin=[min(e[0][j] for e in ends) for j in range(3)];phmin=[[min(e[1][b][j] for e in ends) for j in range(3)] for b in range(3)];phmax=[[max(e[1][b][j] for e in ends) for j in range(3)] for b in range(3)];family=set(map(tuple,g.SEED['keys']));pruned=0
 for kk in sorted(patterns):
  lo,hi,n=g.block_counts(kk);act=[b for b in range(3) if n[b]];m=sum(n[b]*g.WTS[b] for b in act);mm=m+g.WTS[2]*abs(n[1]-n[2])*h;rh=g.radius_interval(mm)[1];spans=[];ranges=[]
  for j in range(3):
   lc=minlc=max(hi[b][j]+phmin[b][j] for b in act);hc=maxhc=min(lo[b][j]+phmax[b][j] for b in act);spans.append(pmin[j]*max(F(0),lc-hc));lower=lc-rh/pmin[j];upper=hc+rh/pmin[j];ranges.append(range((-upper).numerator//(-upper).denominator+1,ceil(1-lower)))
  if any(s**3>=mm for s in spans):pruned+=1;continue
  for q in itertools.product(*ranges):
   nl=[];nh=[]
   for b in range(3):nl+=([lo[b][j]+q[j] for j in range(3)] if n[b] else [0,0,0]);nh+=([hi[b][j]+q[j] for j in range(3)] if n[b] else [-1,-1,-1])
   family.add(tuple(nl+nh))
 import round006_search as eng
 family.update(eng.CLUSTERS);dump('round011d_family.json',{'keys':[list(k) for k in sorted(family)],'known_canonical_patterns':len(patterns),'uniformly_pruned_patterns':pruned,'old639_preserved':True,'cluster_keys':len(eng.CLUSTERS),'total_keys':len(family),'no_new_centers':True})
 print(json.dumps({'direction_sha256':sha(R/'round011d_direction.json'),'s0':str(h),'velocity':out['period_velocity']+out['phase_velocity'][1]+out['phase_velocity'][2],'contacts':[(c['axis'],c['lower_key_id'],c['upper_key_id'],float(F(c['normalized_gap_upper']))) for c in constraints],'family_keys':len(family),'guard':msg},indent=2))
if __name__=='__main__':main()
