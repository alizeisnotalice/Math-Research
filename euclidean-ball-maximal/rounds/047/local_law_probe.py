#!/usr/bin/env python3
"""Exact local record obstruction, including continuous uniform-interval inputs."""
from fractions import Fraction as F
import json
H=F(1,4096);LOW=F(1,4);HIGH=F(1,2)
def rad(k):return F(2)**(2-k)
def atoms(t,e):
 p=31*t/512;m=t/2-p;a=m/3;c=F(1,4)-p/6;z=F(1,4)+p/6
 return [(F(-1,8),2*e),(a,m-e),(c,p),(F(10),1-t/2-e)],z,p

def mass(P,x,r,delta=F(0)):
 if not delta:return sum((w for y,w in P if abs(x-y)<r),F(0))
 return sum((w*max(F(0),min(y+delta,x+r)-max(y-delta,x-r))/(2*delta) for y,w in P),F(0))

def maximal(P,x,delta=F(0)):
 if not delta:
  ds=sorted({abs(x-y) for y,w in P});assert ds[0]>0
  return max(sum((w for y,w in P if abs(x-y)<=d),F(0))/(2*d) for d in ds)
 ds=sorted({abs(x-y+s*delta) for y,w in P for s in (-1,1)});assert ds[0]>0
 return max(mass(P,x,d,delta)/(2*d) for d in ds)

def records(trace):
 r=F(0);out={}
 for k,g in enumerate(trace,1):
  a=max(LOW,r);b=min(HIGH,g)
  if a<b:out[k]=(a,b)
  r=max(r,g)
 return out

def case(t,e,batch):
 assert LOW<=t<=HIGH and 0<e<=H and e<t/64
 P,z,p=atoms(t,e);assert sum(w for y,w in P)==1 and min(w for y,w in P)>0
 # Extreme ball-membership checks certify complete rectangles, not just midpoints.
 for center,last in ((F(0),4),(z,8)):
  for k in range(1,last+1):
   for y,w in P:assert abs(abs(center-y)-rad(k))>2*H
  for q in (center-H,center,center+H):
   assert F(7,5)<maximal(P,q)<F(8,5)
   assert F(413,310)<maximal(P,q,H)<F(59,35)
 for x in (-H,F(0),H):
  for zz in (z-H,z,z+H):
   assert zz-x>rad(4)
   for delta in (F(0),H):
    tx=[mass(P,x,rad(k),delta)/(2*rad(k)) for k in range(1,5)]
    tz=[mass(P,zz,rad(k),delta)/(2*rad(k)) for k in range(1,9)]
    assert tx[0]<=F(1,8) and tz[0]<=F(1,8)
    assert tx[3]==t+2*e and tz[3]==t-2*e
    assert max(tx[:3])==t/2+e and max(tz[:7])==t-2*e and tz[7]==31*t/16
    ix=records(tx);iz=records(tz)
    a,b=ix[4];c,d=iz[8]
    assert (max(a,c),min(b,d))==(max(LOW,t-2*e),min(HIGH,t+2*e))
    if delta:
     joint=sum((w*max(F(0),min(y+delta,x+rad(4),zz+rad(8))-max(y-delta,x-rad(4),zz-rad(8)))/(2*delta) for y,w in P),F(0))
    else:joint=sum((w for y,w in P if abs(y-x)<rad(4) and abs(y-zz)<rad(8)),F(0))
    kernel=joint/(4*rad(4)*rad(8));assert kernel==64*p
 return dict(batch=batch,t=str(t),epsilon=str(e),p=str(p),kernel=str(64*p),density_gap=str(4*e),fixed_rectangle_area=str(4*H*H),atomic_and_continuous_checked=True)

LAWS={
 'uniform':(F(1),F(0),[]),
 'increasing_density':(F(0),F(1),[]),
 'endpoint_atoms':(F(0),F(0),[(LOW,F(1,3)),(HIGH,F(2,3))]),
 'mixed':(F(1,4),F(1,4),[(LOW,F(1,8)),(F(5,16),F(1,8)),(HIGH,F(1,4))])}
def interval_law(law,a,b):
 u,v,at=law;l=max(a,LOW);r=min(b,HIGH);L=HIGH-LOW
 cont=F(0) if l>=r else u*(r-l)/L+v*((r-LOW)**2-(l-LOW)**2)/L**2
 return cont+sum((w for t,w in at if a<=t<b),F(0))
def law_tests():
 rows=[]
 for batch,ms in enumerate(((1024,2048),(4096,16384),(2**20,2**32))):
  for m in ms:
   e=F(1,8*m);gap=4*e
   for name,law in LAWS.items():
    u,v,ats=law;inds={0,m,m//2}
    for t,w in ats:
     v0=(t-LOW)*m/(HIGH-LOW);inds.update((v0.numerator//v0.denominator,min(m,v0.numerator//v0.denominator+1)))
    candidates=[]
    for i in inds:
     t=LOW+F(i,m)*(HIGH-LOW);prob=interval_law(law,t-2*e,t+2*e);candidates.append((prob,t))
    prob,t=max(candidates);assert prob>=F(1,m+1)
    P,z,p=atoms(t,e);fee=64*p*prob
    assert fee/gap**2>=F(31,32)*(2*m)**2/(m+1)
    rows.append(dict(batch=batch,law=name,m=m,t=str(t),epsilon=str(e),probability=str(prob),fee=str(fee),fee_over_gap_squared=str(fee/gap**2),pigeonhole_lower=str(F(31,32)*(2*m)**2/(m+1))))
 return rows

def main():
 rows=[]
 for t in (LOW,F(3,8),HIGH):
  for p in (12,16,20):rows.append(case(t,F(1,2**p),0))
 for i in range(9):rows.append(case(LOW+F(i,32),F(1,2**(12+3*i)),1))
 for i in range(15):rows.append(case(LOW+F((37*i+19)%129,512),F(1,2**(48+3*i)),2))
 receipt=dict(status='passed',case_count=len(rows),batch_counts=[sum(r['batch']==b for r in rows) for b in range(3)],cases=rows,law_tests=law_tests(),
  scope='Exact 1D local certificates and interval-density checks; no general maximal bound',rectangle_halfwidth=str(H),smoothing_radius=str(H),continuous_lower=str(F(413,310)),continuous_upper=str(F(59,35)),main_theorem_proved=False)
 print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
