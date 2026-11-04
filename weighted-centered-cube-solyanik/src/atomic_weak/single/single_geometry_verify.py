#!/usr/bin/env python3
"""Independent interval-integral and input-cost check for sigma=2r delta_1."""
import json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def primitive(x,t):return math.copysign(abs(x)**t/t,x) if x else 0.0
def interval_mass(a,b,t):return primitive(b,t)-primitive(a,t)
def direct(t,r):
 # Halo is the single interval (1-r,1+r), which crosses zero exactly once.
 mass=interval_mass(1-r,1+r,t);input_cost=2*r
 return mass/input_cost,(1-r,1+r),mass,input_cost
def formula(t,r):return ((r-1)**t+(r+1)**t)/(2*t*r)
def main():
 rows=[]
 for t in (1/512,1/64,1/8,1/4,1/3):
  for r in (1.001,17/16,1.2,2,4):
   value,I,mass,cost=direct(t,r);f=formula(t,r);rows.append({'tau':t,'r':r,'halo_interval':I,'direct':value,'formula':f,'mass':mass,'input_cost':cost,'error':abs(value-f),'pass':abs(value-f)<2e-12})
 assert all(z['pass'] for z in rows)
 (ROOT/'single_geometry_summary.json').write_text(json.dumps({'rows':rows,'scope':'single crossing interval and cost 2r'},indent=2));print('single_geometry_rows=%d PASS' % len(rows))
if __name__=='__main__':main()
