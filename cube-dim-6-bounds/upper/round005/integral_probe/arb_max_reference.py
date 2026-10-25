#!/usr/bin/env python3
"""Arb384 direct complete-group H_(1/1000) calibration of one Sobol point."""
from __future__ import annotations
from functools import cmp_to_key
import hashlib, itertools, json, sys
from pathlib import Path
from time import perf_counter
from flint import acb, arb, ctx

HERE=Path(__file__).resolve().parent; FAN=HERE/'../../round003/direct_kernel/fan_integer_input.txt'
def cmpn(a,b): return (a[0]*b[1]>b[0]*a[1])-(a[0]*b[1]<b[0]*a[1])
def eq(a,b): return a[0]*b[1]==b[0]*a[1]
def aval(a): return arb(a[0])/a[1]
def parse():
 t=iter(FAN.read_text().split());n=int(next(t));out=[]
 for _ in range(n):
  d=int(next(t));rs=[]
  for _ in range(6):rs.append((int(next(t)),tuple(int(next(t)) for _ in range(6))))
  out.append((d,rs))
 return out
def jet(x,e,k):
 if k==0:return (x/e).atan()
 i=acb(0,1); z=(i**k/acb(e,x)**k-(-i)**k/acb(e,-x)**k)/(2*k*i)
 return (1 if k&1 else -1)*z.real
def dd5(nodes,e):
 z=sorted(nodes,key=cmp_to_key(cmpn));v=[aval(x) for x in z];a=[[arb(0) for _ in range(6)]for _ in range(6)]
 for j in range(6):a[j][0]=jet(v[j],e,0)
 for k in range(1,6):
  for j in range(k,6):a[j][k]=jet(v[j],e,k) if eq(z[j],z[j-k]) else (a[j][k-1]-a[j-1][k-1])/(v[j]-v[j-k])
 return a[5][5]
def main():
 index=int(sys.argv[1]) if len(sys.argv)>1 else 90;t=iter((HERE/'sobol_4x128.txt').read_text().split());den,n=int(next(t)),int(next(t));xs=[tuple(int(next(t))for _ in range(6))for _ in range(n)];assert 0<=index<n
 ctx.prec=384;e=arb(1)/1000;cones=parse();start=perf_counter();total=arb(0);terms=0;events=0
 for perm in itertools.permutations(range(6)):
  for mask in range(64):
   y=tuple((-1 if(mask>>j)&1 else 1)*xs[index][perm[j]]for j in range(6))
   for det,rs in cones:
    nodes=[];alpha=arb(det)
    for q,r in rs:nodes.append((sum(y[j]*r[j]for j in range(6)),den*q));alpha/=q
    z=sorted(nodes,key=cmp_to_key(cmpn));events+=sum(eq(z[j-1],z[j])for j in range(1,6));total+=alpha*dd5(nodes,e);terms+=1
 result=total/(2*arb.pi())**6
 out={'status':'PASS_ARB384_COMPLETE_GROUP_H_1_OVER_1000','index':index,'point_numerators':list(xs[index]),'denominator':den,'terms':terms,'equal_adjacent_events':events,'H_epsilon_arb_384':result.str(95,more=True),'elapsed_seconds':perf_counter()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'points_sha256':hashlib.sha256((HERE/'sobol_4x128.txt').read_bytes()).hexdigest(),'fan_sha256':hashlib.sha256(FAN.read_bytes()).hexdigest()}
 (HERE/'arb_max_reference_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
