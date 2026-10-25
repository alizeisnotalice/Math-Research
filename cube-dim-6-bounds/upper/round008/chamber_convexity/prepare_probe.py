#!/usr/bin/env python3
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
MAIN=HERE.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 d=1000000007;den=100*d;center=[d,853271829,691287433,517389127,307281763,117361927]
 center=[100*v for v in center];points=[]
 for sign in [-1,0,1]:
  x=center.copy();x[5]+=sign*d
  assert x[0]==den and all(x[j]>x[j+1]for j in range(5))and x[-1]>0
  points.append(x)
 assert all(points[0][i]+points[2][i]==2*points[1][i]for i in range(6))
 (HERE/'points.txt').write_text(f'{den} 3\n'+'\n'.join(' '.join(map(str,x))for x in points)+'\n')
 inputs=[MAIN/'upper/round003/direct_kernel/direct_kernel.cpp',MAIN/'upper/round003/direct_kernel/fan_integer_input.txt',
  MAIN/'upper/round003/direct_kernel/exact_generic_reference.py',HERE/'registration.md']
 out={'center_denominator':den,'center_numerators':center,'direction':[0,0,0,0,0,1],
  'h':'1/100','points_order':['minus','center','plus'],'new_G_points':3,'strict_chamber':True,
  'exact_midpoint_relation':True,'points_sha256':sha(HERE/'points.txt'),
  'input_sha256':{str(p):sha(p)for p in inputs},'source_sha256':sha(Path(__file__))}
 (HERE/'input_manifest.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
