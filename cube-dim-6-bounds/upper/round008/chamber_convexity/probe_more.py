#!/usr/bin/env python3
"""Adaptive finite diagnostic only; stops on the first negative Jensen gap."""
from decimal import Decimal,getcontext
from pathlib import Path
import hashlib,json,subprocess
HERE=Path(__file__).resolve().parent
KERNEL=HERE.parents[2]/'upper/round003/direct_kernel'
getcontext().prec=60
def main():
 old=json.loads((HERE/'probe_results.json').read_text());den=old['points'][1]['denominator'];c=old['points'][1]['numerators']
 out={'status':'DIAGNOSTIC_ONLY','arithmetic':'C++ direct kernel; binary64 on this platform','points':old['points'],
  'tests':[],'max_new_points':15,'max_centers':2,'max_direction_definitions':5}
 directions=[[0,0,0,0,0,1],[0,0,0,1,-1,0],[0,0,0,1,0,0],[0,1,-1,0,0,0],[0,0,1,-1,0,0]]
 def test(center_id,center_index,minus_index,plus_index,direction):
  rows=out['points'];gap=(Decimal(rows[minus_index]['G_direct_long_double'])+Decimal(rows[plus_index]['G_direct_long_double']))/2-Decimal(rows[center_index]['G_direct_long_double'])
  item={'center_id':center_id,'point_indices':[minus_index,center_index,plus_index],
   'direction':direction,'h':'1/100','Jensen_gap_diagnostic':str(gap)}
  out['tests'].append(item)
  if gap<0:out['first_negative']=item
  print('test',len(out['tests']),'gap',gap,flush=True)
  return gap<0
 def evaluate(points):
  p=HERE/'adaptive_points.txt';p.write_text(f'{den} {len(points)}\n'+'\n'.join(' '.join(map(str,x))for x in points)+'\n')
  proc=subprocess.run([str(HERE/'direct_probe'),str(len(points)),str(p)],cwd=KERNEL,text=True,capture_output=True,check=True)
  rows=json.loads(proc.stdout)['points']
  for row in rows:
   assert row['status']=='GENERIC_OK';row['index']=len(out['points']);out['points'].append(row)
  return list(range(len(out['points'])-len(rows),len(out['points'])))
 negative=test(0,1,0,2,directions[0])
 for direction in directions[1:]:
  if negative:break
  points=[[v+sign*(den//100)*a for v,a in zip(c,direction)]for sign in [-1,1]]
  for x in points:assert x[0]==den and all(x[j]>x[j+1]for j in range(5))and x[-1]>0
  ids=evaluate(points);negative=test(0,1,ids[0],ids[1],direction)
 if not negative:
  center=[den,95137181900,79128743300,59138912700,39128176300,19136192700]
  direction=directions[1];points=[]
  for sign in [-1,0,1]:
   x=[v+sign*(den//100)*a for v,a in zip(center,direction)]
   assert x[0]==den and all(x[j]>x[j+1]for j in range(5))and x[-1]>0;points.append(x)
  ids=evaluate(points);negative=test(1,ids[1],ids[0],ids[2],direction)
 out['new_points_computed']=len(out['points']);assert len(out['points'])<=15
 out['direction_definitions']=directions;out['status']='NEGATIVE_DIAGNOSTIC_STOPPED'if negative else'NO_COUNTEREXAMPLE_FOUND_IN_FINITE_BATCH'
 out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 (HERE/'probe_registry.json').write_text(json.dumps(out,indent=2)+'\n')
 if negative:
  item=out['first_negative'];triple=[out['points'][i]for i in item['point_indices']]
  (HERE/'negative_triple.json').write_text(json.dumps({'test':item,'points':triple},indent=2)+'\n')
 print(out['status'],'points',len(out['points']),flush=True)
if __name__=='__main__':main()
