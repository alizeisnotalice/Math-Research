#!/usr/bin/env python3
"""Independent scalar/geometry reconstruction; NO kernel replay."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,os,sys,time
import flint
from flint import arb,fmpz,fmpq,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def A(q):return arb(q.numerator)/q.denominator
def R(n,d):return fmpq(fmpz(n),fmpz(d))
def phi(x):return (12*(x/12).asinh()).exp()
start=time.monotonic()
p=json.loads((HERE/'postprocess_results.json').read_text());curv=json.loads((HERE/'curvature_results.json').read_text());man=json.loads((HERE/'curvature_input_manifest.json').read_text())
for x in p['inputs']+man['inputs']:
 assert sha(x['path'])==x['sha256'],x['path']
assert p['T_count']==4 and p['new_kernel_passes']==0
old=json.loads(Path(p['inputs'][0]['path']).read_text());g=json.loads(Path(p['inputs'][1]['path']).read_text())
G=arb(R(g['Q_over_pi5_numerator'],g['Q_over_pi5_denominator']))/arb.pi()**5
hs=[arb(v['H_arb384'])for v in old['time_nodes']];ts=[F(v['tau'])for v in old['time_nodes']]
masses=[(hs[k]-hs[k+1]).lower()for k in range(7)];assert all(x>0 for x in masses)
cc=[F(v,g['denominator'])for v in g['point_numerators'][1:]]
checks=[]
for row in p['rows']:
 j=row['old_time_index'];T=F(row['T']);assert T==ts[j]
 K=G-hs[j];assert (K-arb(row['K_T_center']['ball'])).contains(0)
 radius=F(row['local_band_radius_m0']);m0=T/4000
 lower=lambda r:sum((masses[k]/phi(A(r/ts[k]))for k in range(j)),arb(0))
 if radius:
  assert lower(radius)>A(m0)
  assert lower(radius+F(1,1000000)).upper()<A(m0)
  bs=[(max(F(0),c-radius),min(F(1),c+radius))for c in cc]
  assert all(bs[i][1]<bs[j][0]or bs[j][1]<bs[i][0]for i in range(5)for j in range(i+1,5))
  vol=F(120)
  for lo,hi in bs:vol*=hi-lo
  assert vol==F(row['exact_covered_face_volume'])<1
  assert all([str(lo),str(hi)]==r for (lo,hi),r in zip(bs,row['clipped_box_free_bounds']))
  assert F(row['exact_uncovered_face_volume'])==1-vol
 else:assert K<A(T/4096)
 checks.append({'T':str(T),'PASS':True,'local_only':True})
ac=json.loads(Path(man['inputs'][0]['path']).read_text());points=ac['points'];RR=[R(x['G_times_pi5_numerator'],x['G_times_pi5_denominator'])for x in points]
JG=arb((RR[0]+RR[2])/2-RR[1])/arb.pi()**5
hh=[arb(x['H_T']['ball'])for x in curv['rows']]
JH=(hh[0]+hh[2])/2-hh[1];JK=JG-JH
assert JK>arb(29219)/10**14 and JK<arb(29221)/10**14
assert (JK-arb(curv['K_Jensen_gap']['ball'])).contains(0)
for i,x in enumerate(curv['rows']):
 assert x['point_numerators']==points[i]['point_numerators'] and x['denominator']==points[i]['denominator']
 assert x['complete_group_terms']==967680 and x['computed_sign_paired_terms']==483840
 assert (arb(RR[i])/arb.pi()**5-hh[i]-arb(x['K_T']['ball'])).contains(0)
assert curv['new_H_passes']==3 and curv['new_G_passes']==0 and curv['exit_code']==0
assert not curv['convexity_proved'] and not curv['new_C6_bound'] and not p['global_m_T_certified']
out={'status':'PASS_POSTPROCESS_AND_THREE_POINT_RECONSTRUCTION','pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,'kernel_replayed':False,'new_kernel_passes':0,'hash_checks':'PASS','old_T_checks':checks,'strict_K_Jensen_interval':['29219/100000000000000','29221/100000000000000'],'J_K':JK.str(100,more=True),'new_H_passes_in_source':3,'new_G_passes_in_source':0,'global_m_T_certified':False,'convexity_proved':False,'new_C6_bound':False,'inputs':[{'path':str(HERE/f),'sha256':sha(HERE/f)}for f in ['postprocess_results.json','curvature_results.json','curvature_input_manifest.json']],'source_sha256':sha(__file__),'python':sys.version,'python_flint':flint.__version__,'arb_bits':ctx.prec,'exit_code':0}
(HERE/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
