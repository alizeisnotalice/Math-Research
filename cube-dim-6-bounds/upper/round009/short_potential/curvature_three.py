#!/usr/bin/env python3
"""Three new H_.032 passes at the exact old A008 triple; no G kernel replay."""
from pathlib import Path
from fractions import Fraction as F
from time import perf_counter
import hashlib,importlib.util,json,os,signal,sys
import flint
from flint import arb,fmpz,fmpq,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent; PROJECT=HERE.parents[2]
A_DIR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper/round008/chamber_convexity')
CERT=A_DIR/'exact_certificate.json'; CHECK=A_DIR/'checked_certificate.json'
SOURCE=PROJECT/'upper/round007/g_time_moments/run_moments.py'
FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
T=F(4,125)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(v):return {'ball':v.str(100,more=True),'lower':v.lower().str(100,more=True),'upper':v.upper().str(100,more=True)}
def rational(n,d):return fmpq(fmpz(n),fmpz(d))
start=perf_counter()
old=json.loads(CERT.read_text());checked=json.loads(CHECK.read_text())
assert old['point_order']==['minus','center','plus'] and checked['whole_segment_in_strict_convex_chamber']
pts=old['points'];assert len(pts)==3
rpoints=[[F(v,p['denominator'])for v in p['point_numerators']]for p in pts]
assert all(x[0]==1 and all(x[i]>x[i+1]>0 for i in range(5))for x in rpoints)
assert all(2*rpoints[1][i]==rpoints[0][i]+rpoints[2][i]for i in range(6))
assert rpoints[2][3]-rpoints[1][3]==F(1,100) and rpoints[2][4]-rpoints[1][4]==-F(1,100)
# Targeted duplicate check across the three accessible upper-result trees.
roots=[PROJECT/'upper',Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper'),Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')]
matches=[];scanned=0
for root in roots:
 for p in root.rglob('*.json'):
  if HERE in p.parents:continue
  scanned+=1;txt=p.read_text()
  if any(s in txt for s in ['853271829','50738912693','52738912707']):
   has_H=('H_arb' in txt or 'H_T_arb' in txt)
   matches.append({'path':str(p),'sha256':sha(p),'contains_H_result_field':has_H})
   assert not has_H, 'Possible prior target H record: '+str(p)
(HERE/'duplicate_scan.json').write_text(json.dumps({'status':'NO_PRIOR_H_AT_TARGET_TRIPLE_FOUND','scanned_json_files':scanned,'roots':list(map(str,roots)),'matched_records':matches,'new_H_T':'4/125','new_spatial_points':3},indent=2)+'\n')
assert sha(FAN)=='7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'
Rs=[rational(p['G_times_pi5_numerator'],p['G_times_pi5_denominator'])for p in pts]
J_R=(Rs[0]+Rs[2])/2-Rs[1]
assert J_R==rational(old['exact_Jensen_gap_times_pi5_numerator'],old['exact_Jensen_gap_times_pi5_denominator']) and J_R<0
Gs=[arb(r)/arb.pi()**5 for r in Rs];JG=arb(J_R)/arb.pi()**5
spec=importlib.util.spec_from_file_location('audited_time_kernel_b009',SOURCE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.FAN=FAN
# No anchor or main call: importing the audited evaluator causes no old-directory writes.
original_dd=m.dd;confluent={'count':0}
def counted_dd(z,e):
 if len(set(z))<len(z):confluent['count']+=1
 return original_dd(z,e)
m.dd=counted_dd
inputs=[CERT,CHECK,A_DIR/'negative_triple.json',A_DIR/'report.md',SOURCE,FAN,HERE/'registration.md',HERE/'postprocess_results.json',HERE/'duplicate_scan.json',Path(__file__)]
manifest={'status':'REGISTERED_BEFORE_NEW_KERNEL','T':str(T),'points':[{'label':lab,'numerators':p['point_numerators'],'denominator':p['denominator']}for lab,p in zip(old['point_order'],pts)],'max_new_H_passes':3,'max_wall_seconds':300,'inputs':[{'path':str(p.resolve()),'sha256':sha(p)}for p in inputs]}
(HERE/'curvature_input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
out={'status':'RUNNING_THREE_H_ONLY','pid':os.getpid(),'T':str(T),'python':sys.version,'python_executable':sys.executable,'python_flint':flint.__version__,'arb_bits':384,'new_G_passes':0,'new_H_passes':0,'rows':[],'G_Jensen_gap_reconstructed':ball(JG),'kernel_source_sha256':sha(SOURCE),'source_sha256':sha(__file__),'manifest_sha256':sha(HERE/'curvature_input_manifest.json')}
def write():
 out['elapsed_seconds']=perf_counter()-start
 (HERE/'curvature_results.json').write_text(json.dumps(out,indent=2)+'\n')
def timeout(signum,frame):raise TimeoutError('registered 300-second cap')
signal.signal(signal.SIGALRM,timeout);signal.alarm(300)
try:
 cones=m.fan();write()
 for label,p,G in zip(old['point_order'],pts,Gs):
  confluent['count']=0
  rr=m.evaluate(p['point_numerators'],p['denominator'],arb(T.numerator)/T.denominator,cones,start,gradient=False)
  H=arb(rr['H_arb384']);assert H>0 and G-H>0 and rr['complete_group_terms']==967680
  row={'label':label,'point_numerators':p['point_numerators'],'denominator':p['denominator'],'G_reconstructed_no_kernel':ball(G),'H_T':ball(H),'K_T':ball(G-H),'complete_group_terms':rr['complete_group_terms'],'computed_sign_paired_terms':rr['computed_sign_paired_terms'],'confluent_terms_handled':confluent['count'],'elapsed_seconds':rr['elapsed_seconds']}
  out['rows'].append(row);out['new_H_passes']+=1;write()
  print(json.dumps({'label':label,'H_T':row['H_T']['ball'],'elapsed_seconds':row['elapsed_seconds']},ensure_ascii=False),flush=True)
 hs=[arb(row['H_T']['ball'])for row in out['rows']]
 JH=(hs[0]+hs[2])/2-hs[1];JK=JG-JH
 out.update({'status':'PASS_THREE_POINT_DIAGNOSTIC_STOPPED','H_Jensen_gap':ball(JH),'K_Jensen_gap':ball(JK),'strict_negative_K_Jensen_gap':bool(JK<0),'strict_positive_K_Jensen_gap':bool(JK>0),'chamber_convexity_refuted':bool(JK<0),'convexity_proved':False,'new_C6_bound':False,'global_m_T_certified':False,'exit_code':0,'stopping':'No additional point or time search permitted by this batch.'})
 write();print(json.dumps({'status':out['status'],'pid':out['pid'],'JH':out['H_Jensen_gap'],'JK':out['K_Jensen_gap'],'elapsed_seconds':out['elapsed_seconds']},indent=2),flush=True)
except BaseException as e:
 out.update({'status':'ERROR_STOPPED','error':repr(e),'exit_code':1});write();raise
finally:signal.alarm(0)
