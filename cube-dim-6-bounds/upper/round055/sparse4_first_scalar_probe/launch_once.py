"""Future one-native launch only after main's separate B055-accepted release."""
import sys
sys.dont_write_bytecode=True
import json,hashlib,subprocess,time,resource
from pathlib import Path
from fractions import Fraction as Q
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bounds(v):
 e=Q(2)**int(v['common_binary_exponent']);return int(v['lower_mantissa'])*e,int(v['upper_mantissa'])*e
def validate_result(n):
 assert n['status']=='COMPLETE_ONE_JOINT_SCALAR_PROBE' and n['precision_bits']==512
 assert n['paired_image_cone_terms']==n['G_DD5_evaluations']==n['H_DD5_evaluations']==25399680
 assert n['G_logical_terms']==n['H_logical_terms']==609592320 and n['total_DD5_evaluations']==50799360
 G,H,K=map(lambda k:bounds(n[k]),['G','H_T','K_T']);assert all(l<=u for l,u in [G,H,K])
 assert max(K[0],G[0]-H[1])<=min(K[1],G[1]-H[0])
 threshold=Q(32,825)/4096;assert threshold==Q(1,105600)
 outcome='UNRESOLVED_NO_AUTOMATIC_RERUN'
 if K[1]<threshold:outcome='EXCLUDES_ONLY_THIS_FIXED_FRAME_TIME_NECESSARY_CONDITION'
 elif K[0]>threshold:outcome='THIS_POINT_NECESSARY_CONDITION_ONLY_PASSED'
 return dict(status='COMPLETE_ONE_FIXED_POINT_COMPARISON',classification=outcome,threshold_exact=str(threshold),zero_node_logical_terms_each_field=n['zero_node_paired_terms']*24,confluent_logical_terms_each_field=n['confluent_paired_terms']*24,global_mT_certified=False,global_g_or_I_certified=False,new_Cn_bound=False)
def main():
 if (P/'full_execution.json').exists():raise SystemExit('Refuse second native attempt')
 release=json.loads((P/'full_release.json').read_text())
 assert release['status']=='MAIN_RELEASED_ONE_NATIVE_AFTER_B055' and release['B055_review_accepted_by_main']
 assert sha(P/'precheck_manifest.json')==release['precheck_manifest_sha256']
 assert sha(Path(release['B055_manifest_path']))==release['B055_manifest_sha256']
 frozen=json.loads((P/'precheck_manifest.json').read_text())
 for f,h in frozen['whitelist_sha256'].items():assert sha(P/f)==h,('frozen local',f)
 for f,h in frozen['external_readonly_sha256'].items():assert sha(Path(f))==h,('frozen external',f)
 s2=json.loads((P/'S2_certificate.json').read_text());assert s2['status']=='PASS_STRICT_OFF_TRUE_S2' and s2['direction_count']==152 and s2['pair_count']==11476 and Q(s2['pair_min_squared'])>0
 assert s2['point_numerators']==[50,50,50,45,27,0] and s2['denominator']==50 and not s2['zero_pairs']
 dep=json.loads((P/'dependency_pins.json').read_text());assert dep['status']=='F050_AND_F052_FINAL_PASS_PINNED'
 build=json.loads((P/'build_record.json').read_text());binary=P/'scratch/sparse4_probe';assert sha(binary)==build['binary_sha256'] and sha(P/'sparse4_probe.cpp')==build['source_sha256']
 assert sha(P/'native_input.txt')==json.loads((P/'input_certificate.json').read_text())['input_sha256']
 cp=P/'full_checkpoint.jsonl';token=P/'full_token.txt';assert not cp.exists() and not token.exists()
 token.write_text('C055_MAIN_RELEASED_FINAL_DEPENDENCIES_CHECKED\n')
 cmd=[str(binary),'full',str(P/'native_input.txt'),str(cp),str(token)]
 rec=dict(attempt=1,phase='LAUNCHING',command=cmd,start_unix=time.time(),timeout_seconds=1800,input_sha256=sha(P/'native_input.txt'),binary_sha256=build['binary_sha256'],source_sha256=build['source_sha256'],main_release_sha256=sha(P/'full_release.json'),precheck_manifest_sha256=sha(P/'precheck_manifest.json'))
 (P/'full_execution.json').write_text(json.dumps(rec,indent=2)+'\n');tic=time.monotonic()
 with (P/'native_full.stdout.json').open('w')as out,(P/'native_full.stderr.txt').open('w')as err:
  child=subprocess.Popen(cmd,stdout=out,stderr=err,start_new_session=True);rec['native_pid']=child.pid;rec['phase']='RUNNING';(P/'full_execution.json').write_text(json.dumps(rec,indent=2)+'\n')
  try:child.wait(timeout=1800);rec['phase']='EXITED'
  except subprocess.TimeoutExpired:
   rec['phase']='TIMEOUT_STOP_ONLY_OWN_NATIVE';child.terminate()
   try:child.wait(timeout=10)
   except subprocess.TimeoutExpired:child.kill();child.wait();rec['needed_SIGKILL']=True
  rec.update(exit_code=child.returncode,process_wall_seconds=time.monotonic()-tic,maxrss_native_platform_units=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,platform=sys.platform,stdout_sha256=sha(P/'native_full.stdout.json'),stderr_sha256=sha(P/'native_full.stderr.txt'))
 (P/'full_execution.json').write_text(json.dumps(rec,indent=2)+'\n')
 if rec['phase']!='EXITED' or rec['exit_code']!=0 or (P/'native_full.stderr.txt').read_text():raise SystemExit('STOP failed native;preserve files;no retry/no finite values')
 n=json.loads((P/'native_full.stdout.json').read_text());comparison=validate_result(n)
 last=json.loads(cp.read_text().splitlines()[-1]);assert last['complete'] and last['paired_image_cone_terms']==25399680
 (P/'point_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n');print(json.dumps(comparison))
if __name__=='__main__':main()
