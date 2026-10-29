"""Round046 portable whole-basis, true-oracle and mixed-MGF probe replay."""
from pathlib import Path
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile,time
H=Path(__file__).resolve().parent
SKIP={'seconds','elapsed_seconds','search_and_recheck_seconds','runtime','manifest','python','numpy'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def compare(a,b):
 if isinstance(a,dict):
  assert set(a)-SKIP==set(b)-SKIP,(set(a),set(b))
  for k in set(a)-SKIP:compare(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):compare(x,y)
 elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-10),(a,b)
 else:assert a==b,(a,b)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--mpmath-vendor',type=Path);args=ap.parse_args();assert not args.output.exists();start=time.time();prov=json.loads((H/'provenance.json').read_text());env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');runs=[]
 if args.mpmath_vendor:env['CUBE_MPMATH_VENDOR']=str(args.mpmath_vendor.resolve())
 for name,hh in prov['payload_SHA256'].items():assert sha(H/name)==hh,name
 for name,hh in prov['reference_SHA256'].items():assert sha(H.parent/name)==hh,name
 gf=H/'round046_gf';m=json.loads((gf/'MANIFEST.json').read_text());plan=json.loads((gf/'FROZEN_PLAN.json').read_text())
 for name,hh in m['source_sha256'].items():assert sha(gf/name)==hh,name
 assert sha(gf/'PREFLIGHT.json')==plan['certificate_sha256'] and sha(gf/'gf_check.c')==plan['C_source_sha256']
 with tempfile.TemporaryDirectory(prefix='cube-round046-')as td:
  root=Path(td)
  for folder in ['round046_gf','round046_oracle','round046_independent','round046_decoupling']:shutil.copytree(H/folder,root/folder)
  labels=root/'source_d10.npz';shutil.copyfile(H.parent/'045/accepted/source_d10.npz',labels)
  def run(folder,script,output,arguments=(),expected=True,delete_output=True):
   out=root/folder/output
   if delete_output:out.unlink(missing_ok=True)
   p=subprocess.run([sys.executable,str(root/folder/script),*map(str,arguments)],cwd=root/folder,env=env,capture_output=True,text=True,timeout=180)
   if p.returncode:raise RuntimeError(script+'\n'+p.stdout[-2000:]+'\n'+p.stderr[-6000:])
   data=json.loads(out.read_text())
   if expected:compare(json.loads((H/folder/output).read_text()),data)
   runs.append(dict(script=folder+'/'+script,passed=True));return data
  # Replay the frozen search first: it verifies the original preflight byte hash.
  search=run('round046_gf','search.py','MANIFEST.json',delete_output=False)
  compare(json.loads((gf/'SEARCH_ATTEMPTS.json').read_text()),json.loads((root/'round046_gf/SEARCH_ATTEMPTS.json').read_text()))
  compare(json.loads((gf/'PYTHON_SPOT_CHECKS.json').read_text()),json.loads((root/'round046_gf/PYTHON_SPOT_CHECKS.json').read_text()))
  preflight=run('round046_gf','certify.py','PREFLIGHT.json')
  # Restore frozen identity bytes for source-input hashes in subsequent replays.
  shutil.copyfile(gf/'MANIFEST.json',root/'round046_gf/MANIFEST.json')
  trace=run('round046_gf','trace_check.py','TRACE_CHECK.json')
  full=run('round046_independent','verify_gf.py','GF_VERIFIED.json')
  safe=run('round046_independent','safe_packet.py','SAFE_PACKET.json')
  small=run('round046_oracle','validate_oracle.py','validation_results.json',('--labels',labels))
  large=run('round046_oracle','validate_oracle.py','large_manifest_results.json',('--manifest',root/'round046_gf/MANIFEST.json'))
  probe=run('round046_decoupling','probe.py','PROBE.json')
  mgf=run('round046_decoupling','log_mgf_probe.py','LOG_MGF.json')
  out=dict(status='PASS',payload_hashes=len(prov['payload_SHA256']),reference_hashes=len(prov['reference_SHA256']),runs=runs,
   frozen_search_attempts=search['attempts_used'],same_winning_basis=search['basis']==m['basis'],degree=full['d'],basis_rank=full['b'],products=full['full_product_count'],ratios=full['full_ratio_count'],integer_independent_full_check=True,
   trace_pairing_rank=trace['trace_pairing_rank'],safe_nodes=safe['nodes_checked'],safe_packets=safe['finite_packets'],
   primitive_clip_exact_cases=small['small_npz']['clipped_one_cell_exact_checks'],large_source_no_K_table=large['no_K_table'],large_three_scale_label_calls=large['label_calls_for_three_scale_query'],binary64_peak_points_distinct=large['binary64_distinct_coordinates_for_five_peak_points'],
   threshold_Walsh_b_scan=[r['b']for r in probe['threshold_Walsh_scan']],Gauss_numerical_cases=probe['gauss_checks'],small_log_MGF_cases=len(mgf['actual_small_field_cases']),
   new_large_n_coordinate_MC=False,FP_outward_probability_certificate=False,new_Cn_log_gain_proved=False,
   seconds=time.time()-start,python=sys.version,numpy=__import__('numpy').__version__)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
