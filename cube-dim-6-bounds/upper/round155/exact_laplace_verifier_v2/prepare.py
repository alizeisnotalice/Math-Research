"""Terra-only v2 preparation: reuse the FIRST v1 packet bytes, never reselect149."""
from pathlib import Path
import hashlib,json,os,sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
V1=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round155/exact_laplace_verifier')
PIN={'case_packet.json':'6587f978228bc7478f25919425d0ed915ff24270786f4a9961027da17b3d495b','source_manifest.json':'3012a7f9a22c838334158228ae2e44c8e2e4f84e10acdaf2fa6039b626e34d9b','registration.json':'c9324fcfac4a4ca1162446c7c1ea2faceff127a2d7961758ace6d90ff6a072df'}
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(v,m):
 if not v:raise ArithmeticError(m)
def write(n,v):
 with(P/n).open('x')as f:json.dump(v,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def main():
 need(not any((P/n).exists()for n in ['case_packet.json','source_manifest.json','registration.json','provenance.json','prepare_receipt.json','parent_release.json','started.json','result.json']),'NO_PREPARE_RETRY_OR_OVERWRITE')
 for n,h in PIN.items():need(sha(V1/n)==h,'ORIGINAL_V1_PREPARE_EVIDENCE')
 need(not any((V1/n).exists()for n in ['started.json','result.json','exit_receipt.json','parent_release.json']),'UNEXPECTED_V1_SCIENCE_EVIDENCE_REQUIRES_MAIN')
 packet=read(V1/'case_packet.json');manifest=read(V1/'source_manifest.json')
 need([c['case_id']for c in packet['selected_cases']]==['case_024','case_015','case_008','case_003','case_000']and packet['cases']==5 and packet['stages']==15,'EXACT_ORIGINAL_FIVE_NO_RESELECTION')
 need(sum(len(c['stages'])for c in packet['selected_cases'])==15,'EXACT_ORIGINAL_STAGE_SET')
 for row in manifest['sources']:need(sha(row['path'])==row['sha256']and Path(row['path']).stat().st_size==row['bytes'],'IMMUTABLE_ORIGINAL_SOURCE')
 for n in ['case_packet.json','source_manifest.json']:
  with(P/n).open('xb')as f:f.write((V1/n).read_bytes());f.flush();os.fsync(f.fileno())
  need(sha(P/n)==PIN[n],'BYTE_IDENTICAL_REUSE')
 provenance={'status':'V2_REUSES_FIRST_V1_SELECTION_BYTE_IDENTICALLY','v1_directory':str(V1),'original_v1_prepare_pins':PIN,'selected_cases_unchanged':True,'reselection_or149_scan':False}
 write('provenance.json',provenance)
 source={n:sha(P/n)for n in ['prepare.py','api_ledger.py','verifier_core.py','run_prepared.py','PROOF.md','OPERATOR.md']}
 write('registration.json',{'status':'V2_PREPARED_UNRUN_PENDING_DUAL_GATE','source_sha256':source,'case_packet_sha256':sha(P/'case_packet.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'provenance_sha256':sha(P/'provenance.json'),'original_v1_prepare_pins':PIN,'cases':5,'stages':15,'workers':1,'nice_minimum':15,'caps':packet['caps'],'counter_budget_semantics':'native_arithmetic caps now include ALL setup/context and construction/arithmetic ATTEMPTS; ALL queries separately counted, not part of3M budget','API_state_semantics':'attempted/completed/raised_exception/inflight_UNKNOWN; process death or interrupted bookkeeping makes exact uncheckpointed totals UNKNOWN','time_limit':None,'automatic_restarts':0,'python_flint_version':'0.9.0','API_runtime_qualification':'unrun until one Terra released batch; no retry','new_density_generations':0,'new_k_integrals':0,'legacy_native_integral_kernel_calls':0})
 write('prepare_receipt.json',{'status':provenance['status'],'case_packet_sha256':sha(P/'case_packet.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'registration_sha256':sha(P/'registration.json'),'new_selection_runs':0,'new_verifier_or_FLINT_calls':0,'new_density_generations':0})
 print(json.dumps({'status':provenance['status'],'registration_sha256':sha(P/'registration.json'),'new_selection_runs':0}))
if __name__=='__main__':main()
