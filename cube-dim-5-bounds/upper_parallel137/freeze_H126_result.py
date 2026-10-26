"""Freeze1679 actual new jets and finite overlay; no mathematical rerun."""
from pathlib import Path
import gzip,hashlib,json,os,time
from common import sha,write
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/next126_H_preflight';B=W/'future_batch0001';P=H/'H126_packet'
def main():
 began=time.monotonic();dest=H/'H126_results_manifest.json';assert not dest.exists() and not P.exists();P.mkdir();prepared=json.load(open(H/'next126_H_prepared_manifest.json'));author=json.load(open(B/'author_validation.json'));assert author['status']=='PASS_all_events_points_sourcecache_qbinding_fulltree_and_exit' and author['new_H_count']==126;payload=[]
 for name in ['state.json','point_manifest.json','journal.jsonl']:
  src=B/name;out=P/(name+'.gz')
  with src.open('rb') as f,out.open('xb') as g:
   with gzip.GzipFile(fileobj=g,mode='wb',filename='',mtime=0,compresslevel=9) as z:
    while b:=f.read(1024*1024):z.write(b)
  h=hashlib.sha256();n=0
  with gzip.open(out,'rb') as f:
   while b:=f.read(1024*1024):h.update(b);n+=len(b)
  assert h.hexdigest()==sha(src) and n==src.stat().st_size;payload.append({'path':str(out),'sha256':sha(out),'bytes':out.stat().st_size,'original_path':str(src),'original_sha256':sha(src),'decompressed_sha256':h.hexdigest(),'original_bytes':n})
 write(P/'compression_receipt.json',{'status':'PASS_actual_decompressed_original_SHA_equality','pid':os.getpid(),'files':payload,'new_H_gradient_log_corner_calls':0});auth=S/'round017-authorizations/primaryH126.json';assert sha(auth)=='592827db259a97f92fa263a52f183bd976190b41c638a42b240c5ae390f249c4';locks={**prepared['source_locks'],str(H/'next126_H_prepared_manifest.json'):sha(H/'next126_H_prepared_manifest.json'),str(auth):sha(auth)};external={str(p):sha(p) for p in B.iterdir() if p.is_file()};repo={n:sha(H/n) for n in ['next126_H_prepared_manifest.json','verify_H126.py','verify_H126_vs_previous.diff','H126_result_report.md','freeze_H126_result.py']};repo.update({str(p.relative_to(H)):sha(p) for p in P.iterdir()});write(dest,{'schema':'round017_H126_finite_result_v1','status':'complete_finite_global_g110_candidate_pending_audits','repository_files':repo,'source_locks':locks,'legacy_code_locks':prepared['legacy_code_locks'],'result_external_files':external,'compressed_payload':payload,'actual_new_H128':126,'actual_new_physical_gradients':504,'precision_bits':256,'segments_per_point':1,'new_passed':126,'new_splits':0,'passed':111201,'pending':0,'splits':58841,'roots':52360,'cache64':184796,'cache128':16606,'faults':0,'obstructions':0,'partial':0,'inflight':0,'cache_hits':0,'actual_workers':15,'all16_production_PIDs_OS_absent':True,'production_wall_seconds':10.118773583002621,'independent_audit_dependency':'new126 numerical/certificate plus all prior sources and primary acceptance','globalg110_qualified':False,'S_credit':0,'new_child_dispatches':0,'automatic_further_dispatch':False});print(json.dumps({'manifest':str(dest),'sha256':sha(dest),'packet_bytes':sum(p.stat().st_size for p in P.iterdir()),'point_manifest_sha256':sha(B/'point_manifest.json'),'state_sha256':sha(B/'state.json'),'author_validation_sha256':sha(B/'author_validation.json'),'pid':os.getpid(),'wall_seconds':time.monotonic()-began}),flush=True)
if __name__=='__main__':main()
