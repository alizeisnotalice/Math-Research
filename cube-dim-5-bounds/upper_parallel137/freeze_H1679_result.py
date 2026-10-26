"""Freeze1679 actual new jets and finite overlay; no mathematical rerun."""
from pathlib import Path
import gzip,hashlib,json,os,time
from common import sha,write
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/next1679_H_preflight';B=W/'future_batch0001';P=H/'H1679_packet'
def main():
 began=time.monotonic();dest=H/'H1679_results_manifest.json';assert not dest.exists() and not P.exists();P.mkdir();prepared=json.load(open(H/'next1679_H_prepared_manifest.json'));author=json.load(open(B/'author_validation.json'));assert author['status']=='PASS_all_events_points_sourcecache_qbinding_fulltree_and_exit' and author['new_H_count']==1679;payload=[]
 for name in ['state.json','point_manifest.json','journal.jsonl']:
  src=B/name;out=P/(name+'.gz')
  with src.open('rb') as f,out.open('xb') as g:
   with gzip.GzipFile(fileobj=g,mode='wb',filename='',mtime=0,compresslevel=9) as z:
    while b:=f.read(1024*1024):z.write(b)
  h=hashlib.sha256();n=0
  with gzip.open(out,'rb') as f:
   while b:=f.read(1024*1024):h.update(b);n+=len(b)
  assert h.hexdigest()==sha(src) and n==src.stat().st_size;payload.append({'path':str(out),'sha256':sha(out),'bytes':out.stat().st_size,'original_path':str(src),'original_sha256':sha(src),'decompressed_sha256':h.hexdigest(),'original_bytes':n})
 write(P/'compression_receipt.json',{'status':'PASS_actual_decompressed_original_SHA_equality','pid':os.getpid(),'files':payload,'new_H_gradient_log_corner_calls':0});auth=S/'round017-authorizations/primaryH1679.json';assert sha(auth)=='ad072f0eaf0025579f47cc1303e91d2fe199d6387cc01f905d13895d4dbfdb12';locks={**prepared['source_locks'],str(H/'next1679_H_prepared_manifest.json'):sha(H/'next1679_H_prepared_manifest.json'),str(auth):sha(auth)};external={str(p):sha(p) for p in B.iterdir() if p.is_file()};repo={n:sha(H/n) for n in ['next1679_H_prepared_manifest.json','verify_H1679.py','verify_H1679_vs_previous.diff','H1679_result_report.md','freeze_H1679_result.py']};repo.update({str(p.relative_to(H)):sha(p) for p in P.iterdir()});write(dest,{'schema':'round017_H1679_finite_result_v1','status':'complete_finite_partial_g110_candidate','repository_files':repo,'source_locks':locks,'legacy_code_locks':prepared['legacy_code_locks'],'result_external_files':external,'compressed_payload':payload,'actual_new_H128':1679,'actual_new_physical_gradients':6716,'precision_bits':256,'segments_per_point':1,'new_passed':1562,'new_splits':117,'passed':110967,'pending':234,'splits':58841,'roots':52360,'cache64':184796,'cache128':16480,'faults':0,'obstructions':0,'partial':0,'inflight':0,'cache_hits':0,'actual_workers':15,'all16_production_PIDs_OS_absent':True,'production_wall_seconds':18.814896750001935,'independent_audit_dependency':'143 new1679 numeric/certificate;141 prior671 plus primary acceptance','globalg110_qualified':False,'S_credit':0,'new_child_dispatches':0,'automatic_further_dispatch':False});print(json.dumps({'manifest':str(dest),'sha256':sha(dest),'packet_bytes':sum(p.stat().st_size for p in P.iterdir()),'point_manifest_sha256':sha(B/'point_manifest.json'),'state_sha256':sha(B/'state.json'),'author_validation_sha256':sha(B/'author_validation.json'),'pid':os.getpid(),'wall_seconds':time.monotonic()-began}),flush=True)
if __name__=='__main__':main()
