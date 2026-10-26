"""Finite result archive and metadata-correction provenance, no mathematics."""
from pathlib import Path
import gzip,hashlib,json,os
from common import sha,write
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/next234_preflight';O=W/'run001';P=H/'next234_packet'
def main():
 dest=H/'next234_results_manifest.json';assert not dest.exists() and not P.exists();P.mkdir();a=json.load(open(W/'author_validation.json'));assert a['new_passed']==108 and a['tree']['queue_and_inflight']==126;payload=[]
 for src,name in [(O/'master_attempt0001_state.json','state.json.gz'),(W/'remaining126_gap_order.jsonl','remaining126.jsonl.gz')]:
  out=P/name
  with src.open('rb') as f,out.open('xb') as g:
   with gzip.GzipFile(fileobj=g,mode='wb',filename='',mtime=0,compresslevel=9) as z:
    while b:=f.read(1024*1024):z.write(b)
  h=hashlib.sha256()
  with gzip.open(out,'rb') as f:
   while b:=f.read(1024*1024):h.update(b)
  assert h.hexdigest()==sha(src);payload.append({'path':str(out),'sha256':sha(out),'bytes':out.stat().st_size,'original_path':str(src),'original_sha256':sha(src),'decompressed_sha256':h.hexdigest()})
 write(P/'compression_receipt.json',{'status':'PASS_actual_decompression_SHA_equality','pid':os.getpid(),'files':payload});m=json.load(open(H/'next234_prepared_manifest.json'));locks=dict(m['source_locks']);locks[str(H/'next234_prepared_manifest.json')]=sha(H/'next234_prepared_manifest.json');auth=S/'round017-authorizations/primary_next234.json';assert sha(auth)=='d23ed6dcd078b613826e68e50d2c3de1ebee1148f2496d6acf3fab3b7d86eca7';locks[str(auth)]=sha(auth);external={str(p):sha(p) for p in O.rglob('*') if p.is_file()}
 for n in ['author_registration.json','author_validation.json','author_metadata_correction.json','remaining126_gap_order.jsonl']:external[str(W/n)]=sha(W/n)
 external.update({str(p):sha(p) for p in (W/'author_metadata_v1').iterdir()});repo={n:sha(H/n) for n in ['next234_prepared_manifest.json','verify_next234.py','next234_result_report.md','freeze_next234_result.py']};repo.update({str(p.relative_to(H)):sha(p) for p in P.iterdir()});write(dest,{'schema':'next234_reuse_result_v1','status':'complete_fixed3904_partial_g110_candidate','repository_files':repo,'source_locks':locks,'result_external_files':external,'compressed_payload':payload,'fixed_pairs':3904,'single_source_tests':2418,'corners':38688,'skip_receipts':1486,'new_passed':108,'passed':111075,'pending':126,'splits':58841,'roots':52360,'cache64':184796,'cache128':16480,'new_H_gradient_kernel_log_exp_calls':0,'new_splits':0,'faults':0,'workers':15,'production_wall_seconds':1.3661286249989644,'author_metadata_correction':'stale displayed constant671 corrected108 from already-asserted actualcount;originals preserved;no math/production rerun','numeric_dependency':'143 plus prior audits and primary acceptance','globalg110_qualified':False,'S_credit':0,'automatic_further_dispatch':False});print(json.dumps({'manifest_sha256':sha(dest),'state_sha256':sha(O/'master_attempt0001_state.json'),'author_sha256':sha(W/'author_validation.json'),'remaining126_sha256':sha(W/'remaining126_gap_order.jsonl'),'packet_bytes':sum(p.stat().st_size for p in P.iterdir())}),flush=True)
if __name__=='__main__':main()
