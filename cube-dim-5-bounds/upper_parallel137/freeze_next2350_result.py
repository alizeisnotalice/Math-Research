"""Freeze finite reuse evidence, no numerical or certificate evaluation."""
from pathlib import Path
import gzip,hashlib,json,os
from common import sha,write
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/next2350_preflight';O=W/'run001';P=H/'next2350_packet'
def main():
 dest=H/'next2350_results_manifest.json';assert not dest.exists() and not P.exists();P.mkdir();m=json.load(open(H/'next2350_prepared_manifest.json'));a=json.load(open(W/'author_validation.json'));assert a['status']=='PASS_sources_all39510_prefixes_single_source_metadata_fulltree_and_exit';payload=[]
 for src,name in [(O/'master_attempt0001_state.json','state.json.gz'),(W/'remaining1679_gap_order.jsonl','remaining1679.jsonl.gz')]:
  p=P/name
  with src.open('rb') as f,p.open('xb') as g:
   with gzip.GzipFile(fileobj=g,mode='wb',filename='',mtime=0,compresslevel=9) as z:
    while b:=f.read(1024*1024):z.write(b)
  h=hashlib.sha256();length=0
  with gzip.open(p,'rb') as f:
   while b:=f.read(1024*1024):h.update(b);length+=len(b)
  assert h.hexdigest()==sha(src) and length==src.stat().st_size;payload.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'original_path':str(src),'original_sha256':sha(src),'decompressed_sha256':h.hexdigest(),'original_bytes':length})
 write(P/'compression_receipt.json',{'status':'PASS_actual_decompression_equal_original_bytes','pid':os.getpid(),'files':payload,'H_log_corner_calls':0});external={str(p):sha(p) for p in sorted(O.rglob('*')) if p.is_file()};external.update({str(W/n):sha(W/n) for n in ['author_registration.json','author_validation.json','remaining1679_gap_order.jsonl']});locks=dict(m['source_locks']);locks[str(H/'next2350_prepared_manifest.json')]=sha(H/'next2350_prepared_manifest.json');authorization=S/'round016-authorizations/primary_next2350.json';assert sha(authorization)=='771379a5c7d195de5e524a5efbc0fa4ac72ca9a0977e5eeea1096a73bd2fb637';locks[str(authorization)]=sha(authorization);repo={n:sha(H/n) for n in ['next2350_prepared_manifest.json','verify_next2350.py','next2350_result_report.md','freeze_next2350_result.py']};repo.update({str(p.relative_to(H)):sha(p) for p in P.iterdir()});write(dest,{'schema':'next2350_reuse_result_v1','status':'complete_fixed39510_partial_g110_candidate','repository_files':repo,'source_locks':locks,'result_external_files':external,'compressed_payload':payload,'fixed_pairs':39510,'single_source_tests':29727,'exact_corners':475632,'skipped_pairs':9783,'new_passed':671,'passed':109405,'pending':1679,'splits':58724,'roots':52360,'cache64':184796,'cache128':14801,'new_H_gradient_kernel_log_exp_calls':0,'new_splits':0,'faults':0,'workers':15,'production_wall_seconds':2.6572561250068247,'new4274_numeric_dependency':'audit138 plus primary acceptance','independent_integer_dependency':'separate review of fixed39510/new671 proofs','globalg110_qualified':False,'S_credit':0,'automatic_further_dispatch':False});print(json.dumps({'manifest':str(dest),'sha256':sha(dest),'packet_bytes':sum(p.stat().st_size for p in P.iterdir()),'state_sha256':sha(O/'master_attempt0001_state.json'),'author_sha256':sha(W/'author_validation.json'),'remaining1679_sha256':sha(W/'remaining1679_gap_order.jsonl')}),flush=True)
if __name__=='__main__':main()
