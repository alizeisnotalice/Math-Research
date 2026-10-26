"""New finite-result archive only, no numerical or certificate calculation."""
from pathlib import Path
import gzip,hashlib,json,os,time
from common import sha,write
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/residual_preflight';P=H/'residual_packet'
def main():
 began=time.monotonic();dest=H/'residual_results_manifest.json';assert not dest.exists() and not P.exists();P.mkdir();prepared=json.load(open(H/'residual_prepared_manifest.json'));locks=dict(prepared['source_locks']);locks[str(H/'residual_prepared_manifest.json')]=sha(H/'residual_prepared_manifest.json');results={};compressed=[];receipts=[]
 for batch in (1,2):
  b=W/f'future_batch{batch:04d}';audit=json.load(open(b/'author_validation.json'));assert audit['status']=='PASS_all_events_points_sourcecache_qbinding_fulltree_and_exit';receipts.append(audit)
  for file in sorted(b.iterdir()):
   if file.is_file():results[str(file)]=sha(file)
  for name in ('state.json','point_manifest.json','journal.jsonl'):
   raw=b/name;packed=P/f'batch{batch:04d}_{name}.gz'
   with raw.open('rb') as src,packed.open('xb') as dst:
    with gzip.GzipFile(fileobj=dst,mode='wb',filename='',mtime=0,compresslevel=9) as z:
     while block:=src.read(1024*1024):z.write(block)
   hasher=hashlib.sha256();count=0
   with gzip.open(packed,'rb') as src:
    while block:=src.read(1024*1024):hasher.update(block);count+=len(block)
   assert hasher.hexdigest()==results[str(raw)] and count==raw.stat().st_size;compressed.append({'path':str(packed),'sha256':sha(packed),'bytes':packed.stat().st_size,'original_path':str(raw),'original_sha256':results[str(raw)],'decompressed_sha256':hasher.hexdigest(),'original_bytes':count})
 auth={str(S/'round016-authorizations/primary_residual137_batch1.json'):'0c96da764d5cb29ee32304dbd56b8555619e3dab5555c16d1fc230bd13bf49ec',str(S/'round016-authorizations/primary_residual137_batch2.json'):'5a3af9d28373bf6d10a0d989e17331f28ea1a3313039d2fb4cb43d1b6348f508'}
 for p,h in auth.items():assert sha(p)==h
 locks.update(auth);write(P/'compression_receipt.json',{'status':'PASS_actual_decompression_SHA_equal_originals','pid':os.getpid(),'files':compressed,'new_H_gradient_kernel_log_exp_corner_calls':0,'wall_seconds':time.monotonic()-began});write(P/'source_refs.json',{'source_prepared_manifest':str(H/'residual_prepared_manifest.json'),'source_prepared_manifest_sha256':sha(H/'residual_prepared_manifest.json'),'result137_manifest':str(H/'result_manifest.json'),'result137_manifest_sha256':sha(H/'result_manifest.json'),'source_seed_123':str(S/'upper_continue123/seed.json'),'source_seed_123_sha256':'52fe7f9f4ecd8d4c16c1c4d48423f78bc09821770fb8e3c5d8899a9dfa4f691b','source_locks':locks,'authorizations':auth,'finite_inputs':prepared['finite_batches'],'scope':'new4274 physicalEPS128 jets and finite overlays only;immutable source cache not copied'})
 repo={n:sha(H/n) for n in ['residual_prepared_manifest.json','verify_residual_batch.py','residual_result_report.md','freeze_residual_result.py']};repo.update({str(p.relative_to(H)):sha(p) for p in sorted(P.iterdir())});write(dest,{'schema':'residual137_twofinite_result_v1','status':'complete_original4274_inputs_partial_g110_candidate','repository_files':repo,'source_locks':locks,'legacy_code_locks':prepared['legacy_code_locks'],'result_external_files':results,'compressed_files':compressed,'author_receipts':[{'batch':r['batch'],'path':str(W/f'future_batch{r["batch"]:04d}'/'author_validation.json'),'sha256':sha(W/f'future_batch{r["batch"]:04d}'/'author_validation.json')} for r in receipts],'new_H128':4274,'new_physical_gradients':17096,'new_passed':3099,'new_splits':1175,'passed':108734,'pending':2350,'splits':58724,'cache64':184796,'cache128':14801,'root_count':52360,'faults':0,'obstructions':0,'inflight':0,'actual_workers_per_batch':15,'production_wall_seconds':[39.079243541986216,9.920766291004838],'new_child_dispatches':0,'globalg110_qualified':False,'S_credit':0,'no_automatic_further_dispatch':True});print(json.dumps({'manifest':str(dest),'sha256':sha(dest),'packet_bytes':sum(p.stat().st_size for p in P.iterdir()),'compressed_payload_bytes':sum(r['bytes'] for r in compressed),'compressed_files':len(compressed),'result_external_files':len(results),'pid':os.getpid(),'wall_seconds':time.monotonic()-began}),flush=True)
if __name__=='__main__':main()
