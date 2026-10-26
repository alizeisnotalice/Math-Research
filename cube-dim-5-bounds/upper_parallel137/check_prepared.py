"""Input-object equality and synthetic checkpoint tests; no corner call."""
from pathlib import Path
import ast,json,os,tempfile,time
from common import sha,digest,write,CORE_PATH,CORE_SHA
from shard_worker import recover
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_parallel137');S=W.parent/'upper_crossanchor131'
def main():
 began=time.monotonic();plan=json.load(open(W/'shard_plan.json'));pairs=list(map(json.loads,open(S/'candidate_pairs.jsonl')));children=json.load(open(S/'child_inputs.json'));need={p['source_id'] for p in pairs};sources={r['source_id']:r for r in map(json.loads,open(S/'source_pool.jsonl')) if r['source_id'] in need};seen=set();owners={};copies=0
 for spec in plan['shards']:
  assert sha(spec['input_path'])==spec['input_sha256'];data=json.load(open(spec['input_path']));sid=spec['shard_id'];expected=[p for p in pairs if p['child_index']%15==sid];assert data['pairs']==expected
  assert data['children']=={str(i):children[i] for i in sorted({p['child_index'] for p in expected})}
  assert data['sources']=={k:sources[k] for k in sorted({p['source_id'] for p in expected})};copies+=len(data['sources'])
  for p in expected:
   assert p['pair_index'] not in seen;seen.add(p['pair_index']);assert owners.setdefault(p['child_index'],sid)==sid
 assert seen==set(range(89426)) and set(owners)==set(range(5427)) and len(sources)==10549
 assert sha(CORE_PATH)==CORE_SHA
 for file in H.glob('*.py'):ast.parse(file.read_text(),filename=str(file))
 # Test receipts only: these invented margins never represent a mathematical certificate.
 toy=[{'pair_index':i,'child_index':0,'child_node':[0,'toy'],'source_id':str(i),'curvature_rule':'toy'} for i in range(2)];fp={'synthetic':True};cert={'source_id':'0','child_node':[0,'toy'],'curvature_rule':'toy','checked_points':16,'passed':True,'minimum_log_margin':'1'}
 with tempfile.TemporaryDirectory(prefix='synthetic_protocol_',dir=W) as td:
  td=Path(td);rp=td/'rows.jsonl';cp=td/'checkpoint.json';r1={'seq':1,'previous':'0'*64,'global_pair_index':0,'child_index':0,'source_id':'0','kind':'evaluated','certificate':cert};r1['hash']=digest(r1);raw1=(json.dumps(r1,sort_keys=True,separators=(',',':'))+'\n').encode();rp.write_bytes(raw1)
  first=recover(rp,cp,toy,fp,True);write(cp,first)
  r2={'seq':2,'previous':r1['hash'],'global_pair_index':1,'child_index':0,'source_id':'1','kind':'skipped','reason':'child already certified by one full16corner source'};r2['hash']=digest(r2);raw2=(json.dumps(r2,sort_keys=True,separators=(',',':'))+'\n').encode();rp.write_bytes(raw1+raw2)
  both=recover(rp,cp,toy,fp,True);assert both['next_pair']==2 and both['evaluated']==1 and both['skipped']==1 and both['qualified_children']=={'0':cert}
  rejects=0
  for wrong_fp,resume in [({'synthetic':False},True),(fp,False)]:
   try:recover(rp,cp,toy,wrong_fp,resume)
   except AssertionError:rejects+=1
  rp.write_bytes(raw1+raw2[:-1])
  try:recover(rp,cp,toy,fp,True)
  except AssertionError:rejects+=1
  assert rejects==3 and rp.read_bytes()==raw1+raw2[:-1]
 result={'status':'PASS_input_equality_and_synthetic_recovery_only','pid':os.getpid(),'pairs_exactly_once':len(seen),'children_unique_worker':len(owners),'unique_selected_sources':len(sources),'source_copies':copies,'all_original_source_and_child_records_equal':True,'stale_checkpoint_valid_extra_rows_recovered_without_retesting':True,'first_pass_causes_later_same_child_skip':True,'wrong_fingerprint_fresh_overwrite_truncated_tail_rejected':True,'actual_source_corner_calls':0,'new_H_gradient_kernel_log_exp_calls':0,'wall_seconds':time.monotonic()-began};write(W/'prepared_check_receipt.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
