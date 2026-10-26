"""Project necessary newq/dual-margin data; deterministic gzip, no math rerun."""
from pathlib import Path
import gzip,hashlib,json,os,time
from fractions import Fraction as F
from util155 import sha,digest,atomic
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_retarget155');O=W/'run001';P=H/'packet'
class Stream:
 def __init__(self,name):self.name=name;self.i=0;self.parts=[];self.open()
 def open(self):self.i+=1;self.path=P/f'{self.name}-{self.i:04d}.jsonl.gz';self.f=open(self.path,'xb');self.z=gzip.GzipFile(fileobj=self.f,mode='wb',filename='',mtime=0,compresslevel=6);self.h=hashlib.sha256();self.n=self.rows=0
 def row(self,x):
  b=(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode();self.z.write(b);self.h.update(b);self.n+=len(b);self.rows+=1
  if self.f.tell()>=85000000:self.close();self.open()
 def close(self):
  if self.f.closed:return
  self.z.close();self.f.close();h=hashlib.sha256();n=rows=0
  with gzip.open(self.path,'rb') as f:
   for b in f:h.update(b);n+=len(b);rows+=1
  assert h.hexdigest()==self.h.hexdigest() and n==self.n and rows==self.rows and self.path.stat().st_size<100000000;self.parts.append({'path':str(self.path.relative_to(H)),'sha256':sha(self.path),'bytes':self.path.stat().st_size,'uncompressed_sha256':h.hexdigest(),'uncompressed_bytes':n,'rows':rows,'decode_SHA_checked':True})
def main():
 began=time.monotonic();assert not P.exists() and not (H/'result_manifest.json').exists();P.mkdir();author=json.load(open(O/'author_validation.json'));assert author['status']=='PASS_all_source_qrecord_leafrecord_hash_prefix_algebra_and_exit_gates';q=Stream('new_q');leaves=Stream('leaf_margins');qn=ln=below128=0;min128=None
 for i in range(15):
  for e in map(json.loads,open(O/f'shard{i:02d}/q_updates.jsonl')):
   r=e['body'];q.row({'canonical_id':r['canonical_id'],'q_new':r['q_new'],'q_new_rule':r['q_new_rule'],'source_event_hash':e['hash'],'source_body_sha256':digest(r)});qn+=1
   if r['epsilon']=='1/128':v=F(r['q_new']);below128+=v<256;min128=v if min128 is None else min(min128,v)
  for e in map(json.loads,open(O/f'shard{i:02d}/leaf_results.jsonl')):
   r=e['body'];leaves.row({**{k:r[k] for k in ('leaf_index','node','canonical_id','margin110_new','margin115','worst_corner_bits','passed115')},'source_event_hash':e['hash'],'source_body_sha256':digest(r)});ln+=1
 q.close();leaves.close();assert qn==107530 and ln==111201
 for name in ['common_log_certificate.json','full_domain_candidate.json','target115_pending.json','actual_attempt0001_summary.json','author_validation.json']:
  x=json.load(open(O/name));atomic(P/name,x,True)
 docs=[];m=json.load(open(H/'prepared_manifest.json'))
 for p in m['source_locks']:
  if '/upper_curvature153/' in p or '/audit_curvature154/' in p:docs.append({'original_path':p,'sha256':sha(p),'utf8':Path(p).read_bytes().decode('utf-8')})
 atomic(P/'theory_sources.json',docs,True);auth=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/round018-authorizations/primary155.json');assert sha(auth)=='e1a0007f2c2e4a7f80a255223412d1cc532a8cb69a1db5fff8f6c3a90b703573';atomic(P/'projection_and_decode_receipt.json',{'pid':os.getpid(),'new_q_rows':qn,'leaf_rows':ln,'q_parts':q.parts,'leaf_parts':leaves.parts,'all_part_files_under100MB':True,'actual_decode_SHA_all_passed':True,'current_q_formula_log_corner_H_reexecution':0,'EPS128_q_below256_objects':below128,'minimum_new_EPS128_q':str(min128),'source151_dependency':'upper_global151/manifest.json SHA2ad3625ea95d29447850e995465afc2480d4b5d61574180af23a080a4b0c98a8;oldcoeff/jets/mapping not duplicated','Delta115_exact_reconstruction':'source151 Delta110 minus packet commonlog upper;raw fullqparams retained immutable externally','wall_seconds':time.monotonic()-began},True)
 names=['prepared_manifest.json','author155.py','archive155.py','result155.md'];repo={n:sha(H/n) for n in names};repo.update({str(p.relative_to(H)):sha(p) for p in P.iterdir()});external={str(p):sha(p) for p in O.rglob('*') if p.is_file()};manifest={'schema':'155_joint153_g115_complete_fixed_candidate_v1','repository_files':repo,'read_only_source_locks':m['source_locks'],'actual_authorization':{'path':str(auth),'sha256':sha(auth)},'math_outputs_manifest_sha256':sha(O/'math_outputs_manifest.json'),'raw_external_math_and_prefix_files':external,'new_q_parts':q.parts,'new_leaf_margin_parts':leaves.parts,'new_q_objects':qn,'leaf16corner_cases':ln,'target115_passed':109961,'target115_pending':1240,'complete_old110_domain_unchanged':True,'complete_new110_margin_domain':True,'candidate_full_floor':json.load(open(O/'full_domain_candidate.json'))['candidate_floor'],'new_log_objects':1,'new_H_gradient_or_split_calls':0,'oldq_corner_baseline_runs':0,'formal_g_or_C5_updated':False,'archive_current_math_reexecutions':0,'all_files_under100MB':all(p.stat().st_size<100000000 for p in P.iterdir()),'source151_oldnumeric_and_coeff_not_copied':True};atomic(H/'result_manifest.json',manifest,True);print(json.dumps({'result_manifest_sha256':sha(H/'result_manifest.json'),'packet_bytes':sum(p.stat().st_size for p in P.iterdir()),'maximum_packet_file_bytes':max(p.stat().st_size for p in P.iterdir()),'new_q_rows':qn,'leaf_rows':ln,'EPS128_newq_below256_objects':below128,'pid':os.getpid(),'wall_seconds':time.monotonic()-began}),flush=True)
if __name__=='__main__':main()
