"""Pure exact-field projection, deterministic compression and decode-SHA checks."""
from pathlib import Path
from fractions import Fraction as F
import ast,base64,gzip,hashlib,importlib.util,json,os,sys,time,resource
D=Path('/Users/zhengzhihao/Desktop/T');S=D/'cube-finite-dim-scratch';H=Path(__file__).resolve().parent;R=H.parent/'upper_parallel137';A=S/'upper_retarget149';W=S/'upper_global151';MODEL='n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()
def write(p,x):
 assert not Path(p).exists();Path(p).write_bytes(canonical(x)+b'\n')
class Stream:
 def __init__(self,name):self.name=name;self.parts=[];self.seq=0;self.count=0;self.open()
 def open(self):
  self.seq+=1;self.path=H/'data'/f'{self.name}-{self.seq:04d}.jsonl.gz';self.raw=open(self.path,'xb');self.z=gzip.GzipFile(fileobj=self.raw,mode='wb',filename='',mtime=0,compresslevel=6);self.h=hashlib.sha256();self.n=self.rows=0
 def row(self,obj):
  b=canonical(obj)+b'\n';self.z.write(b);self.h.update(b);self.n+=len(b);self.rows+=1;self.count+=1
  if self.raw.tell()>=85000000:self.close();self.open()
 def close(self):
  if self.raw.closed:return
  self.z.close();self.raw.close();assert self.path.stat().st_size<100000000
  h=hashlib.sha256();n=rows=0
  with gzip.open(self.path,'rb') as f:
   for b in f:h.update(b);n+=len(b);rows+=1
  assert h.hexdigest()==self.h.hexdigest() and n==self.n and rows==self.rows
  self.parts.append({'path':str(self.path.relative_to(H)),'sha256':sha(self.path),'bytes':self.path.stat().st_size,'uncompressed_sha256':self.h.hexdigest(),'uncompressed_bytes':n,'records':rows,'actual_decode_SHA_equal_written_projection':True})
def main():
 began=time.monotonic();os.nice(10);assert not (H/'data').exists() and not (H/'manifest.json').exists();(H/'data').mkdir();assert (W/'layout_registration.json').exists();source=H.parent/'upper_retarget149/prepared_manifest.json';assert sha(source)=='bcbd1f26218bc63018170463fac1b69ffbdad0bda003632249cf678ad0434f32';m=json.load(open(source));locks=dict(m['source_locks']);locks[str(source)]=sha(source)
 for p,h in m['prepared_external_files'].items():assert sha(p)==h
 for p,h in locks.items():assert sha(p)==h
 write(W/'packing_registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'layout_sha256':sha(H/'layout.json'),'registered_before_projection_and_large_copy':True,'source149_manifest_sha256':sha(source),'new_H_gradient_log_corner_calls':0,'checker_execution':False,'nice':os.getpriority(os.PRIO_PROCESS,0)})
 sys.path.insert(0,str(R));spec=importlib.util.spec_from_file_location('readonly151_H126_codec',R/'H126_run.py');codec=importlib.util.module_from_spec(spec);spec.loader.exec_module(codec);seed=codec.source_seed({});lastpath=S/'upper_parallel137/next126_H_preflight/future_batch0001/state.json';last=json.load(open(lastpath))
 for k,v in last['new_cache_by_epsilon']['1/128'].items():assert k not in seed['cache_by_epsilon']['1/128'];seed['cache_by_epsilon']['1/128'][k]=v
 cache=seed['cache_by_epsilon'];assert len(cache['1/64'])==184796 and len(cache['1/128'])==16606
 needed=set();rawshas={};objects=0
 for o in map(json.loads,open(A/'canonical_objects.jsonl')):
  e,key=o['typed_cache_key'];record=cache[e][key];assert digest(record)==o['raw_record_sha256'] and record['values'][0]==o['H_interval'];needed.add((e,key));rawshas[(e,key)]=o['raw_record_sha256'];objects+=1
 assert objects==107530;supportnodes=set()
 for row in map(json.loads,open(A/'source_provenance.jsonl')):
  if row['curvature_proof'] is not None:supportnodes.add(tuple(row['curvature_proof']['node']))
 assert len(supportnodes)==3812;supports={};basepath=S/'upper_mixed098/final_state.json';base=json.load(open(basepath));locks[str(basepath)]=sha(basepath)
 for leaf in base['passed']:
  node=tuple(leaf['node'])
  if node in supportnodes:
   c=leaf['certificate'];target=c.get('target',leaf.get('target'));epsilon=c.get('epsilon',leaf['leaf_epsilon']);assert leaf['leaf_epsilon']=='1/64' and F(target)==F(21,200000) and F(epsilon)==F(1,64) and c['passed'];key=c.get('point_key',c.get('pair_cache_key',[None,None])[1]) or c.get('anchor');assert isinstance(key,str);oldrecord=base['cache_by_epsilon']['1/64'][key];record=cache['1/64'][key];assert digest(oldrecord)==digest(record);needed.add(('1/64',key));rawshas[('1/64',key)]=digest(record)
   if 'penalty_coefficient' in c:penalty=c['penalty_coefficient'];qsource='original penalty_coefficient'
   else:assert c['quadratic']=='Q_full=I4/16' and F(epsilon)==F(1,64);penalty=str(1/(16*F(epsilon)**2));qsource='exact original Q_full=I4/16 andepsilon1/64;scaled coefficient1/(16epsilon²)'
   view={k:c[k] for k in ('Delta_lower','loggradient_directed','minimum_log_margin','checked_points','passed')};view.update({'target':target,'epsilon':'1/64','point_key':key,'penalty_coefficient':penalty})
   supports[node]={'node':leaf['node'],'typed_cache_key':['1/64',key],'raw_record_sha256':digest(record),'source_leaf_sha256':digest(leaf),'source_state_sha256':locks[str(basepath)],'certificate':view,'original_certificate':c,'original_certificate_sha256':digest(c),'field_normalization':{'point_key':'point_key/pair_cache_key/originalanchor,unrounded','penalty':qsource}}
 assert set(supports)==supportnodes;del base
 # Source-stage membership only, not a numerical or certificate re-evaluation.
 origin={};stages={};documents={};auditbase=D/'Math-Research-cube5-audit123-sol/cube-dim-5-bounds';legacyaudit=D/'Math-Research-cube5-lower-audit-cont/cube-dim-5-bounds'
 def document(p):
  p=Path(p);assert p.exists();key=sha(p);locks[str(p)]=key
  if key not in documents:documents[key]={'document_id':key,'original_path':str(p),'original_sha256':key,'encoding':'utf-8','text':p.read_bytes().decode('utf-8')}
  return key
 def stage(name,statepath,epsilon,field,refs):
  statepath=Path(statepath);state=json.load(open(statepath));records=state[field] if field=='new_cache' else state[field][epsilon];sourcehash=sha(statepath);locks[str(statepath)]=sourcehash;stages[name]={'production_state':str(statepath),'production_state_sha256':sourcehash,'cache_field':field,'epsilon':epsilon,'numeric_qualification_documents':[document(p) for p in refs]}
  for key,record in records.items():
   k=epsilon,key
   if k in needed:assert k not in origin and digest(cache[epsilon][key])==digest(record);origin[k]=name
 stage('H64_new070',S/'upper_cover070/batch0005/state.json','1/64','new_cache',[auditbase/'audit_kernel081/report.md',auditbase/'audit_kernel081/manifest.json'])
 stage('H64_new089',S/'upper_g105089/batch0009/state.json','1/64','new_cache',[auditbase/'audit_g105090/report.md',auditbase/'audit_g105090/manifest.json'])
 stage('H128_new107',S/'upper_g110107/batch0001/state.json','1/128','new_cache_by_epsilon',[legacyaudit/'audit_g110116/report.md',legacyaudit/'audit_g110116/manifest.json'])
 stage('H128_new123',S/'upper_continue123/batch0002/state.json','1/128','new_cache_by_epsilon',[auditbase/'audit_g123_sol/report.md',auditbase/'audit_g123_sol/manifest.json'])
 stage('H128_new137',S/'upper_parallel137/residual_preflight/future_batch0002/state.json','1/128','new_cache_by_epsilon',[auditbase/'audit_residual138/report.md',auditbase/'audit_residual138/manifest.json'])
 stage('H128_new1679',S/'upper_parallel137/next1679_H_preflight/future_batch0001/state.json','1/128','new_cache_by_epsilon',[auditbase/'audit_next143/report.md',auditbase/'audit_next143/manifest.json'])
 stage('H128_new126',lastpath,'1/128','new_cache_by_epsilon',[auditbase/'audit_final147/report.md',auditbase/'audit_final147/manifest.json'])
 stages['H128_old098']={'production_state':str(basepath),'production_state_sha256':locks[str(basepath)],'cache_field':'cache_by_epsilon','epsilon':'1/128','numeric_qualification_documents':[document(auditbase/'audit_mixed100/report.md'),document(auditbase/'audit_mixed100/manifest.json')]}
 kernel=H.parent/'upper_cover046/arb_cover.py';table=D/'Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz';assert sha(table)=='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9';locks[str(table)]=sha(table);kernel_id=document(kernel);model_ref=document(H.parent/'upper_continue123/binding.py');stages['H64_inherited052']={'production_state':str(S/'upper_cover052/batch0043/state.json'),'production_state_sha256':'065382625da1645e267025ef7c26d292a29a62b65cb7ab97f49dfa9c06995dbf','epsilon':'1/64','cache_field':'cache','numeric_qualification_scope':'inherited frozen046 protocol/cache;081 explicitly did not independently recompute old178722','numeric_qualification_documents':[document(H.parent/'upper_cover046/report.md'),document(H.parent/'upper_cover046/registration.json'),kernel_id,document(auditbase/'audit_kernel081/report.md')]}
 for p in [D/'Math-Research-cube5-astra-lower/cube-dim-5-bounds/audit_frame075/report.md',D/'Math-Research-cube5-astra-lower/cube-dim-5-bounds/audit_frame075/manifest.json',D/'Math-Research-cube5-auditcurvature129-sol/cube-dim-5-bounds/audit_curvature129/report.md',D/'Math-Research-cube5-auditcurvature129-sol/cube-dim-5-bounds/audit_curvature129/manifest.json',D/'Math-Research-cube-dim-5/cube-dim-5-bounds/upper_alpha1_029/report.md',D/'Math-Research-cube-dim-5/cube-dim-5-bounds/upper_alpha1_029/manifest.json']:document(p)
 streams={name:Stream(name) for name in ('coefficients','jets','leaves','provenance','supports','documents')};jetids={};origin_counts={};projection_counts={'coefficients':0,'jets':0,'leaves':0,'provenance':0,'supports':0}
 for e,key in sorted(needed):
  record=cache[e][key];stageid=origin.get((e,key),'H128_old098' if e=='1/128' else 'H64_inherited052');origin_counts[stageid]=origin_counts.get(stageid,0)+1;assert record['bits']>0 and len(record['values'])==5
  if 'epsilon' in record:assert F(record['epsilon'])==F(e)
  if 'model' in record:assert record['model']==MODEL
  jid=digest({'epsilon':e,'point_key':key,'raw_record_sha256':rawshas[(e,key)]});jetids[(e,key)]=jid;minimal={'jet_id':jid,'model':MODEL,'epsilon':e,'point':record['point'],'bits':record['bits'],'values':record['values'],'raw_record_sha256':rawshas[(e,key)],'numeric_source_id':stageid,'raw_model_field_present':'model' in record,'model_tag_source':'original record' if 'model' in record else model_ref,'raw_epsilon_field_present':'epsilon' in record,'epsilon_tag_source':'original record' if 'epsilon' in record else 'immutable source outertypedcache 1/64'};minimal['minimal_projection_sha256']=digest(minimal);streams['jets'].row(minimal);projection_counts['jets']+=1
 for o in map(json.loads,open(A/'canonical_objects.jsonl')):
  e,key=o['typed_cache_key'];streams['coefficients'].row({'canonical_id':o['canonical_id'],'core':o['core'],'jet_id':jetids[(e,key)],'raw_record_sha256':o['raw_record_sha256']});projection_counts['coefficients']+=1
 for r in map(json.loads,open(A/'leaf_mapper.jsonl')):
  streams['leaves'].row({'leaf_index':r['leaf_index'],'node':r['node'],'canonical_id':r['canonical_id'],'source_alias_index':r['source_alias_index'],'epsilon':r['epsilon'],'box_lower':r['box_lower'],'box_upper':r['box_upper'],'closed_box':r['closed_box'],'claimed_minimum_log_margin':r['saved_g110_margin']});projection_counts['leaves']+=1
 for r in map(json.loads,open(A/'source_provenance.jsonl')):
  streams['provenance'].row({k:v for k,v in r.items() if k not in ('saved_g110_margin','target110')});projection_counts['provenance']+=1
 for n,r in sorted(supports.items()):
  e,key=r['typed_cache_key'];streams['supports'].row({**r,'jet_id':jetids[(e,key)]});projection_counts['supports']+=1
 documents[sha(table)]={'document_id':sha(table),'original_path':str(table),'original_sha256':sha(table),'encoding':'base64-original-gzip','data':base64.b64encode(table.read_bytes()).decode()}
 for key in sorted(documents):streams['documents'].row(documents[key])
 for stream in streams.values():stream.close()
 assert projection_counts['coefficients']==107530 and projection_counts['leaves']==projection_counts['provenance']==111201 and projection_counts['supports']==3812
 metadata={'schema':'portable_global151_minimal_v1','model':MODEL,'model_source_document':model_ref,'target':'11/100000','epsilon_components':['1/64','1/128'],'counts':projection_counts,'source149_counts':{'objects':107530,'leaves':111201,'roots':52360,'splits':58841},'numeric_stages':stages,'numeric_origin_counts':origin_counts,'kernel_document_id':kernel_id,'bodytable_document_id':sha(table),'qualification_status':'candidate;148/primary acceptance not asserted or imported from a live future','numeric_boundary':'minimal original interval projections linked to original fullraw SHA;attempts omitted so fullraw bytes are not reconstructed','independent_numeric_claim_for_inherited052':False,'no_current_H_gradient_log_corner_evaluation':True,'checker_currently_executed':False};write(H/'metadata.json',metadata)
 for p,h in locks.items():assert sha(p)==h
 files={'layout.json':sha(H/'layout.json'),'metadata.json':sha(H/'metadata.json'),'build151.py':sha(__file__),'check151.py':sha(H/'check151.py')};datasets={k:v.parts for k,v in streams.items()};allparts=[p for rows in datasets.values() for p in rows];total=sum(p['bytes'] for p in allparts);manifest={'schema':'upper_global151_data_manifest_v1','status':'pure_data_portable_g110_candidate_no_new_mathematics','files':files,'datasets':datasets,'source149_manifest':str(source),'source149_manifest_sha256':sha(source),'source149_original_data_hashes':{str(A/n):sha(A/n) for n in ['canonical_objects.jsonl','leaf_mapper.jsonl','source_provenance.jsonl']},'read_only_source_locks':locks,'counts':projection_counts,'compressed_data_bytes':total,'maximum_single_compressed_file_bytes':max(p['bytes'] for p in allparts),'all_files_under100MB':all(p['bytes']<100000000 for p in allparts),'all_parts_actual_decode_SHA_verified':True,'new_H_gradient_kernel_log_corner_calls':0,'checker_or_numeric_reproduction_run':False,'global_g_or_C5_acceptance_claim':False,'production_data_PIDs':{'packager':os.getpid()},'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};write(H/'manifest.json',manifest);write(W/'packing_receipt.json',{'manifest_sha256':sha(H/'manifest.json'),'compressed_bytes':total,'maximum_file_bytes':manifest['maximum_single_compressed_file_bytes'],'actual_decode_SHA_all_passed':True,'checker_run':False,'new_H_log_corner_calls':0,'pid':os.getpid(),'wall_seconds':manifest['wall_seconds']});print(json.dumps({'manifest':str(H/'manifest.json'),'manifest_sha256':sha(H/'manifest.json'),'compressed_bytes':total,'maximum_file_bytes':manifest['maximum_single_compressed_file_bytes'],'counts':projection_counts,'numeric_origin_counts':origin_counts,'pid':os.getpid(),'wall_seconds':manifest['wall_seconds'],'peak_rss_bytes':manifest['peak_rss_bytes'],'H_log_corner_calls':0}),flush=True)
if __name__=='__main__':main()
