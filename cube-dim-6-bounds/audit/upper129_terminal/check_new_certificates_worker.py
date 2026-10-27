"""A finite shard of only C129-new A057 corner endpoints, independent Arb768 reconstruction."""
import sys,os,json,hashlib,itertools,math,time,traceback
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
for n in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[n]='1'
from pathlib import Path
from fractions import Fraction as Q
from flint import arb,fmpq,ctx
ctx.prec=768
P=Path(__file__).resolve().parent;A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round129/recovered_continuation')
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ck(v,label):
 if not v:raise AssertionError(label)
def canon(v):return(json.dumps(v,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()
def dyadic(x):
 ck(set(x)=={'lower_mantissa','upper_mantissa','common_binary_exponent'},'canonical dyadic fields');n=[]
 for key in ['lower_mantissa','upper_mantissa','common_binary_exponent']:
  s=x[key];v=int(s);ck(str(v)==s,'canonical integer');n.append(v)
 l,u,e=n;ck(l<=u,'ordered dyadic');return(l*Q(2)**e,u*Q(2)**e)
def read_selected(ref):
 ck(sha(ref['path'])==ref['sha256'],'frozen selected source hash');x=rd(ref['path'])
 for k in ref.get('selector',[]):x=x[k]
 return x
REG=rd(A/'registration_129.json');levelpins=REG['binaries'];ORIGINAL=rd(REG['inputs']['117_inputs_path']);SOURCE_PID=rd(A/'started_v129.json')['pid']

def verify_source(job):
 if job['source_kind']=='FROZEN_BASELINE':
  m=job['frozen_source_metadata'];jet=read_selected(job['jet_ref']);ck(jet['center']==m['expected']['center']and jet['window']==m['expected']['window']and jet['provenance']['native_sha256']==m['native_sha256'],'unchanged accepted baseline jet only')
  return jet,m['expected'],m['artifact_key'],False
 state=read_selected(job['state_ref']);k=job['semantic_key'];e=state['expected'];ck(state['status']=='AVAILABLE'and state['semantic_key']==k and state['artifact_key']=='C129.'+k and len(state['levels'])==1,'new source state')
 ck(sha(state['native_path'])==state['native_sha256']and sha(state['jet_path'])==job['jet_sha256']and sha(state['receipt_path'])==job['receipt_sha256'],'new raw/jet/receipt hashes')
 m=rd(state['native_path']);jet=rd(state['jet_path']);receipt=rd(state['receipt_path']);intent=read_selected({'path':job['intent_path'],'sha256':job['intent_sha256']});process=read_selected({'path':job['process_path'],'sha256':job['process_sha256']});bits=state['accepted_precision_bits'];b=levelpins[str(bits)]
 semantic={'center':[str(Q(z))for z in e['center']],'a':str(Q(e['window']['a'])),'T':str(Q(e['window']['T'])),'formula_family_id':'C016-cauchy-frame32-outward-precision-family-v1','fan_sha256':next(h for p,h in REG['runtime_and_fan_pins'].items()if p.endswith('fan_integer_input.txt')),'model':e['model'],'model_definition':e['model_definition'],'quadratic':e['quadratic']}
 ck(hashlib.sha256(canon(semantic)).hexdigest()==k,'new semantic physical identity')
 ck(intent['expected']==e and intent['semantic_key']==receipt['semantic_key']==k and intent['precision_bits']==receipt['precision_bits']==bits and receipt['exit_code']==0,'new intent/receipt metadata')
 nums=[Q(z)for z in e['center']];den=math.lcm(*(z.denominator for z in nums));pointints=[int(z*den)for z in nums];w=e['window'];cmd=[b['path'],'full',next(p for p in REG['runtime_and_fan_pins']if p.endswith('fan_integer_input.txt')),str(Q(w['a']).numerator),str(Q(w['a']).denominator),str(Q(w['T']).numerator),str(Q(w['T']).denominator),str(den),*map(str,pointints)]
 ck(intent['command']==receipt['command']==process['command']==cmd,'exact original native command')
 ck(receipt['stdout_sha256']==receipt['raw_sha256']==state['native_sha256']and receipt['intent_sha256']==job['intent_sha256']and receipt['process_sha256']==job['process_sha256'],'original raw command chain')
 ck(set(m)=={'status','precision_bits','full_group_terms','computed_sign_pairs','confluent_terms','zero_node_terms','unique_node_values','kernel_seconds','band','gradient','p_a','p_T','Euler_xgrad_plus5B_plusTpT_minusapa'},'exact native fields')
 ck(e['model']=='body32_unit_cauchy_v1'and e['quadratic']=='frame32'and jet['provenance']['gradient_gate_sha256']==ORIGINAL['F061_manifest_sha256'],'original model/gradient gate')
 ck(process['wrapper_pid']==SOURCE_PID and process['pid']==state['levels'][0]['pid'],'original source PID')
 ck(m['status']=='all_real_nodes_supported'and m['precision_bits']==bits==128 and m['full_group_terms']==967680 and m['computed_sign_pairs']==483840,'new complete native family')
 ck(len(m['gradient'])==6 and all(type(m[z])is int and 0<=m[z]<=967680 and m[z]%2==0 for z in ['confluent_terms','zero_node_terms'])and 1<=m['unique_node_values']<=2903040 and math.isfinite(m['kernel_seconds'])and m['kernel_seconds']>=0,'new native schema/counts')
 F=dyadic(m['band']);g=[dyadic(v)for v in m['gradient']];pa=dyadic(m['p_a']);pt=dyadic(m['p_T']);ne=dyadic(m['Euler_xgrad_plus5B_plusTpT_minusapa']);ck(F[0]>0 and pa[0]>0 and pt[0]>0 and ne[0]<=0<=ne[1],'positive band/endpoints/native Euler')
 terms=[(5*F[0],5*F[1]),(Q(w['T'])*pt[0],Q(w['T'])*pt[1]),(-Q(w['a'])*pa[1],-Q(w['a'])*pa[0])]+[tuple(sorted(Q(c)*z for z in v))for c,v in zip(e['center'],g)]
 ck(sum(v[0]for v in terms)<=0<=sum(v[1]for v in terms),'independent exact six-gradient Euler')
 intervals=[F,*g,pa,pt,ne];report=[]
 for j,(l,u)in enumerate(intervals):
  width=u-l;limit=Q(1,65000*(1<<40))*max(Q(1),abs(l),abs(u));ck(width<=limit,'all10 original strict width');report.append({'index':j,'width':str(width),'limit':str(limit)})
 ck(receipt['width_acceptable']is True and receipt['width_report']=={'widths':report,'failed_indices':[]}and state['levels'][0]['validated_intervals']==[[str(l),str(u)]for l,u in intervals] and state['levels'][0]['width_report']==receipt['width_report'],'new saved exact width report')
 ck(jet['F_interval']==list(map(str,F))and jet['face_gradient_intervals']==[list(map(str,v))for v in g[1:]],'new jet raw value/five physical gradients')
 ck(all(jet[k]==e[k]for k in ['model','model_definition','center','window'])and jet['provenance']['synthetic_fixture']is False and jet['provenance']['precision_bits']==bits and jet['provenance']['native_sha256']==state['native_sha256']and jet['provenance']['receipt_sha256']==job['receipt_sha256']and jet['provenance']['binary_sha256']==b['sha256']and jet['provenance']['source_sha256']==b['source_sha256'],'new true jet provenance')
 return jet,e,state['artifact_key'],True

def verify(job):
 jet,source,artifact,isnew=verify_source(job);entry=job['certificate']
 if entry is None:return{'new_native_records':int(isnew),'frozen_old_jets':0,'certificates':0,'corners':0}
 ck(sha(entry['path'])==entry['sha256'],'saved new certificate hash');saved=rd(entry['path']);e=entry['expected']
 ck(entry['native_artifact_key']==artifact and e['center']==source['center']and e['window']==source['window'],'new geometry matching true source jet')
 ck(e['target']=='1/65000'and e['quadratic']==saved['quadratic']=='frame32'and saved['model']=='body32_unit_cauchy_v1'and saved['window']==e['window']and saved['target']==e['target']and saved['precision_bits']==384 and saved['conditional_on_jet_source']==jet['source_proof'],'new A057 source/target/model')
 center=list(map(Q,e['center']));box=[tuple(map(Q,b))for b in e['face_box']];L,U=map(Q,jet['F_interval']);grad=[tuple(map(Q,g))for g in jet['face_gradient_intervals']];a=Q(e['window']['a']);ck(L>0 and len(box)==len(grad)==5 and len(saved['corners'])==32,'wholebox dimensions')
 lo=[];hi=[]
 for bits,row in zip(itertools.product((0,1),repeat=5),saved['corners']):
  d=[box[i][bit]-center[i+1]for i,bit in enumerate(bits)];dot=sum(min(g[0]*z,g[1]*z)for g,z in zip(grad,d));f=min(U,max(L,dot));q=sum(z*z for z in d)/32;exponent=dot/f-q/a**2
  ck(row['corner_bits']==list(bits)and list(map(Q,row['displacement']))==d and Q(row['gradient_dot_lower'])==dot and Q(row['minimizing_center_value'])==f and Q(row['quadratic'])==q and Q(row['exponent'])==exponent,'independent new exact A057 corner')
  value=arb(fmpq(f.numerator,f.denominator))*arb(fmpq(exponent.numerator,exponent.denominator)).exp();ck(value.is_finite(),'finite new Arb768');ends=[]
  for z in(value.lower(),value.upper()):m,p=z.man_exp();ends.append(Q(int(m))*Q(2)**int(p))
  l,u=Q(row['certificate_lower']),Q(row['certificate_upper']);ck(l<=ends[0]<=ends[1]<=u and row['passes_target']==(l>=Q(1,65000)),'saved new384 encloses independent768')
  lo.append(l);hi.append(u)
 ck(Q(saved['tangent_bound_lower'])==min(lo)and Q(saved['tangent_bound_upper'])==min(hi)and saved['status']==('PASS_CONDITIONAL_ON_MATCHING_CENTER_JET'if min(lo)>=Q(1,65000)else'UNRESOLVED')and saved['global_g_certificate']is False and saved['kernel_calls']==0,'new wholebox extrema/scope')
 return{'new_native_records':int(isnew),'frozen_old_jets':int(not isnew),'certificates':1,'corners':32}

def main():
 shard=int(sys.argv[1]);ck(0<=shard<15,'finite shard ID');release=rd(P/'parent_release.json');ck(release.get('authorized')is True and release['worker_script_sha256']==sha(__file__),'parent audit launch required');prepared=rd(P/'prepared_new_jobs.json');spec=prepared['shards'][shard];ck(sha(spec['path'])==spec['sha256'],'registered shard bytes');ck(os.getpriority(os.PRIO_PROCESS,0)>=10,'nice>=10')
 stop=P/'STOP_REQUESTED.json';outdir=P/'worker_results';outdir.mkdir(exist_ok=True)
 with(outdir/f'start_{shard:02d}.json').open('x')as f:json.dump({'pid':os.getpid(),'shard':shard,'automatic_retry':False},f)
 counter={'jobs':0,'new_native_records':0,'frozen_old_jets':0,'certificates':0,'corners':0};start=time.monotonic();failure=None
 try:
  for line in Path(spec['path']).read_text().splitlines():
   ck(not stop.exists(),'another shard failed; no automatic retry');j=json.loads(line);r=verify(j)
   counter['jobs']+=1
   for k,v in r.items():counter[k]+=v
   if counter['jobs']%100==0:(outdir/f'progress_{shard:02d}.json').write_text(json.dumps({'shard':shard,'pid':os.getpid(),**counter})+'\n')
 except BaseException:
  failure=traceback.format_exc()
  try:
   with stop.open('x')as f:json.dump({'shard':shard,'failure':failure},f)
  except FileExistsError:pass
 finally:(outdir/f'result_{shard:02d}.json').write_text(json.dumps({'status':'PASS'if failure is None else'FAILED_NO_RETRY','shard':shard,'pid':os.getpid(),**counter,'failure':failure,'CPU_seconds':time.process_time(),'wall_seconds':time.monotonic()-start,'old_baseline_native_or_Arb_replay':0,'new_precision_bits':768},indent=2)+'\n')
 if failure:raise SystemExit(1)
if __name__=='__main__':main()
