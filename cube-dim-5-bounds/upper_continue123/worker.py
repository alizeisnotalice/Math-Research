"""Fresh H128 jet wrapper, parameterized by the reviewed g110 target.

No pool or H is evaluated on import. Numerical FirstJetKernel is098 unchanged.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,os,resource,time
from binding import api,KERNEL_PATH,KERNEL_SHA,MODEL
from cache_model import pair_key
EPSILON=F(1,128);KERNEL=None;RUN_TARGET=None
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def atomic(path,obj,fresh=False):
 p=Path(path)
 if fresh:
  with open(p,'x') as f:json.dump(obj,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 else:
  q=p.with_suffix(p.suffix+'.tmp')
  with open(q,'w') as f:json.dump(obj,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
  os.replace(q,p)
def initialize(table_path,table_sha,target,ready_queue):
 global KERNEL,RUN_TARGET
 RUN_TARGET=F(target);assert RUN_TARGET==F(11,100000)
 assert sha(table_path)==table_sha and sha(KERNEL_PATH)==KERNEL_SHA
 os.nice(5) # Yield CPU scheduling priority to the concurrent integral workers.
 KERNEL=api.FirstJetKernel(json.loads(gzip.decompress(Path(table_path).read_bytes())),EPSILON)
 assert KERNEL.epsilon==EPSILON
 ready_queue.put({'pid':os.getpid(),'epsilon':str(EPSILON),'target':str(RUN_TARGET),'kernel_sha256':KERNEL_SHA,'vertices':280,'cones':3504,'nice':os.getpriority(os.PRIO_PROCESS,0)})
def validate_task(task):
 target=F(task['target']);assert target==RUN_TARGET==F(11,100000)
 assert F(task['epsilon'])==EPSILON and task['pair_cache_key']==list(pair_key(EPSILON,task['point']))
 guard=task.get('same_point_H64_guard')
 if task['candidate']!=task['original_H64_candidate']:assert guard is None,'changed anchor cannot use old-point H64 monotonicity guard'
 if guard is not None:assert guard['epsilon']=='1/64' and pair_key(EPSILON,guard['point'])[1]==task['pair_cache_key'][1]
 return target
def evaluate_job(task):
 target=validate_task(task);assert KERNEL is not None and KERNEL.eligible(task['point'])['eligible']
 assert not Path(task['point_checkpoint']).exists(),'Existing checkpoint requires fingerprint inspection and separate recovery authorization'
 began=time.monotonic();attempts=[];pi=0;partial=None;calls=0
 def save(stage,fresh=False,**extra):
  atomic(task['point_checkpoint'],{'stage':stage,'model':MODEL,'epsilon':str(EPSILON),'target':str(target),'pair_cache_key':task['pair_cache_key'],'point':task['point'],
   'precision_index':pi,'partial':partial,'attempts':attempts,'calls':calls,'pid':os.getpid(),'kernel_sha256':KERNEL_SHA,**extra},fresh)
 def result(kind,**extra):return {'kind':kind,'epsilon':str(EPSILON),'target':str(target),'pair_cache_key':task['pair_cache_key'],'pid':os.getpid(),'calls':calls,'attempts':attempts,'seconds':time.monotonic()-began,**extra}
 try:
  save('registered_before_H128_g110',fresh=True)
  while pi<3:
   bits=(256,384,512)[pi]
   got=KERNEL.evaluate(task['point'],bits,partial,time.monotonic()+1.0,128);calls+=1
   if not got['complete']:
    partial=got['partial'];save('cone_prefix_H128_g110')
    if Path(task['stop_file']).exists():return result('partial',precision_index=pi,partial=partial)
    continue
   partial=None
   values=got['values'];packed=list(map(api.pack,values));quality=api.numeric_quality(values)
   attempt={k:got[k] for k in ('bits','seconds','cancellation_ratio','minimum_integer_projection_gap','common_denominator','all_cones')};attempt['numeric_quality']=quality;attempts.append(attempt)
   counter=values[0].is_finite() and api.endpoint(values[0],True)<target
   separated=values[0].is_finite() and api.endpoint(values[0])>target
   if counter or quality and separated:
    record={'model':MODEL,'epsilon':str(EPSILON),'qualification_target':str(target),'point':task['point'],'pair_cache_key':task['pair_cache_key'],'values':packed,'bits':bits,
     'candidate':task['candidate'],'geometry_tests':task['geometry_tests'],'attempts':attempts,'kernel_sha256':KERNEL_SHA}
    guard=task.get('same_point_H64_guard')
    if guard is not None and F(record['values'][0][1])<F(guard['H_interval'][0]):
     save('same_point_cross_epsilon_fault',record=record);return result('same_point_cross_epsilon_fault',record=record)
    kind='H128_g110_target_obstruction' if counter else 'complete';save(kind,record=record)
    return result(kind,record=record,peak_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
     scope='fixed H128 at target only;not a G/C5 refutation' if counter else 'new H128 physical jet')
   attempts[-1]['values']=packed;pi+=1;partial=None;save('numerical_precision_retry_H128_g110')
  save('numerical_diagnostic');return result('numerical_diagnostic',precision_index=pi,partial=partial)
 except BaseException as exc:
  save('actual_worker_fault',exception=repr(exc));return result('worker_error',exception=repr(exc),precision_index=pi,partial=partial)
