"""Prepared epsilon-aware worker functions; no pool or H runs on import."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,os,resource,time
import epsilon_kernel as api
from cache_model import pair_key
EPSILON=F(1,128);TARGET=F(21,200000)
KERNEL=None;KERNEL_SHA=None
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def atomic(path,obj):
 p=Path(path);q=p.with_suffix(p.suffix+'.tmp')
 with open(q,'w') as f:json.dump(obj,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(q,p)
def initialize(table_path,table_sha,kernel_sha,ready_queue):
 global KERNEL,KERNEL_SHA
 assert sha(table_path)==table_sha and sha(Path(api.__file__))==kernel_sha
 KERNEL_SHA=kernel_sha;KERNEL=api.FirstJetKernel(json.loads(gzip.decompress(Path(table_path).read_bytes())),EPSILON)
 assert KERNEL.epsilon==EPSILON
 ready_queue.put({'pid':os.getpid(),'epsilon':'1/128','kernel_sha256':kernel_sha,'vertices':280,'cones':3504})
def evaluate_job(task):
 assert F(task['epsilon'])==EPSILON and task['pair_cache_key']==list(pair_key(EPSILON,task['point']))
 assert KERNEL is not None and KERNEL.eligible(task['point'])['eligible']
 began=time.monotonic();attempts=[];pi=0;partial=None;calls=0
 def save(stage,**extra):
  atomic(task['point_checkpoint'],{'stage':stage,'epsilon':'1/128','pair_cache_key':task['pair_cache_key'],'point':task['point'],
   'precision_index':pi,'partial':partial,'attempts':attempts,'calls':calls,'pid':os.getpid(),'kernel_sha256':KERNEL_SHA,**extra})
 def result(kind,**extra):return {'kind':kind,'epsilon':'1/128','pair_cache_key':task['pair_cache_key'],'pid':os.getpid(),'calls':calls,'attempts':attempts,'seconds':time.monotonic()-began,**extra}
 try:
  save('registered_before_H128')
  while pi<3:
   bits=(256,384,512)[pi]
   # Cooperative cone snapshot only; there is no program wall cap.
   got=KERNEL.evaluate(task['point'],bits,partial,time.monotonic()+1.0,128);calls+=1
   if not got['complete']:
    partial=got['partial'];save('cone_prefix_H128')
    if Path(task['stop_file']).exists():return result('partial',precision_index=pi,partial=partial)
    continue
   values=got['values'];packed=list(map(api.pack,values));quality=api.numeric_quality(values)
   attempt={k:got[k] for k in ('bits','seconds','cancellation_ratio','minimum_integer_projection_gap','common_denominator','all_cones')};attempt['numeric_quality']=quality;attempts.append(attempt)
   counter=values[0].is_finite() and api.endpoint(values[0],True)<TARGET
   separated=values[0].is_finite() and api.endpoint(values[0])>TARGET
   if counter or quality and separated:
    record={'epsilon':'1/128','point':task['point'],'values':packed,'bits':bits,'candidate':task['candidate'],'geometry_tests':task['geometry_tests'],'attempts':attempts,'kernel_sha256':KERNEL_SHA}
    # This comparison is only a model/normalization consistency guard, never a
    # substitution of an H64 value or gradient for an H128 quantity.
    if task.get('old_H64_interval') and F(record['values'][0][1])<F(task['old_H64_interval'][0]):
     save('cross_epsilon_monotonicity_fault',record=record);return result('cross_epsilon_monotonicity_fault',record=record)
    kind='H128_target_obstruction' if counter else 'complete';save(kind,record=record)
    return result(kind,record=record,peak_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
     scope='H_1/128 target only;not a full G or C5 refutation' if counter else 'fresh H128 jet,not H64 cache reuse')
   attempts[-1]['values']=packed;pi+=1;partial=None;save('numerical_precision_retry_H128')
  save('numerical_diagnostic');return result('numerical_diagnostic',precision_index=pi,partial=partial)
 except BaseException as exc:
  save('actual_worker_fault',exception=repr(exc));return result('worker_error',exception=repr(exc),precision_index=pi,partial=partial)
