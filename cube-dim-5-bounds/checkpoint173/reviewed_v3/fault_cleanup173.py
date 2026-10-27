"""Bounded cleanup ONLY after an observed fault; never limits normal numerical work."""
from pathlib import Path
import time,os
from util173 import atomic,sha,read
def prefix_snapshot(stage_root):
 rows=[]
 for p in sorted(Path(stage_root).glob('shard*/checkpoint.json')):
  try:r=read(p);rows.append({'path':str(p),'sha256':sha(p),'next_index':r.get('next_index'),'last_hash':r.get('last_hash'),'status':r.get('status')})
  except BaseException as e:rows.append({'path':str(p),'snapshot_error':repr(e)})
 return rows

def cleanup_fault_processes(processes,stop,start,stage_root,receipt_path,fault_trigger,grace_seconds=5.0,terminate_seconds=2.0,kill_seconds=2.0):
 assert fault_trigger,'cleanup deadlines require an observed fault;normal batches have no deadline'
 assert grace_seconds>=0 and terminate_seconds>0 and kill_seconds>0
 began=time.monotonic();stop.set();start.set();receipt={'kind':'fault_only_bounded_cleanup','master_pid':os.getpid(),'fault_trigger':fault_trigger,'normal_work_wall_cap':None,'grace_seconds':grace_seconds,'terminate_seconds':terminate_seconds,'kill_seconds':kill_seconds,'prefix_before_cleanup':prefix_snapshot(stage_root),'actions':[],'status':'fault_diagnostic_saved_before_grace'};atomic(receipt_path,receipt)
 def live():return [p for p in processes if p.is_alive()]
 def wait_for_exit(seconds):
  deadline=time.monotonic()+seconds
  while live() and time.monotonic()<deadline:
   for p in processes:p.join(timeout=0)
   time.sleep(min(.02,max(0,deadline-time.monotonic())))
 wait_for_exit(grace_seconds);receipt['prefix_after_grace']=prefix_snapshot(stage_root)
 for p in live():
  receipt['actions'].append({'action':'terminate','pid':p.pid});p.terminate()
 receipt['status']='terminate_sent_after_grace';atomic(receipt_path,receipt);wait_for_exit(terminate_seconds)
 for p in live():
  receipt['actions'].append({'action':'kill','pid':p.pid});p.kill()
 receipt['status']='kill_sent_if_terminate_unresponsive';atomic(receipt_path,receipt);wait_for_exit(kill_seconds)
 evidence=[]
 for p in processes:
  p.join(timeout=0);evidence.append({'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive(),'actual_join_called':True})
 receipt.update({'status':'all_fault_children_joined' if not live() else 'OS_unresponsive_after_kill_failed_closed','workers':evidence,'all_joined':not live(),'prefix_after_cleanup':prefix_snapshot(stage_root),'wall_seconds_after_fault':time.monotonic()-began});atomic(receipt_path,receipt);return receipt
