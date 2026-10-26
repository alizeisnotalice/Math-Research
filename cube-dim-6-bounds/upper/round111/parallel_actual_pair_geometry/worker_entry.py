"""One authorized k/slab worker; stdout is private binary IPC, no file writes."""
import sys
sys.dont_write_bytecode=True
import hashlib,json,os,time,traceback
from pathlib import Path
from protocol import read_frame,send_frame
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 boot,_=read_frame(sys.stdin.buffer);assert boot['command']=='BOOT'
 rel=json.loads((P/'parent_release.json').read_text());reg=json.loads((P/'registration.json').read_text())
 assert boot['release_sha256']==sha(P/'parent_release.json')and rel['authorized']is True and rel['round']=='111'
 for name in ['math_core.py','protocol.py','worker_entry.py']:assert sha(P/name)==reg['code_sha256'][name]
 assert rel['registration_sha256']==sha(P/'registration.json')
 job=next(j for j in reg['jobs']if j['id']==boot['job_id']);assert job['k']in(2,3,4,5)
 for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:assert os.environ.get(name)=='1'
 assert sha(P/'geometry_certificate.json')==boot['geometry_sha256']
 geometry=json.loads((P/'geometry_certificate.json').read_text())['by_k'][str(job['k'])]
 from math_core import CellEngine
 engine=CellEngine(job,geometry);st=None;phase=1;previous_snapshot=None;start=time.monotonic();phase_start=None
 try:
  st=engine.initial_state();assert st['initial_leaves']==job['initial_leaves']
  while True:
   phase_start=time.monotonic();entry={'leaves':len(st['leaves']),'cell_bounds':st['cell_bounds'],'nextid':st['nextid'],'refinements':len(st['refinements'])}
   def progress(state):
    full=phase==2 and state['splits']%4096==0
    data={'kind':'partial_state'if full else'progress','job_id':job['id'],'phase':phase,'leaves':len(state['leaves']),'cell_bounds':state['cell_bounds'],'lower':state['lower'],'upper':state['upper'],'elapsed_seconds':time.monotonic()-start}
    if full:data['state']=engine.snapshot(state)
    send_frame(sys.stdout.buffer,data)
   target=reg['phase_limits'][str(phase)]['leaves_per_job']
   engine.run_to(st,target,reg['max_total_depth'],lambda:False,progress)
   st['completed_phases'].append({'phase':phase,'target_leaves':target,'entry_counts':entry,'exit_leaves':len(st['leaves']),'exit_cell_bounds':st['cell_bounds'],'outcome':st['last_outcome'],'continued_from_transport_sha256':previous_snapshot})
   message={'kind':'phase_done','job_id':job['id'],'phase':phase,'state':engine.snapshot(st),'leaves':len(st['leaves']),'cell_bounds':st['cell_bounds'],'lower':st['lower'],'upper':st['upper'],'elapsed_seconds':time.monotonic()-start,'phase_elapsed_seconds':time.monotonic()-phase_start,'CPU_process_seconds_before_transport':time.process_time(),'pid':os.getpid(),'threads':1,'continued_from_transport_sha256':previous_snapshot}
   previous_snapshot=send_frame(sys.stdout.buffer,message)
   command,_=read_frame(sys.stdin.buffer)
   if command['command']=='STOP':return
   assert command['command']=='CONTINUE'and phase==1
   assert command['previous_snapshot_transport_sha256']==previous_snapshot
   assert command['next_phase']==2;phase=2
 except Exception:
  try:send_frame(sys.stdout.buffer,{'kind':'failure','job_id':job['id'],'phase':phase,'pid':os.getpid(),'failure':traceback.format_exc(),'state':None if st is None else engine.snapshot(st)})
  except Exception:pass
  raise
if __name__=='__main__':main()
