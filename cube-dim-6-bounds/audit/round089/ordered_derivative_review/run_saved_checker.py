#!/usr/bin/env python3
"""One frozen saved-proof checker invocation; only output and guards adapted."""
from pathlib import Path
import hashlib,json,resource,signal,sys,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
DEST=HERE/'saved_checker_review_receipt.json';RUN=HERE/'saved_checker_execution_receipt.json'
if DEST.exists() or RUN.exists():raise FileExistsError('Single registered checker invocation already recorded')
CPU=time.process_time();WALL=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(15,15))
def wall(*_):raise TimeoutError('180-second total wall cap including cleanup')
signal.signal(signal.SIGALRM,wall);signal.alarm(180)
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-orderstat-astra/cube-dim-6-bounds/lower/round089/ordered_a1_derivative')
SRC=AUTHOR/'check_saved_complete.py';code=SRC.read_text()
sha=hashlib.sha256(SRC.read_bytes()).hexdigest()
assert sha=='cc56c5caccf8e3be882a10c6fbc59572698781559839511cd5ff726683046aa1'
gate=json.loads((HERE/'source_gate.json').read_text())
for p,h in gate['source_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
patches=[('resource.setrlimit(resource.RLIMIT_CPU,(20,20))','resource.setrlimit(resource.RLIMIT_CPU,(15,15))'),('time.process_time()-CPU<18.5','time.process_time()-CPU<14'),("with (H/'checker_complete_receipt.json').open('x') as f:","with Path("+repr(str(DEST))+").open('x') as f:")]
for before,after in patches:
 assert code.count(before)==1;code=code.replace(before,after,1)
adapted_sha=hashlib.sha256(code.encode()).hexdigest()
sys.path.insert(0,str(AUTHOR))
exec(compile(code,str(SRC),'exec'),{'__name__':'__main__','__file__':str(SRC)})
# Source's proof equations are unchanged. Re-arm the remaining total wall guard for receipt cleanup.
signal.alarm(max(1,int(180-(time.monotonic()-WALL))))
r=json.loads(DEST.read_text());assert r['status']=='pass_saved_certificate'
for p,h in gate['source_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
receipt={'status':'pass_once_frozen_saved_checker','source_checker_sha256':sha,'adapted_execution_text_sha256':adapted_sha,'only_adaptations':[{'from':a,'to':b} for a,b in patches],'source_math_equations_changed':False,'output_receipt_path':str(DEST),'output_receipt_sha256':hashlib.sha256(DEST.read_bytes()).hexdigest(),'registered_audit_phase':3,'actual_checker_invocations':1,'new_parameter_configurations':0,'new_integration_states':0,'author_integral_generation_replays':0,'old087_DP_replays':0,'CPU_seconds_including_checker_save_and_input_recheck':time.process_time()-CPU,'wall_seconds_including_checker_save_and_input_recheck':time.monotonic()-WALL,'all20_frozen_sources_unchanged':True}
with RUN.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({'audit_phase':3,'status':receipt['status'],'CPU_seconds_including_save':time.process_time()-CPU,'wall_seconds_including_save':time.monotonic()-WALL,'new_integration_states':0},indent=2))
signal.alarm(0)
