"""Authorized bounded full old-tree consumption, separate <=4096 registrations."""
from pathlib import Path
import json,hashlib,subprocess,sys,os,time
H=Path(__file__).resolve().parent
S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_path112_astra')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,r):Path(p).write_text(json.dumps(r,indent=2)+'\n')
if __name__=='__main__':
 t=time.monotonic();p=H/'continuation_registration.json'
 if not p.exists():
  reg=read(H/'scalar_registration.json');src=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-alt-cont/cube-dim-5-bounds/upper_adaptive103');a=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/upper_pathmean111');audit=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-cover/cube-dim-5-bounds/audit_adaptive113');assert sha(src/'manifest.json')=='8abb73101ed1b328c0ee658a835dc22792379bbfa4297e8a842fda6276d36cb6';assert sha(a/'manifest.json')=='daa8051cae6ff6b4e5b0c65cc82f01943791f1af55567d93f86fea6d4b9d3083';assert sha(audit/'manifest.json')=='1532819753a6c6c5b82053f4ee6d316e97e9b6c71ed3901291199cf7050b3d54'
  inputs=[src/'manifest.json',src/'registration.json',src/'sourceg105_primary_acceptance.json',a/'manifest.json',audit/'manifest.json',H/'production_v2.py',H/'nodes.json.gz',H/'factor_catalog.json',H/'pilot_receipt.json'];write(p,{'inputs':[{'path':str(x),'sha256':sha(x)} for x in inputs],'full_nodes':97972,'remaining_after_pilot':97716,'max_per_batch':4096,'max_new_batches':24,'budget_E':read(src/'registration.json')['E076_upper'],'g105':'21/200000','new_H_R_jet_log_exp_asinh':0,'source_tree_fixed':True,'max_workers':1,'code_sha256':sha(__file__),'pid':os.getpid(),'unix':time.time()})
 reg=read(p);assert sha(__file__)==reg['code_sha256'];receipt={'status':'running','pid':os.getpid(),'registration_sha256':sha(p),'processes':[]};write(H/'continuation_receipt.json',receipt)
 for _ in range(24):
  if (H/'first_strict32.json').exists():receipt['status']='stopped_at_first_strict32_pending_review';break
  for x in reg['inputs']:assert sha(x['path'])==x['sha256']
  c=subprocess.Popen([sys.executable,str(H/'production_v2.py'),'block']);receipt['processes'].append({'pid':c.pid,'began_unix':time.time()});write(H/'continuation_receipt.json',receipt);status=c.wait();receipt['processes'][-1]['exit_code']=status
  if status:receipt['status']='stopped_on_real_fault';break
  if (H/'first_strict32.json').exists():receipt['status']='stopped_at_first_strict32_pending_review';break
  rr=sorted(H.glob('block*_receipt.json'))
  if rr and read(rr[-1])['processed_total']==97972:receipt['status']='all97972_complete';break
 receipt['seconds']=time.monotonic()-t;write(H/'continuation_receipt.json',receipt);print(json.dumps(receipt),flush=True)
