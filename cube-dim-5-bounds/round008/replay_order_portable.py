"""Hash-only by default; explicit --run rebuilds one recorded batch in a clean temp tree."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,tempfile,shutil,ast
P=Path(__file__).resolve().parents[1]
ENTRIES={'order001':('upper_order001','order_certificate.py'),'order002':('upper_order002','top2_certificate.py'),'order003':('upper_order003','top3_pilot.py'),'adaptive004':('upper_adaptive004','adaptive_pilot.py'),'resolvent005':('upper_resolvent005','resolvent_pilot.py')}
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('batch',choices=ENTRIES);parser.add_argument('--run',action='store_true');args=parser.parse_args()
 chosen,entry=ENTRIES[args.batch];needed={chosen,'upper_order002'} if chosen not in {'upper_order001','upper_order002'} else {chosen};files=[]
 for directory in sorted(needed):
  receipt=json.loads((P/directory/'completion_receipt.json').read_text())
  for row in receipt['white_list']:
   f=P/directory/row['file'];assert hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256']
   if f.suffix=='.py':ast.parse(f.read_text(),filename=str(f))
   files.append(f)
 assert (P/chosen/entry).is_file()
 if not args.run:
  print(json.dumps({'status':'PASS_locked_layout_and_AST_only','batch':args.batch,'files':len(files),'calculation_executed':False,'run_omits_saved_checkpoints_and_regenerates_omitted_NPZ':True}));return
 # The frozen checkpoint records excluded reproducible NPZ files. Do not copy it into a fresh run.
 with tempfile.TemporaryDirectory(prefix='cube5-order-new-replay-') as tmp:
  target=Path(tmp)/'cube-dim-5-bounds'
  for f in files:
   if 'checkpoint' in f.name:continue
   to=target/f.relative_to(P);to.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,to)
  proc=subprocess.run([sys.executable,'-B',str(target/chosen/entry)],cwd=tmp,check=False)
 assert not Path(tmp).exists()
 print(json.dumps({'batch':args.batch,'new_execution_exit_code':proc.returncode,'frozen_sources_modified':False,'temporary_tree_released':True}))
 if proc.returncode:raise SystemExit(proc.returncode)
if __name__=='__main__':main()
