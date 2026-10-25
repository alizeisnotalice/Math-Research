"""Portable, SHA-locked Astra001/002 entry; default never reruns proofs.

python -B round006/replay_astra.py [--repo-root /clone] [--run 1|2]
Temporary source copies are removed; frozen source and results stay unchanged.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,shutil,tempfile
HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',type=Path,default=HERE.parents[1]);ap.add_argument('--run',choices=['1','2']);args=ap.parse_args();repo=args.repo_root.resolve();proj=repo/'cube-dim-5-bounds'
 locked={}
 for batch in ['astra_upper001','astra_upper002']:
  rec=json.loads((proj/batch/'verification.json').read_text())
  for raw,h in rec['source_sha256'].items():
   if 'cube-dim-5-bounds/' in raw:
    rel=Path('cube-dim-5-bounds')/raw.split('cube-dim-5-bounds/',1)[1]
   elif Path(raw).name=='AGENTS.md':rel=Path('AGENTS.md')
   else:raise ValueError(raw)
   assert '..' not in rel.parts and not rel.is_absolute();path=repo/rel
   assert sha(path)==h, str(rel);locked[path]=h
 with tempfile.TemporaryDirectory(prefix='cube5-astra-replay-') as name:
  stage=Path(name)
  for batch in ['astra_upper001','astra_upper002']:
   (stage/batch).mkdir()
   for fn in ['verify.py','verification.json','registration.md','report.md']:
    shutil.copy2(proj/batch/fn,stage/batch/fn)
  module=None
  for batch in ['astra_upper001','astra_upper002']:
   path=stage/batch/'verify.py';spec=importlib.util.spec_from_file_location('portable_'+batch,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.MAIN=repo
   if batch.endswith('002'):m.old.MAIN=repo
   if args.run==batch[-1]:module=m
  if module is not None:module.main()
 for path,h in locked.items():assert sha(path)==h
 result={'status':'PASS_SHA_mapping_and_safe_import','locked_sources':len(locked),'full_batch_replayed':args.run,'temporary_copies_removed':True,'frozen_originals_unchanged':True}
 print(json.dumps(result,indent=2))
 if args.run is None:(HERE/'replay_astra_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
