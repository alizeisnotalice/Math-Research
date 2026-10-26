#!/usr/bin/env python3
"""Preregister the next actual finite prefix; performs no scientific evaluation."""
from pathlib import Path
import sys,json,hashlib,datetime
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
if (P/'positive_candidate_LOCK.json').exists():
    print(json.dumps({'stop':'positive_candidate_LOCK','lock':load(P/'positive_candidate_LOCK.json')}));raise SystemExit(0)
state=load(P/'state.json');common=load(P/'registration.json');assert sha(P/'exact_cells.py')==common['source_sha256']
first=state['next_cell'];total=state['total_cells']
if first>=total:print(json.dumps({'stop':'all_cells_complete','covered':first}));raise SystemExit(0)
count=min(int(sys.argv[1]) if len(sys.argv)>1 else 64,64,total-first);assert count>0
last=first+count-1;tag=f'{first:04d}_{last:04d}';name=f'batch_{tag}.registration.json';path=P/'batches'/name
if path.exists():
    old=load(path);assert old['expected_state_sha256']==sha(P/'state.json') and old['source_sha256']==common['source_sha256']
else:
    previous=sorted((P/'batches').glob('batch_*.receipt.json'))
    row={'status':'PREREGISTERED_next_actual_exact_PWA_prefix','registered_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':common['source_sha256'],'expected_state_sha256':sha(P/'state.json'),'interface_sha256':sha(P/'interface.json.gz'),'first_cell':first,'last_cell':last,'new_cells_limit':count,'descriptor_dual_union_objects':count*91,'previous_receipt_sha256':sha(previous[-1]) if previous else None,'cell_artifact':f'cells_{tag}.json.gz','receipt':f'batch_{tag}.receipt.json','scope':'oneexactA/Bcomputationpernewcell;nointerpolationorparametergrid','wall_limit':None,'nice_minimum':10,'threads':1,'stop':'newpositiveLOCKoractualfault;noautomaticrestart'}
    path.write_text(json.dumps(row,separators=(',',':'))+'\n')
print(json.dumps({'registration':name,'first_cell':first,'last_cell':last,'count':count,'registration_sha256':sha(path)}))
