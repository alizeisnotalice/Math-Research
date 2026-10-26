from pathlib import Path
from fractions import Fraction as F
from math import ceil
import gzip,json,time,os
import engine as E
H=Path(__file__).resolve().parent;S=E.S;N=1<<160

def main():
 t=time.monotonic();reg=json.loads((H/'registration.json').read_text());state=json.loads(gzip.decompress((S/'state.json.gz').read_bytes()));meta=state['reenvelope088'];assert meta['next_rank']==56139 and len(meta['blocks'])==14;verification=json.loads((H/'verification.json').read_text());assert verification['status']=='pass' and verification['nodes_verified']==56139;rows=[];receipts=[]
 for b in range(14):
  tag=f'{b:04d}';r=json.loads((H/('block'+tag+'_receipt.json')).read_text());assert r['status']=='complete' and r['workers_joined'];assert E.sha(H/('block'+tag+'_cells.json.gz'))==r['certificates_sha256'];receipts.append(r);rows+=json.loads(gzip.decompress((H/('block'+tag+'_cells.json.gz')).read_bytes()))
 assert [r['rank'] for r in rows]==list(range(56139));after=receipts[-1]['after'];internal={k:n for k,n in state['nodes'].items() if n['children']};E.gz(H/'internal_hierarchy.json.gz',{'source082_manifest_sha256':reg['source082_manifest_sha256'],'internal_nodes':internal,'physical_frontier_source082_file':'../upper_adapt082/batch0001_frontier.json.gz','physical_frontier_count':55819,'global_budget':after});samples=sorted({0,4095,56138}|{r['rank'] for r in sorted(rows,key=lambda r:r['old_self_U']-r['new_self_U'],reverse=True)[:8]});E.gz(H/'samples.json.gz',[rows[r] for r in samples]);active=[]
 for r in receipts+[verification]:
  for pid in r['worker_pids']:
   try:os.kill(pid,0);active.append(pid)
   except ProcessLookupError:pass
 assert not active;copyg=Path(reg['g_acceptance_path']).read_bytes();assert __import__('hashlib').sha256(copyg).hexdigest()==reg['g_acceptance_sha256'];(H/'globalg_acceptance012.json').write_bytes(copyg)
 summary={'status':'all56139_zero_H_reenvelopes_complete_verified','node_count':56139,'physical_frontier':55819,'geometry_and_edges_unchanged':True,'lower_bounds_unchanged':True,'processed_programs':14,'counts_per_program':[r['count'] for r in receipts],'improved_self_U_nodes':sum(r['improved_nodes'] for r in receipts),'new_H_calls':0,'new_anchors':0,'new_stencils_or_DD':0,'moment_source':'085 face25/24,mean25/48 Emax²','g_qualified':'1/10000','g105_used':False,'g_acceptance_sha256':reg['g_acceptance_sha256'],'budget':after,'I_upper_simple':str(F(ceil(F(after['I_upper'])*10**12),10**12)),'UH_upper_simple':str(F(ceil(F(after['UH_truncated_integral'])*10**12),10**12)),'production_program_walls':[r['wall_seconds'] for r in receipts],'production_wall_sum':sum(r['wall_seconds'] for r in receipts),'full_independent_readback_seconds':verification['seconds'],'all_worker_pids_gone':True,'state_sha256':E.sha(S/'state.json.gz'),'internal_hierarchy_sha256':E.sha(H/'internal_hierarchy.json.gz'),'samples_sha256':E.sha(H/'samples.json.gz'),'seconds':time.monotonic()-t,'pid':os.getpid(),'source_sha256':E.sha(__file__)};E.atomic(H/'summary.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
