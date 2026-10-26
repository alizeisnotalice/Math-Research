from pathlib import Path
from fractions import Fraction as F
from math import ceil,floor
import gzip,json,time,os
import engine as E
H=Path(__file__).resolve().parent;S=E.S;N=1<<160

def main():
 began=time.monotonic();reg=json.loads((H/'registration.json').read_text());state=json.loads(gzip.decompress((S/'state.json.gz').read_bytes()));meta=state['range092'];assert meta['next_rank']==56139 and len(meta['blocks'])==14;v=json.loads((H/'verification.json').read_text());assert v['status']=='pass' and v['nodes']==56139;receipts=[];rows=[]
 for b in range(14):
  tag=f'{b:04d}';r=json.loads((H/('block'+tag+'_receipt.json')).read_text());assert r['status']=='complete' and r['workers_joined'];assert E.sha(H/('block'+tag+'_cells.json.gz'))==r['certificates_sha256'];receipts.append(r);rows+=json.loads(gzip.decompress((H/('block'+tag+'_cells.json.gz')).read_bytes()))
 assert [r['rank'] for r in rows]==list(range(56139));budget=receipts[-1]['after'];E.gz(H/'internal_hierarchy.json.gz',{'internal_nodes':{k:n for k,n in state['nodes'].items() if n['children']},'physical_frontier_source':'../upper_adapt082/batch0001_frontier.json.gz','budget':budget});front=sorted(state['frontier'],key=lambda k:-(state['nodes'][k]['selfU']-state['nodes'][k]['selfL']));E.gz(H/'remaining_priority.json.gz',[{'key':k,'level':state['nodes'][k]['level'],'indices':state['nodes'][k]['indices'],'weight':state['nodes'][k]['weight'],'gap_units':state['nodes'][k]['selfU']-state['nodes'][k]['selfL']} for k in front[:1024]]);samplekeys=sorted({0,4095,56138}|{r['rank'] for r in sorted(rows,key=lambda r:r['old_self_U']-r['new_self_U'],reverse=True)[:8]});E.gz(H/'samples.json.gz',[rows[k] for k in samplekeys]);active=[]
 for r in receipts+[v]:
  for pid in r['worker_pids']:
   try:os.kill(pid,0);active.append(pid)
   except ProcessLookupError:pass
 assert not active;lower=F(budget['LH_truncated_integral']);upper=F(budget['UH_truncated_integral']);I=F(budget['I_upper']);summary={'status':'all56139_new_chain_ratio_range_mean_chord_complete_verified','nodes':56139,'physical_frontier':55819,'geometry_anchors_edges_unchanged':True,'ratio_count':89,'ratio_bits':192,'programs':14,'counts':[r['count'] for r in receipts],'upper_improved_nodes':sum(r['improved_upper_nodes'] for r in receipts),'lower_improved_nodes':sum(r['improved_lower_nodes'] for r in receipts),'new_H_calls':0,'new_points_DD_stencils':0,'old_exp_asinh_calls':0,'budget':budget,'LH_lower_simple':str(F(floor(lower*10**12),10**12)),'UH_upper_simple':str(F(ceil(upper*10**12),10**12)),'I_upper_simple':str(F(ceil(I*10**12),10**12)),'max_remaining_gap_key':front[0],'max_remaining_gap':str(F(state['nodes'][front[0]]['selfU']-state['nodes'][front[0]]['selfL'],N)),'production_walls':[r['wall_seconds'] for r in receipts],'production_wall_sum':sum(r['wall_seconds'] for r in receipts),'independent_readback_seconds':v['seconds'],'all_workers_gone':True,'diagnostic_fault_count':len(list((S/'diagnostics').glob('*.json.gz'))),'state_sha256':E.sha(S/'state.json.gz'),'internal_hierarchy_sha256':E.sha(H/'internal_hierarchy.json.gz'),'priority_sha256':E.sha(H/'remaining_priority.json.gz'),'source088_manifest_sha256':reg['source088_manifest_sha256'],'source091_manifest_sha256':reg['source091_manifest_sha256'],'seconds':time.monotonic()-began,'pid':os.getpid(),'source_sha256':E.sha(__file__)};E.atomic(H/'summary.json',summary);(H/'globalg_acceptance012.json').write_bytes((H.parent/'upper_reenvelope088/globalg_acceptance012.json').read_bytes());print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
