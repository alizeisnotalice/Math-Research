"""Independent readback of all exact budget powers; noH and no roots recalculated."""
from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib,time,os
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive103');N=1<<160

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
def main():
 t=time.monotonic();reg=read(H/'budget_g105_registration.json');assert sha(H/'budget_g105.py')==reg['consumer_sha256'] and sha(H/'registration.json')==reg['production_registration_sha256'];gp=Path(reg['g105_acceptance_path']);assert sha(gp)==reg['g105_acceptance_sha256'];assert read(gp)['qualified_g']=='21/200000';p=read(H/'registration.json');E=F(p['E076_upper']);g=F(21,200000);rows=read(H/'g105_prefix_budgets.json.gz');evs=[json.loads(x) for x in (S/'parent_events.jsonl').read_text().splitlines()];assert len(rows)==len(evs)==2048;previousU=None;positives=[]
 for row,ev in zip(rows,evs):
  assert row['parent_count']==ev['seq'] and row['event_hash']==ev['event_hash'];U=F(ev['global_U_units'],N);assert F(row['UH'])==U and (previousU is None or U<=previousU);previousU=U;B,I,C=map(F,[row['B_upper'],row['I_upper'],row['cost_upper']]);assert B>E and (B-E)**5>=U**4 and I**4>=B**5 and C**4>=(5*B/g)**5;delta=32**4*g**5-5**5*B**5;assert delta==F(row['exact_32_gate_margin']) and (delta>0)==row['strict32'];assert row['valid_g']=='21/200000' and row['no_tube_added']
  if delta>0:positives.append(row['parent_count'])
 result=read(H/'g105_budget_result.json');assert result['status']=='complete' and result['all_complete_parent_prefixes_checked']==2048 and result['final']==rows[-1];assert result['first_exact_positive_event']==(positives[0] if positives else None);assert result['production_kept_g100'] and p['valid_g']=='1/10000';assert rows[-1]['UH']==read(H/'run_receipt.json')['after']['UH_truncated_integral']
 if positives:
  snap=read(H/'g105_first_strict32.json');sp=Path(snap['state_path']);assert sha(sp)==snap['state_sha256'];st=read(sp);meta=st['adaptive103'];assert meta['event_sequence']==snap['parent_count']>=positives[0] and meta['event_hash']==snap['event_hash'];nodes=st['nodes'];assert F(sum(nodes[k]['effU'] for k in st['roots']),N)==F(snap['budget']['UH']);ml=max(nodes[k]['level'] for k in st['frontier']);assert sum(nodes[k]['weight']*2**(4*(ml-nodes[k]['level'])) for k in st['frontier'])==2**(4*ml)
 else:assert not (H/'g105_first_strict32.json').exists()
 out={'status':'pass','all2048_complete_parent_prefixes_exact_powers':True,'g105':'21/200000','first_exact_positive_event':positives[0] if positives else None,'final':rows[-1],'production_g100_locks_unchanged':True,'source_g105_acceptance_sha256':reg['g105_acceptance_sha256'],'new_H_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-t,'source_sha256':sha(__file__)};(H/'g105_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='final'},indent=2))
if __name__=='__main__':main()
