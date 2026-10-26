#!/usr/bin/env python3
"""Assemble completed initial layer only; never starts production/refinement."""
from pathlib import Path
from fractions import Fraction as F
import json,gzip,hashlib,time,os,resource,sys,signal
from math import ceil,floor
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_global080');N=1<<160

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rq(x):
 x=F(str(x));return fmpq(x.numerator,x.denominator)
def ep(x,up=False):return str((x.upper() if up else x.lower()).fmpq())
def alarm(sig,frame):raise TimeoutError('finalize590s')
def main():
 start=time.monotonic();ctx.prec=384;state=json.loads((S/'state.json').read_text());assert state['next_rank']==52360 and state['computed_weight']==1048576;rows=[];receipts=[];readbacks=[]
 for block in state['blocks']:
  bid=f"{block['id']:04d}";bp=S/'blocks'/bid;assert sha(bp/'ledger.json.gz')==block['ledger_sha256'];assert sha(bp/'receipt.json')==block['receipt_sha256'];r=json.loads((bp/'receipt.json').read_text());assert r['workers_joined'] and r['status']=='complete';receipts.append(r);vr=json.loads((H/('block'+bid+'_verification.json')).read_text());assert vr['status']=='pass' and vr['full_bands_checked'] and vr['ledger_sha256']==block['ledger_sha256'];readbacks.append(vr);rows+=json.loads(gzip.decompress((bp/'ledger.json.gz').read_bytes()))
 rows.sort(key=lambda r:r['rank']);assert [r['rank'] for r in rows]==list(range(52360));assert len({r['rootid'] for r in rows})==52360;assert sum(r['weight'] for r in rows)==1048576
 sums=[sum(r['weight']*r['bounds'][k] for r in rows) for k in ('integral_lower_i','integral_upper_i')];assert sums==[state['sum_lower_weighted_units'],state['sum_upper_weighted_units']]
 cols=['rank','rootid','weight','L_i','U_i','M_i','integral_lower_i','integral_upper_i','relative_baseline_upper_i','positive_bands','taylor_bands','exp_bands','global_cap_used','evidence_sha256','mathematical_sha256'];compact=[]
 for r in rows:
  b=r['bounds'];path=S/'evidence'/f"{r['rank']//1024:02d}"/f"{r['rank']:05d}.json.gz";assert sha(path)==r['evidence_sha256'];ev=json.loads(gzip.decompress(path.read_bytes()));ev.pop('runtime');mathsha=hashlib.sha256(json.dumps(ev,sort_keys=True,separators=(',',':')).encode()).hexdigest();compact.append([r['rank'],r['rootid'],r['weight'],*[b[k] for k in cols[3:-2]],r['evidence_sha256'],mathsha])
 payload={'columns':cols,'denominator':str(N),'rows':compact};(H/'initial_cells.json.gz').write_bytes(gzip.compress(json.dumps(payload,separators=(',',':')).encode(),mtime=0));UH=F(sums[1],N);LH=F(sums[0],N);cap=json.loads((H/'cap_scalar.json').read_text());E=F(cap['E076_norm']['upper']);I=(arb(rq(UH))**(arb(4)/5)+arb(rq(E)))**(arb(5)/4);target=arb(rq('1/25000'));UHbudget=(target**(arb(4)/5)-arb(rq(E)))**(arb(5)/4)
 worst=sorted(rows,key=lambda r:r['weight']*(r['bounds']['integral_upper_i']-r['bounds']['integral_lower_i']),reverse=True)[:32];ratio=sorted(rows,key=lambda r:r['bounds']['relative_baseline_upper_i'],reverse=True)[:8]
 sample_ranks=sorted({0,32767,32768,52359}|{r['rank'] for r in worst[:8]}|{r['rank'] for r in rows if r['bounds']['global_cap_used']});samples=[]
 for rank in sample_ranks:
  ev=json.loads(gzip.decompress((S/'evidence'/f'{rank//1024:02d}'/f'{rank:05d}.json.gz').read_bytes()));ev.pop('runtime');samples.append(ev)
 (H/'endpoint_samples.json.gz').write_bytes(gzip.compress(json.dumps({'selection':'first/last/bridge endpoints plus worst8 absolute gaps and all globalcap cases','samples':samples},sort_keys=True,separators=(',',':')).encode(),mtime=0))
 summary={'status':'complete_initial_level5_only','all_full_band_readbacks_pass':True,'full_band_readback_count':sum(r['positive_band_count'] for r in readbacks),'cached_endpoint_records_independently_checked':sum(r['cached_endpoint_records_checked'] for r in readbacks),'full_band_readback_wall_sum':sum(r['seconds'] for r in readbacks),'cell_count':52360,'weighted_volume':'1','sum_weights':1048576,'H_integral_lower':str(LH),'H_integral_upper':str(UH),'H_integral_lower_simple':str(F(floor(LH*10**12),10**12)),'H_integral_upper_simple':str(F(ceil(UH*10**12),10**12)),'I_upper_simple':str(F(ceil(F(ep(I,True))*10**12),10**12)),'H_integral_lower_display':float(LH),'H_integral_upper_display':float(UH),'I_upper':ep(I,True),'I_upper_display':float(I),'E076_norm_upper':str(E),'UH_sufficient_budget_for_I_le_1over25000_lower':ep(UHbudget),'UH_sufficient_budget_display':float(UHbudget),'UH_above_sufficient_budget':bool(arb(rq(UH))>UHbudget),'target_I':'1/25000','new_global_g':False,'new_Cn_claim':False,'no_tube_added':True,'new_H_count':sum(r['new_H_count'] for r in receipts),'reused_H1over64_count':sum(r['reused_H_count'] for r in receipts),'new_jet_count':0,'new_cone_count':sum(r['runtime']['new_cone_count'] for r in rows),'production_program_count':len(receipts),'production_program_walls':[r['wall_seconds'] for r in receipts],'production_wall_sum':sum(r['wall_seconds'] for r in receipts),'production_child_user_seconds':sum(r['children_user_seconds'] for r in receipts),'production_child_system_seconds':sum(r['children_system_seconds'] for r in receipts),'max_program_wall':max(r['wall_seconds'] for r in receipts),'all_workers_joined':all(r['workers_joined'] for r in receipts),'global_cap_used_count':sum(r['bounds']['global_cap_used'] for r in rows),'gate_1p05_pass_count':sum(F(r['bounds']['relative_baseline_upper_i'],N)<=F(21,20) for r in rows),'gate_1p1_pass_count':sum(F(r['bounds']['relative_baseline_upper_i'],N)<=F(11,10) for r in rows),'worst32_absolute_gap_ranks':[r['rank'] for r in worst],'worst8_relative_baseline_ranks':[r['rank'] for r in ratio],'evidence_tree_sha256':hashlib.sha256(''.join(r['evidence_sha256'] for r in rows).encode()).hexdigest(),'endpoint_sample_count':len(samples),'endpoint_samples_sha256':sha(H/'endpoint_samples.json.gz'),'compact_cells_sha256':sha(H/'initial_cells.json.gz'),'state_sha256':sha(S/'state.json'),'seconds':time.monotonic()-start,'pid':os.getpid(),'source_sha256':sha(__file__),'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (H/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');assert summary['new_H_count']==72*52360 and summary['reused_H1over64_count']==52360;print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
