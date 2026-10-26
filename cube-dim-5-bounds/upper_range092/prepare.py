from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib,time
import ratio091 as R
H=Path(__file__).resolve().parent;OLD=H.parent/'upper_reenvelope088';S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_range092')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def gz(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(gzip.compress(serial(x).encode(),mtime=0))
def main():
 t=time.monotonic();assert not (H/'nodes.json.gz').exists();op=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_reenvelope088/state.json.gz');assert sha(op)==json.loads((OLD/'summary.json').read_text())['state_sha256'];st=json.loads(gzip.decompress(op.read_bytes()));assert len(st['nodes'])==56139 and len(st['frontier'])==55819;srcnodes=json.loads(gzip.decompress((OLD/'nodes.json.gz').read_bytes()));source={r['key']:r for r in srcnodes};certs={}
 for p in sorted(OLD.glob('block*_cells.json.gz')):
  for r in json.loads(gzip.decompress(p.read_bytes())):certs[r['key']]=r
 ts=json.loads((OLD/'registration.json').read_text())['time_endpoints'];xs=sorted({F(r['h_infinity'])/F(a) for r in certs.values() for a in ts});ratio_reg={'id':'upper_range092.ratio_catalog','source091_sha256':sha(H/'ratio091.py'),'source088_manifest_sha256':sha(OLD/'manifest.json'),'x_values':list(map(str,xs)),'bits':192,'H_calls':0,'old_ratio_calls':0};(H/'ratio_registration.json').write_text(json.dumps(ratio_reg,indent=2)+'\n');catalog={}
 for x in xs:
  lo,hi=R.ratio_interval(x,192);s1l,s1u=R._sqrt_interval(1+x*x/16,192);s2l,s2u=R._sqrt_interval(1+x*x/256,192);catalog[str(x)]={'lower':str(lo),'upper':str(hi),'sqrt1_lower':str(s1l),'sqrt1_upper':str(s1u),'sqrt2_lower':str(s2l),'sqrt2_upper':str(s2u)}
 (H/'ratio_catalog.json').write_text(serial(catalog)+'\n');print('CATALOG',len(catalog),flush=True);front=set(st['frontier']);order=sorted(st['nodes'],key=lambda k:(k not in front,-(st['nodes'][k]['selfU']-st['nodes'][k]['selfL']),st['nodes'][k]['level'],st['nodes'][k]['indices']));rows=[]
 for rank,k in enumerate(order):
  n=st['nodes'][k];c=certs[k];src=source[k];assert n['selfL']==c['old_self_L'] and n['selfU']==c['new_self_U'];raw=Path(src['source_evidence_path']);assert sha(raw)==src['source_evidence_sha256'];ev=json.loads(gzip.decompress(raw.read_bytes()));assert ev['indices']==n['indices'];ends=ev['endpoints'];assert [x['t'] for x in ends]==ts;packet={'key':k,'point':ev['point'] if 'point' in ev else ev['anchor'],'endpoints':ends,'original_raw_sha256':src['source_evidence_sha256']};pp=S/'bands'/f'{rank//1024:02d}'/f'{rank:05d}.json.gz';gz(pp,packet)
  rows.append({'rank':rank,'key':k,'level':n['level'],'indices':n['indices'],'weight':n['weight'],'old_self_L':n['selfL'],'old_self_U':n['selfU'],'old_rangeL_i':c['rangeL_i'],'old_rangeU_i':c['rangeU_i'],'old_M_i':c['new_mean_i'],'M2':c['M2'],'mu_infinity':c['mu_infinity'],'h_infinity':c['h_infinity'],'physical_frontier':k in front,'packet_path':str(pp),'packet_sha256':sha(pp),'original_raw_path':str(raw),'original_raw_sha256':src['source_evidence_sha256'],'source088_rank':c['rank'],'source088_record_sha256':hashlib.sha256(serial(c).encode()).hexdigest()})
  if (rank+1)%8192==0:print('PACKED',rank+1,round(time.monotonic()-t,3),flush=True)
 gz(H/'nodes.json.gz',rows);gz(S/'initial_state.json.gz',st);receipt={'status':'complete','nodes':len(rows),'ratio_count':len(catalog),'initial_state_sha256':sha(S/'initial_state.json.gz'),'source088_state_sha256':sha(op),'nodes_sha256':sha(H/'nodes.json.gz'),'catalog_sha256':sha(H/'ratio_catalog.json'),'ratio_registration_sha256':sha(H/'ratio_registration.json'),'source_ratio_sha256':sha(H/'ratio091.py'),'endpoint_packets_external':True,'old_exp_asinh_H_calls':0,'seconds':time.monotonic()-t,'source_sha256':sha(__file__)};(H/'prepare_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
