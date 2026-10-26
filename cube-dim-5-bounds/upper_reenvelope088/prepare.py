from pathlib import Path
import gzip,json,hashlib,time
H=Path(__file__).resolve().parent;O=H.parent/'upper_global080';A=H.parent/'upper_adapt082';S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_reenvelope088')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def gz(p,obj):p.write_bytes(gzip.compress(json.dumps(obj,sort_keys=True,separators=(',',':')).encode(),mtime=0))
def main():
 t=time.monotonic();assert not (H/'nodes.json.gz').exists();stpath=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adapt082/state.json.gz');ast=json.loads(gzip.decompress(stpath.read_bytes()));assert sha(stpath)==json.loads((A/'summary.json').read_text())['state_sha256'];assert len(ast['nodes'])==56139 and len(ast['frontier'])==55819;old=json.loads(gzip.decompress((O/'initial_cells.json.gz').read_bytes()));rootcert={}
 for xs in old['rows']:
  r=dict(zip(old['columns'],xs));rootcert[r['rootid']]=r
 acerts={}
 for p in sorted(A.glob('batch*_cells.json.gz')):
  for r in json.loads(gzip.decompress(p.read_bytes())):acerts[r['key']]=r
 front=set(ast['frontier']);order=sorted(ast['nodes'],key=lambda k:(k not in front,-(ast['nodes'][k]['selfU']-ast['nodes'][k]['selfL']),ast['nodes'][k]['level'],ast['nodes'][k]['indices']));rows=[]
 for rank,k in enumerate(order):
  n=ast['nodes'][k]
  if n['source']=='080':
   c=rootcert[n['rootid']];rr=c['rank'];path=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_global080/evidence')/f'{rr//1024:02d}'/f'{rr:05d}.json.gz';expected=c['evidence_sha256'];assert n['selfU']==n['weight']*c['integral_upper_i'] and n['selfL']==n['weight']*c['integral_lower_i']
  else:
   c=acerts[k];path=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adapt082/evidence')/(k.replace(':','_').replace(',','-')+'.json.gz');expected=c['evidence_sha256'];assert n['selfU']==n['weight']*c['bounds']['integral_upper_i'] and n['selfL']==n['weight']*c['bounds']['integral_lower_i']
  rows.append({'rank':rank,'key':k,'level':n['level'],'indices':n['indices'],'weight':n['weight'],'source_evidence_path':str(path),'source_evidence_sha256':expected,'old_self_U':n['selfU'],'old_self_L':n['selfL'],'physical_frontier':k in front})
 assert len({r['key'] for r in rows})==56139;gz(H/'nodes.json.gz',rows);gz(S/'initial_state.json.gz',ast);(H/'prepare_receipt.json').write_text(json.dumps({'source082_state_sha256':sha(stpath),'copied_initial_state_sha256':sha(S/'initial_state.json.gz'),'nodes_sha256':sha(H/'nodes.json.gz'),'count':len(rows),'physical_frontier':len(front),'ordering':'physical leaves first, then descending exact old self weighted gap, level/indices; historical internal nodes last','seconds':time.monotonic()-t,'source_sha256':sha(__file__)},indent=2)+'\n');print((H/'prepare_receipt.json').read_text())
if __name__=='__main__':main()
