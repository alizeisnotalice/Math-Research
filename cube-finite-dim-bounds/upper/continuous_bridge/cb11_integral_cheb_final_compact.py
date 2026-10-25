"""Lossless post-run compaction. No evaluator calls and no raw-file modifications."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,shutil,subprocess,sys
OUT=Path(__file__).resolve().parent
STAGE=OUT/'cb11_integral_cheb_final_staging'
BACKUP=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/round011o-run/bridge-raw-archive')
TOP=['config','checker_sha256','part','parts','patches']
KEYS=['patch','L','U','mean_H_upper','IH_contribution_upper','IH_contribution_exact']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def copy(p,q):
 if q.exists():assert sha(p)==sha(q)
 else:shutil.copy2(p,q)
def fraction_from_pack(a):
 assert int(a[2])==0
 return F(int(a[0]))*F(2)**int(a[1])
def run():
 assert not STAGE.exists(),'refuse to overwrite final staging'
 STAGE.mkdir();BACKUP.mkdir(parents=True,exist_ok=True)
 historical_manifest_sha=sha(OUT/'cb11_integral_cheb_manifest.json');checker_sha=sha(OUT/'cb11_integral_cheb_check.py')
 provenance={'status':'post-run lossless compaction; original bytes preserved','raw_archive':str(BACKUP),'raw_parts':{},'compact_parts':{},'historical_pilot_manifest_sha256':historical_manifest_sha,'checker_sha256':checker_sha,'patches':8256,'no_new_node_evaluations':True,'parent_execution':{'session':80449,'parts':12,'all_exit_codes':0,'log_directory':'/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/round011o-run'},'part_statistics':[]}
 raw_sum=F(0);compact_sum=F(0);seen=set();cachehits=[];raw_bytes=0;compact_bytes=0
 for part in range(12):
  name=f'cb11_integral_cheb_part_{part:02}.json';src=OUT/name;r=json.loads(src.read_text());assert r['part']==part and r['parts']==12 and r['checker_sha256']==checker_sha
  copy(src,BACKUP/name);provenance['raw_parts'][name]={'source':str(src),'archive':str(BACKUP/name),'sha256':sha(src),'bytes':src.stat().st_size};raw_bytes+=src.stat().st_size
  cr={k:r[k] for k in TOP if k!='patches'};cr['patches']=[]
  expected={(i,j) for i in range(128) for j in range(i+1) if (i*(i+1)//2+j)%12==part};assert {tuple(x['patch']) for x in r['patches']}==expected and len(r['patches'])==688
  for x in r['patches']:
   y={k:x[k] for k in KEYS};assert all(y[k]==x[k] for k in KEYS);cr['patches'].append(y)
   p=tuple(x['patch']);assert p not in seen;seen.add(p)
   v=F(x['IH_contribution_exact']);assert v==fraction_from_pack(x['IH_contribution_upper']);raw_sum+=v;compact_sum+=F(y['IH_contribution_exact'])
   if x['cache_hit']:cachehits.append(p)
  target=STAGE/name;save(target,cr)
  # Reparse serialized compact evidence and require every retained field unchanged.
  rr=json.loads(target.read_text());assert all(rr[k]==r[k] for k in TOP if k!='patches') and rr['patches']==[{k:x[k] for k in KEYS} for x in r['patches']]
  provenance['compact_parts'][name]={'sha256':sha(target),'bytes':target.stat().st_size};compact_bytes+=target.stat().st_size
  provenance['part_statistics'].append({'part':part,'patches':len(r['patches']),'wall_seconds':r['wall_seconds'],'node_evaluations':sum(x['nodes'] for x in r['patches'] if not x['cache_hit']),'cached_patches':sum(bool(x['cache_hit']) for x in r['patches'])})
 assert len(seen)==8256 and seen=={(i,j) for i in range(128) for j in range(i+1)} and raw_sum==compact_sum
 assert set(cachehits)=={(127,57),(127,127),(127,0),(0,0)} and len(cachehits)==4
 aggregate=json.loads((OUT/'cb11_integral_cheb_all_result.json').read_text());assert F(aggregate['IH_upper_exact'])==raw_sum
 for name in ['cb11_integral_cheb_check.py','cb11_integral_cheb_selfcheck.py','cb11_cones.json']:copy(OUT/name,STAGE/name)
 for name in ['cb11_integral_cheb_all_result.json','cb11_integral_cheb_all_selfcheck.json','cb11_integral_cheb_four_cache.json']:
  copy(OUT/name,BACKUP/name)
 provenance['optional_four_cache']={'archive':str(BACKUP/'cb11_integral_cheb_four_cache.json'),'sha256':sha(OUT/'cb11_integral_cheb_four_cache.json'),'used_cached_patches':4,'repository_required':False}
 for script,args in [('cb11_integral_cheb_selfcheck.py',['--parts','12']),('cb11_integral_cheb_check.py',['--stage','aggregate','--parts','12'])]:
  p=subprocess.run([sys.executable,str(STAGE/script),*args],capture_output=True,text=True,check=True)
  (STAGE/(script+'.stdout.txt')).write_text(p.stdout)
 compact_aggregate=json.loads((STAGE/'cb11_integral_cheb_all_result.json').read_text())
 numeric_fields=['status','I_certified_locally','config','patch_count','IH_upper_exact','IH_upper','Etime','I_upper','I_upper_display','triangle_weighted_area_exact','checker_sha256']
 assert all(compact_aggregate[k]==aggregate[k] for k in numeric_fields)
 for part in range(12):assert sha(OUT/f'cb11_integral_cheb_part_{part:02}.json')==provenance['raw_parts'][f'cb11_integral_cheb_part_{part:02}.json']['sha256']
 assert sha(OUT/'cb11_integral_cheb_manifest.json')==historical_manifest_sha and sha(OUT/'cb11_integral_cheb_check.py')==checker_sha
 provenance.update({'raw_bytes':raw_bytes,'compact_bytes':compact_bytes,'raw_exact_sum_equals_compact':True,'IH_upper_exact':str(raw_sum),'packed_final_I_identical':True,'compact_selfcheck':'PASS8256','compact_aggregate':'PASS exact numerical fields identical','full_nodes_total_including_four_pilot':8256*1089,'new_nodes_in_full_parent_batch':(8256-4)*1089})
 save(OUT/'cb11_integral_cheb_final_provenance.json',provenance)
 print(json.dumps({k:v for k,v in provenance.items() if k not in ['raw_parts','compact_parts','part_statistics']},indent=2))
if __name__=='__main__':run()
