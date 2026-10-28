"""Read-only independent log-domain winner/tie audit; never samples."""
from pathlib import Path
import math,json,hashlib,argparse
import numpy as np
HERE=Path(__file__).resolve().parent
CELL=HERE/'cell'
if not CELL.exists():CELL=HERE.parent/'round037_quantile_cell'
final=json.loads((CELL/'RESULTS_stable_ESS.json').read_text())
revision=json.loads((CELL/'ESS_revision.json').read_text())
for n,k in [('probe.py','original_source_SHA256'),('registration.json','registration_SHA256'),('RESULTS.json','original_RESULTS_SHA256'),('replay_stable_ess.py','fix_script_SHA256'),('ESS_FIX.patch','minimal_patch_SHA256')]:
 assert hashlib.sha256((CELL/n).read_bytes()).hexdigest()==revision[k]
def logadd(a,b):
 if a==-math.inf:return b
 if b==-math.inf:return a
 return max(a,b)+math.log1p(math.exp(-abs(a-b)))
def stable_scan(m,ld,ess_min):
 order=np.argsort(m)[::-1];ls=lss=-math.inf;best=-math.inf;answer=None
 for i in order:
  lw=-float(ld[i]);ls=logadd(ls,lw);lss=logadd(lss,2*lw)
  le=2*ls-lss;score=float(m[i])+ls-math.log(len(ld))
  if le>=math.log(ess_min)-1e-10 and score>best:
   best=score;answer=(math.exp(score),float(m[i]),math.exp(le))
 return answer
rows=[];mc=me=0.;tie_flags=0;count=0
for case,saved in zip(final['cases'],final['same_data_protocol_replay_checks']):
 path=CELL/saved['local_raw_logdata'];assert hashlib.sha256(path.read_bytes()).hexdigest()==saved['raw_logdata_SHA256']
 with np.load(path,allow_pickle=False)as data:m=data['logM'];ld=data['logD']
 for i,(name,vals)in enumerate(case['models'].items()):
  for key,miness in [('empirical_unrestricted',1),('empirical_ESS_ge64',64)]:
   target=vals[key];scan=stable_scan(m[:,i],ld,miness)
   assert abs(scan[0]-target['coefficient'])<1e-10
   mask=m[:,i]>=target['log_threshold'];lw=-ld[mask];top=float(lw.max())
   ls=top+math.log(math.fsum(math.exp(float(z)-top)for z in lw))
   lss=2*top+math.log(math.fsum(math.exp(2*(float(z)-top))for z in lw))
   coef=math.exp(target['log_threshold']+ls-math.log(len(ld)));ess=math.exp(2*ls-lss)
   dc=abs(coef-target['coefficient']);de=abs(ess-target['ESS']);mc=max(mc,dc);me=max(me,de)
   tie_flags+=int(dc>1e-10 or de>1e-7);count+=1
 assert tie_flags==0
 rows.append({'n':case['metadata']['n'],'R':case['metadata']['R'],'saved_log_hash_passed':True})
out={'scope':'Independent stdlib logaddexp scan and fsum complete-tie recomputation from saved data only; no MC or source evaluation; no FP interval certificate.',
 'numpy':np.__version__,'winner_count':count,'max_complete_tie_coefficient_difference':mc,'max_complete_tie_ESS_difference':me,'tie_omission_flags':tie_flags,'cases':rows}
parser=argparse.ArgumentParser();parser.add_argument('--output');args=parser.parse_args()
if args.output:
 p=Path(args.output)
 if p.exists():raise RuntimeError('Refusing to overwrite an audit output')
 p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='cases'},indent=2))
