from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,sys
p=Path(__file__).resolve().parent
reg=json.loads((p/'registration.json').read_text());rows=[]
cfile=p.parent/'lower_round003/a2_equals_a3_results.json'
clo,chi=map(Q,json.loads(cfile.read_text())['ratio_interval'])
source=hashlib.sha256((p/'parity.py').read_bytes()).hexdigest()
for spec in reg['candidates']:
 name=spec['name'];inn=json.loads((p/(name+'_results.json')).read_text());out=json.loads((p/(name+'_outer_results.json')).read_text())
 assert inn['source_sha256']==out['source_sha256']==source
 assert Q(inn['ratio'])<=Q(out['ratio'])<clo<chi
 assert inn['union_engines_agree'] and out['union_engines_agree']
 for mode,result in [('inner',inn),('outer',out)]:
  for w in result['witnesses']:
   r,m=Q(w['radius']),Q(w['mass'])
   assert (32*r**5<m) if mode=='inner' else (32*r**5>=m)
   assert w['distinct_positive_atoms']>0
  assert 8*Q(result['volume'])==Q(result['ratio'])
 rows.append({'name':name,'ratio_inner':inn['ratio'],'ratio_outer':out['ratio'],'exact_positive_gap':str(clo-Q(out['ratio']))})
receipt={'status':'VERIFIED_ALL_FOUR_FINITE_FAMILY_EXCLUSIONS','python':sys.executable,'c_source_sha256':hashlib.sha256(cfile.read_bytes()).hexdigest(),'source_sha256':source,'candidate_count':len(rows),'rows':rows,'files':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(p.iterdir()) if f.is_file() and f.name!='completion_receipt.json'}}
(p/'completion_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(receipt['status']);print('All exact inner <= outer < C_lower comparisons hold.')
