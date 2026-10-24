"""Read completed global checker outputs; no model or certificate rerun."""
from pathlib import Path
from hashlib import sha256
import json,re
root=Path(__file__).resolve().parent;out=root/'full_checker_summary.json'
if out.exists():raise SystemExit('Summary already exists.')
summary={'scope':'Sufficient enclosure failures are not negative point proofs. Both implementations share Arb and analytic formulas.'}
for stage in ('first','second'):
 run=json.loads((root/(stage+'.run.json')).read_text());assert run['status'] in ('terminal','timed_out')
 item={'run':run};path=root/(stage+'_report.json')
 if path.exists():
  raw=path.read_bytes();d=json.loads(raw);item.update(report_sha256=sha256(raw).hexdigest(),report_bytes=len(raw))
  assert item['report_sha256']==run['report_sha256']
  if stage=='first':
   item['summary']={k:d[k] for k in ('mass_half','U','A_cells')}
   for kind in ('K','A'):
    fails=d[kind+'_fail'];item['summary'][kind+'_fail_count']=len(fails)
    item['summary'][kind+'_first_failures']=fails[:3]
    if kind=='A':item['summary'][kind+'_lowest_reported_failures']=sorted(fails,key=lambda x:float(x[-1]))[:3]
  else:
   assert d['spec_sha256']==run['candidate_sha256']
   item['summary']={k:d[k] for k in ('spec_sha256','claim','provenance','args','c1','c2','c3','sec')}
   c4=d['c4'];item['summary']['c4']={k:v for k,v in c4.items() if k!='results'}
   failed=[x for x in c4['results'] if not x['ok']]
   item['summary']['c4_nfail']=sum(x.get('nfail',len(x.get('fails',[]))) for x in c4['results'])
   item['summary']['failed_parent_count']=len(failed)
   item['summary']['failed_parent_counts']=[{'p':x['p'],'q':x['q'],'nfail':x.get('nfail'),'maxdepth':x.get('maxdepth')} for x in failed]
   item['summary']['first_failed_parents']=failed[:2]
   match=re.search(r'C1 K >= 0:.*failures (\d+)',(root/'second.log').read_text());assert match
   item['summary']['c1_total_failures_from_log']=int(match.group(1))
   item['summary']['c1_saved_witness_count']=len(d['c1']['fails']);item['summary']['c1_witness_cap']=50
   item['summary']['three_constrained_parents']=[{k:v for k,v in x.items() if k!='fails'} for x in c4['results'] if (x['p'],x['q']) in ((7,5),(8,6),(10,6))]
 else:item['report_missing']=True
 summary[stage]=item
out.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:{'exit_code':v['run']['exit_code'],'status':v['run']['status'],'report_missing':v.get('report_missing',False)} for k,v in summary.items() if isinstance(v,dict)}))
