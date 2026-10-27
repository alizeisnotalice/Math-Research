#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json,sys,time,traceback
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
lock=json.loads((H/'input_lock.json').read_text()); reg=json.loads((H/'registration.json').read_text())
for name in ('registration.json','geometry.py','candidate.py','theory_event_certificate.json','theory_manifest.json'):
 assert sha(H/name)==lock['local_sha256'][name],name
s=importlib.util.spec_from_file_location('geometry',H/'geometry.py');g=importlib.util.module_from_spec(s);s.loader.exec_module(g)
name=sys.argv[1];event=next(e for e in reg['events'] if e['name']==name)
oldpath=H.parent/'round002/candidate_96.json';assert sha(oldpath)==lock['r002_candidate_96_sha256']
old=json.loads(oldpath.read_text())['candidate'];ref=F(old['union_interval'][1]);inherit=F(old['finite_L1']['lower'])
t=time.time();res={'event':event,'comparison_upper':str(ref),'inherited_finite_lower':str(inherit),'input_lock_sha256':sha(H/'input_lock.json')}
try:
 rr,tail,records,orbit=g.build()
 q=[g.root(v*g.D**8) for v in event['roots']]
 for v,x in zip(event['roots'],q): assert x**8<=v*g.D**8<(x+1)**8
 spacing=(2*sum(q),2*(sum(q)+len(q)))
 res['event_root_lower_numerators']=list(map(str,q));res['result']=g.calculate(name,spacing,rr,tail,records,orbit,ref,inherit);res['status']='complete'
except g.Exhausted:res['status']='memo limit reached; incomplete'
except Exception as exc:res['status']='program failure; incomplete';res['error']=repr(exc);res['traceback']=traceback.format_exc()
res['wall_seconds_total']=time.time()-t
with (H/f'candidate_{name}.json').open('x') as f:json.dump(res,f,indent=2);f.write('\n')
print(json.dumps({'name':name,'status':res['status'],'winner':res.get('result',{}).get('winner',False)}),flush=True)
