#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction
import hashlib,importlib.util,json,sys,time,traceback
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
s=importlib.util.spec_from_file_location('geometry',H/'geometry_r001.py');g=importlib.util.module_from_spec(s);s.loader.exec_module(g)
reg=json.loads((H/'registration.json').read_text());m=int(sys.argv[1]);assert m in reg['candidates']
oldpath=H.parent/'round001/results.json'; lock=json.loads((H/'input_lock.json').read_text()); assert sha(H/'registration.json')==lock['registration_sha256'] and sha(H/'geometry_r001.py')==lock['geometry_r001_sha256'] and sha(oldpath)==lock['r001_results_sha256']; old=json.loads(oldpath.read_text()); ref=Fraction(old['trials'][0]['union_interval'][1]); inherited=Fraction(old['trials'][0]['finite_L1']['lower'])
started=time.time();result={'m':m,'registration_sha256':sha(H/'registration.json'),'geometry_sha256':sha(H/'geometry_r001.py'),'r001_results_sha256':sha(oldpath),'comparison_upper':str(ref),'inherited_lower':str(inherited)}
try:
 rr,tail,records,orbit=g.build(); result['candidate']=g.calculate(m,rr,tail,records,orbit,ref,inherited);result['status']='complete'
except g.Exhausted:result['status']='memo limit reached; incomplete'
except Exception as e:result['status']='program failure; incomplete';result['error']=repr(e);result['traceback']=traceback.format_exc()
result['wall_seconds_total']=time.time()-started
(H/f'candidate_{m}.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'m':m,'status':result['status'],'winner':result.get('candidate',{}).get('winner',False),'wall_seconds_total':result['wall_seconds_total']},flush=True))
