#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json,sys,time,traceback,resource
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
lock=json.loads((H/'input_lock.json').read_text());reg=json.loads((H/'registration.json').read_text())
for n,h in lock['files'].items():assert sha(H/n)==h,n
old=H.parent/'round003/candidate_left.json';assert sha(old)==reg['baseline_r003_left_sha256'];baseline=json.loads(old.read_text())['result'];ref=F(baseline['union_interval'][1]);inherit=F(baseline['finite_L1']['lower'])
s=importlib.util.spec_from_file_location('geom',H/'geometry.py');g=importlib.util.module_from_spec(s);s.loader.exec_module(g)
rank=int(sys.argv[1]);event=next(e for e in reg['selected'] if e['rank']==rank);a=event['selected_alias'];terms=[(a['p'],F(1,2)),(a['q'],F(1,2))] if a['kind']=='sum' else [(a['p'],F(1,2) if a['kind']=='r' else F(1))]
t=time.time();res={'event':event,'comparison_upper':str(ref),'inherited_finite_lower':str(inherit),'input_lock_sha256':sha(H/'input_lock.json')}
try:
 roots=[]
 for m,c in terms:
  q=g.root(m*g.D**8); exact=q**8==m*g.D**8; assert q**8<=m*g.D**8<(q+1)**8
  roots.append({'radicand':m,'coefficient':str(c),'lower_numerator':str(q),'upper_numerator':str(q if exact else q+1),'perfect_eighth_power':exact})
 low=4*sum(c*int(x['lower_numerator']) for (_,c),x in zip(terms,roots));high=4*sum(c*int(x['upper_numerator']) for (_,c),x in zip(terms,roots));assert low.denominator==high.denominator==1
 rr,tail,records,orbit=g.build();res['event_roots']=roots;res['result']=g.calculate(str(rank),(int(low),int(high)),rr,tail,records,orbit,ref,inherit);res['status']='complete'
except g.Exhausted:res['status']='memo limit reached; incomplete'
except Exception as e:res['status']='program failure; incomplete';res['error']=repr(e);res['traceback']=traceback.format_exc()
res['elapsed_seconds']=time.time()-t;res['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
with (H/f'candidate_{rank:02}.json').open('x') as f:json.dump(res,f,indent=2);f.write('\n')
print(json.dumps({'rank':rank,'status':res['status'],'winner':res.get('result',{}).get('winner',False),'peak_rss_bytes':res['peak_rss_bytes']}),flush=True)
