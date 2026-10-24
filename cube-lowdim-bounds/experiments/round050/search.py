from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import sys,json,subprocess,time,os,signal,copy
seed,checker=map(Path,sys.argv[1:3]);out=Path(sys.argv[3]);out.mkdir(exist_ok=True)
assert not (out/'result.json').exists()
assert sha256(seed.read_bytes()).hexdigest()=='e2ca1616be95afa6680e3ccf157cc16fdd39026ad0f5fa79ebd32e5214779888'
assert sha256(checker.read_bytes()).hexdigest()=='abb311c73157e35f384ee2aca6aded4772cc13ea8f70650b8499400621a80e23'
c=json.loads(seed.read_text());px=F(c['Px']);py=F(c['Py']);assert len(c['pts'])==50
old=F(482167717377842991882283179477574801392641890014291050323,283174368139841919615520431313000857449005056000000000000)
target=F(33977370204434849545366297954842256187719839393,19950077409778055009062309259909267456000000000)
start=time.monotonic();rows=[];status='completed';best=None
for h in map(F,['-1/500','-1/1000','1/1000','1/500']):
 for m in map(F,['-1/50','-1/100','1/100','1/50']):
  remain=120-(time.monotonic()-start)
  if remain<=0:status='batch_timeout';break
  v=copy.deepcopy(c)
  for a in v['pts'][2:]:a[0]=str(px/2+(1+h)*(F(a[0])-px/2));a[2]=str(F(a[2])*(1+m))
  assert all(0<=F(x)<px and 0<=F(y)<py and F(w)>0 for x,y,w in v['pts'])
  n=len(rows);p=out/('candidate%02d.json'%n);p.write_text(json.dumps(v,indent=1)+'\n');log=out/('candidate%02d.log'%n)
  t=time.monotonic();r={'index':n,'h':str(h),'m':str(m),'config_sha256':sha256(p.read_bytes()).hexdigest()}
  with log.open('x') as f:
   proc=subprocess.Popen([sys.executable,'-I',str(checker),str(p),'--Dr',str(2**40)],stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
   try:r['exit_code']=proc.wait(timeout=min(10,remain))
   except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait();r['exit_code']=proc.returncode;status='timeout'
  r['seconds']=time.monotonic()-t
  if r['exit_code']!=0:status='failed_or_timeout';rows.append(r);break
  a=json.loads(p.with_suffix('.check2.json').read_text());q=F(int(a['ratio_num']),int(a['ratio_den']));r.update(ratio=str(q),ratio_display=float(q),rectangles=a['rectangles'],above_public=q>old,above_claimed58=q>target)
  rows.append(r)
  if best is None or q>F(best['ratio']):best=r
 if status!='completed':break
assert sha256(checker.read_bytes()).hexdigest()=='abb311c73157e35f384ee2aca6aded4772cc13ea8f70650b8499400621a80e23'
result={'status':status,'budget_seconds':120,'per_call_seconds':10,'elapsed_seconds':time.monotonic()-start,'completed':sum('ratio' in r for r in rows),'rows':rows,'best':best,'public':str(old),'claimed58_unverified':str(target),'independent_certification_required':bool(best and best['above_claimed58'])}
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
