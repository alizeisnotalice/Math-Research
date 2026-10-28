import os,time,resource,json,hashlib,collections
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
mp.mp.dps=100;st=time.monotonic();cpu=time.process_time();src=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round070/compute');here=Path(__file__).resolve().parent;sha=lambda b:hashlib.sha256(b).hexdigest()
configs=[('pilot','pilot_summary.json','registration.json','h2.py',set(range(14))),('production','production_summary.json','registration.json','h2.py',set(range(14,256))),('width_surrogate','width_summary.json','width_registration.json','finite_width.py',set(range(256))),('width_natural','width_e_summary.json','width_e_registration.json','finite_width_e.py',set(range(256))),('frontier4','frontier4_summary.json','frontier4_registration.json','frontier4.py',{(i,j) for i in range(20) for j in range(64)})]
results=[]
for name,sn,rn,code,indices in configs:
 sb=(src/sn).read_bytes();rb=(src/rn).read_bytes();s=json.loads(sb);r=json.loads(rb);assert s['registration_sha256']==sha(rb) and r['source_sha256']==sha((src/code).read_bytes());seen=set();cpu_sum=0;rss=0;statuses=collections.Counter();counts=collections.Counter();pids=set()
 for item in s['receipts']:
  raw=(src/item['file']).read_bytes();assert sha(raw)==item['sha256'];d=json.loads(raw);idx=(d['il'],d['iq']) if name=='frontier4' else d['i'];assert idx not in seen;seen.add(idx);assert d['status']==item['status'] and d['source_sha256']==r['source_sha256'];statuses[d['status']]+=1;pids.add(d['pid']);cpu_sum+=d['cpu_seconds'];rss=max(rss,d['peak_rss_bytes'])
  if name=='frontier4':assert d['counts']==item['counts'] and sum(d['counts'].values())==64*32*4;counts.update(d['counts'])
  elif name.startswith('width'):assert len(d['rows'])==item['states']==16;counts['states']+=16
  else:assert d['states']==item['states']==11250 and len(d['independent_checks'])==item['checks']==5;counts['states']+=d['states'];counts['crosschecks']+=len(d['independent_checks'])
 assert seen==indices and len(seen)==len(s['receipts']);monfile=src/({'width_surrogate':'width_monitor.jsonl','width_natural':'width_e_monitor.jsonl'}.get(name,name+'_monitor.jsonl'));mon=[json.loads(x) for x in monfile.read_text().splitlines()] if monfile.exists() else [];good=[x for x in mon if 'worker_count' in x];results.append({'stage':name,'summary_sha256':sha(sb),'registration_sha256':sha(rb),'tasks':len(seen),'statuses':dict(statuses),'counts':dict(counts),'recorded_worker_cpu_seconds':cpu_sum,'max_worker_peak_rss_bytes':rss,'distinct_pid_count':len(pids),'reported_wall_seconds':s['wall_seconds'],'monitor_samples':len(mon),'successful_monitor_samples':len(good),'maximum_monitor_cpu_percent':max((x.get('total_cpu_percent',x.get('cpu_percent',0)) for x in good),default=None)})
# Independent direct integration at one rational point: normalization and source-domain split.
q=F(9,10);s=F(3,2);h=F(1,1000);eps=F(1,10000);rho=F(27,5);w=F(91,100);R=s*q-h;a=(1-eps)*q;g=R-a;assert 0<g<R and abs(w)<=s*q
v=lambda x:mp.mpf(x.numerator)/x.denominator
rh,rm,am,gm,wm=map(v,(rho,R,a,g,w));cuts=sorted({-rh,rh,-rh+gm,rh-gm,max(-rh,min(rh,wm-rh-am)),max(-rh,min(rh,wm+rh+am))})
def f(z):
 ell=min(mp.mpf(1),(rm+rh+am-abs(z))/(2*rm));return (-mp.log(ell))**2/ell if abs(z-wm)<rh+am else mp.mpf(0)
direct=mp.fsum(mp.quad(f,[u,vv]) for u,vv in zip(cuts,cuts[1:]))/(2*rh);m=max(F(0),abs(w)-a);Q0=-mp.log1p(-gm/(2*rm));Q1=-mp.log1p(-v(g-m)/(2*rm)) if g>m else mp.mpf(0);closed=rm*(Q0**3+Q1**3)/(3*rh);err=abs(direct-closed);assert err<mp.mpf('1e-90')
out={'scope':'receipts, resource metadata, and one independent numerical H2 identity; no operator theorem or exhaustive arithmetic replay','stages':results,'single_H2_test':{'q':str(q),'s':str(s),'h':str(h),'eps':str(eps),'rho':str(rho),'w':str(w),'direct':mp.nstr(direct,60),'closed':mp.nstr(closed,60),'absolute_error':mp.nstr(err,10),'precision_decimal':100},'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'cpu_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-st,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(here/'check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
