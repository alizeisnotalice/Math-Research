#!/usr/bin/env python3
import json,hashlib,time,resource,os,sys
from fractions import Fraction as F
from pathlib import Path
sys.set_int_max_str_digits(0)
st=time.monotonic();cpu=time.process_time();base=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper');source=base/'round039/compute';old=base/'round038/compute';prior=Path('/Users/zhengzhihao/Desktop/T/cube-r030-num-audit/cube-dim-7-bounds/audit/round038/phasea/results');proof=Path('/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round039/joint_cover/PROOF.md')
def digest(raw):return hashlib.sha256(raw).hexdigest()
def Q(x):return F(x)
regraw=(source/'registration.json').read_bytes();reg=json.loads(regraw);assert digest(regraw)=='cc24240275ff240cad9fe18665fb817554de321a3c019cf6270ee803a910a639';assert digest(proof.read_bytes())=='b0719f3c04f0fb88fe794763e36e916636be58a0572db139db67ffb6e979c7d1'
rows=[]
files=sorted(source.glob('scan_*.json'));assert len(files)==6
for f in files:
 raw=f.read_bytes();r=json.loads(raw);case=r['case'];caseraw=(old/case).read_bytes();d=json.loads(caseraw);prevraw=(prior/(Path(case).stem+'.audit.json')).read_bytes();prev=json.loads(prevraw)
 assert prev['status']=='pass' and digest(caseraw)==prev['input_sha256']==r['input_sha256']==reg['input_hashes'][case]
 C=[tuple(c) for c in d['maximal_cliques']];h=list(map(Q,d['witness_budgets']));q=list(map(Q,d['q_after']));N=len(h);assert len(C)==prev['cliques'];i=r['best_pair']['i'];j=r['best_pair']['j']; assert 0<=i<j<N
 adjacent=any(i in c and j in c for c in C);assert r['domain']!='nonedge' or not adjacent
 w=[sum((q[k] for k in c),F()) for c in C];L=max(w);S=sum(q,F());E=sum((h[k]*q[k]**2 for k in range(N)),F());U=Q(d['U']);assert sum((h[k]*q[k] for k in range(N)),F())==U
 b=r['best_pair'];P=h[i]*q[i]+h[j]*q[j];assert P>0 and P==Q(b['P']);groups={}
 for c,wt in zip(C,w):
  vi=int(i in c);vj=int(j in c);tag=f'{vi}{vj}';v=wt-vi*q[i]-vj*q[j];groups[tag]=max(groups.get(tag,v),v)
 assert groups=={k:Q(v) for k,v in b['class_max_outside'].items()}
 lines={tag:(F(int(tag[0]))-int(tag[1])*h[i]/h[j],outside+int(tag[1])*P/h[j]) for tag,outside in groups.items()}
 assert lines=={x['class']:(Q(x['slope']),Q(x['intercept'])) for x in b['branches']}
 xmax=P/h[i];pts={F(0),xmax};lin=list(lines.values())
 for k,(a,z) in enumerate(lin):
  for aa,zz in lin[:k]:
   if a!=aa:
    x=(zz-z)/(a-aa)
    if 0<=x<=xmax:pts.add(x)
 assert pts==set(map(Q,b['candidates']))
 val=lambda x:max(a*x+z for a,z in lin)
 m=min(map(val,pts));assert m==Q(b['min_load'])
 lo,hi=F(0),xmax
 for a,z in lin:
  if a>0:hi=min(hi,(m-z)/a)
  elif a<0:lo=max(lo,(m-z)/a)
  else:assert z<=m
 assert 0<=lo<=hi<=xmax and [lo,hi]==list(map(Q,b['optimal_interval']))
 if h[i]>h[j]:x=hi
 elif h[i]<h[j]:x=lo
 else:x=max(lo,min(hi,P/(h[i]+h[j])))
 y=(P-h[i]*x)/h[j];assert x==Q(b['new_i']) and y==Q(b['new_j']) and min(x,y)>=0
 after=q.copy();after[i]=x;after[j]=y;assert after==list(map(Q,r['q_after'])) and sum((h[k]*after[k] for k in range(N)),F())==U
 L2=max(sum((after[k] for k in c),F()) for c in C);S2=sum(after,F());E2=sum((h[k]*after[k]**2 for k in range(N)),F())
 assert L==Q(r['L_before']) and L2==Q(r['L_after'])==m==L and S==Q(r['S_before']) and S2==Q(r['S_after']) and S2<S and E==Q(r['E_before']) and E2==Q(r['E_after']) and E2>E
 assert S2-S==Q(b['deltaS']) and E2-E==Q(b['deltaE']) and r['change_class']=='sum_q' and r['no_iteration']
 rows.append({'case':case,'domain':r['domain'],'pair':[i,j],'adjacent':adjacent,'case_sha256':digest(caseraw),'scan_sha256':digest(raw),'prior_audit_sha256':digest(prevraw),'class_count':len(groups),'breakpoint_count':len(pts),'checks':'pass','relative_sum_q_reduction':float((S-S2)/S),'relative_energy_increase':float((E2-E)/E)})
out={'scope':'six chosen pairs, not full-pair scan or global-best verification','registration_sha256':digest(regraw),'proof_sha256':digest(proof.read_bytes()),'audit_script_sha256':digest(Path(__file__).read_bytes()),'rows':rows,'pid':os.getpid(),'wall_seconds':time.monotonic()-st,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('audit_six.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
