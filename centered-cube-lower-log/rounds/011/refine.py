from pathlib import Path
import json,numpy as np
P=Path(__file__).resolve().parent;old=json.loads((P/'results.json').read_text())['rows'];n=2048;R=11;N=32;p=.9
rng=np.random.default_rng(2026102807);u=rng.random((1024,n));records=[]
for name,rs,ws in [('single',[R],[1]),('adjacent',[R,R+1],[.5,.5])]:
 h=next(r['log_threshold'] for r in old if r['n']==n and r['name']==name)
 freq={}
 for r,w in zip(rs,ws):
  for j in range(1,N):freq[j*r]=freq.get(j*r,0)+2*p*w*(1-j/N)
 m=np.array(sorted(freq));coef=np.array([freq[k] for k in m]);vals=[]
 for G,T in [(2048,257),(8192,1025)]:
  C=np.zeros((len(u),G));inds=(u*G).astype(int);np.add.at(C,(np.repeat(np.arange(len(u)),n),inds.ravel()),1)
  x=(np.arange(G)+.5)/G;L=np.linspace(1.02,1.48,T)
  avg=1+(np.sin(np.pi*L[:,None]*m)/(np.pi*L[:,None]*m)*coef)@np.cos(2*np.pi*m[:,None]*x)
  v=(C@np.log(avg).T).max(axis=1);vals.append(v);q=np.mean(v>h)
  records.append(dict(name=name,center_grid=G,scale_grid=T,frozen_h=h,hits=int(np.sum(v>h)),estimate=float(np.exp(h)*q),se=float(np.exp(h)*np.sqrt(q*(1-q)/len(v)))))
 records[-1]['paired_max_log_response_change']=float(np.max(abs(vals[1]-vals[0])))
 records[-1]['paired_threshold_disagreements']=int(np.sum((vals[1]>h)!=(vals[0]>h)))
(P/'refinement.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
