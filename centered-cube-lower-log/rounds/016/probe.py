from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent
src=P.parent/'014/results.json'
c=np.array(json.loads(src.read_text())['witness']['counts'],dtype=np.longdouble)
f=np.full(64,16,dtype=np.longdouble);f[0]-=1
labels=np.array([[i//16,(i//4)%4,i%4] for i in range(64)],dtype=np.longdouble)
assert c.sum()==f.sum()==1023
S=((labels-1.5)/np.sqrt(np.longdouble(1.25))).sum(axis=1);S-=np.dot(f,S)/1023
assert abs(np.dot(c-f,S))<1e-12 and abs(np.dot(c-f,S*S))<1e-12
rows=[]
for n in [10000,1000000,100000000]:
 v=np.longdouble(16)*np.log(np.longdouble(n));a=np.sqrt(v/n)
 for t in [.5,1,2]:
  z=np.expm1(t*a*S);diff=n*(np.log1p(np.dot(c,z)/1023)-np.log1p(np.dot(f,z)/1023))
  rows.append(dict(n=n,t=t,v=float(v),difference=float(diff),cubic_scale=float(n*a**3),ratio=float(diff/(n*a**3))))
(P/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
