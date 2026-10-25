"""Only the four preregistered nonzero exact path parameters; full unions."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
import numpy as np
import round011_geometry as g
from round006_union import volume
R=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):(R/p).write_text(json.dumps(x,indent=2)+'\n')
def parameters(s,d):
 a=list(map(F,d['period_velocity']));v=[list(map(F,q)) for q in d['phase_velocity']];delta=F(d['delta'])
 return {'periods':[str(g.P[j]*(1+a[j]*s)) for j in range(3)],'phases':[[str(g.PH[b][j]+v[b][j]*s) for j in range(3)] for b in range(3)],'weights':[str(g.WTS[0]),str(g.WTS[1]+delta*s),str(g.WTS[2]-delta*s)]}
def boxes(inp,keys):
 p=np.array([float(F(x)) for x in inp['periods']]);ph=np.array([[float(F(x)) for x in row] for row in inp['phases']]);w=np.array([float(F(x)) for x in inp['weights']]);arr=np.asarray(keys,dtype=np.int64);lo=arr[:,:9].reshape(-1,3,3);hi=arr[:,9:].reshape(-1,3,3);cnt=np.prod(np.maximum(0,hi-lo+1),axis=2);mass=cnt@w;r=np.cbrt(mass)/2
 low=(lo+ph)*p;high=(hi+ph)*p;low[cnt==0]=np.inf;high[cnt==0]=-np.inf;lower=np.maximum(high.max(axis=1)-r[:,None],0);upper=np.minimum(low.min(axis=1)+r[:,None],p);good=np.all(upper>lower,axis=1)
 return np.stack([lower[good],upper[good]],axis=2),arr[good].tolist()
def main():
 d=json.loads((R/'round011d_direction.json').read_text());target=json.loads((R/'round011d_target.json').read_text());family=json.loads((R/'round011d_family.json').read_text());keys=family['keys'];assert set(map(tuple,g.SEED['keys']))<=set(map(tuple,keys));t=F(target['ratio_interval'][1]);expected=sha(R/'round011d_direction.json')
 out={'direction_sha256':expected,'family_sha256':sha(R/'round011d_family.json'),'input_sha256':sha(R/'round008_certificate_input.json'),'target_sha256':target['sha256'],'target_full_upper':str(t),'target_display':float(t),'baseline_zero_reused_no_union_query':True,'objective':'complete finite receiver-box union / exact constant W','full_W':str(g.W),'old639_preserved':True,'queries':[],'stop':None}
 start=time.monotonic()
 for step in d['steps'][1:]:
  s=F(step);inp=parameters(s,d);assert sum(map(F,inp['weights']))==g.W;pstart=time.monotonic();b,active=boxes(inp,keys);val=volume(b)/float(g.W);row={'s':str(s),'parameters':inp,'full_union_ratio_float':val,'float_excess_over_full022_upper':val-float(t),'active_boxes':len(b),'elapsed_seconds':time.monotonic()-pstart};out['queries'].append(row);out['elapsed_seconds']=time.monotonic()-start;dump('round011d_results.json',out);print(json.dumps(row),flush=True)
  if val>float(t)+1e-10:
   saved=sorted(set(map(tuple,active))|set(map(tuple,g.SEED['keys'])));inp['keys']=[list(k) for k in saved];dump('round011d_certificate_input.json',inp);out['candidate_input_sha256']=sha(R/'round011d_certificate_input.json');out['candidate_s']=str(s);out['stop']='first numerical candidate; exact certification before further queries';dump('round011d_results.json',out);return
 out['stop']='four fixed nonzero steps exhausted; no numerical candidate';out['elapsed_seconds']=time.monotonic()-start;dump('round011d_results.json',out)
if __name__=='__main__':main()
