"""High precision exploration of uniform arc circuits; not interval proofs."""
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path
import mpmath as mp
import platform,json,time
BITS=list(product([-1,1],repeat=4))
def run(job):
 ts,zs=job;mp.mp.dps=70;t=mp.mpf(ts);zeta=mp.mpf(zs);a=1-t;c=2*a*a-1;theta=mp.pi/4+zeta*t
 co=mp.cos(theta);si=mp.sin(theta);rb=(a-mp.sqrt(c*mp.tan(theta)))/2;sb=(a-mp.sqrt(c/mp.tan(theta)))/2
 d=mp.mpf(1)/64;h=k=mp.mpf(1)/4096;r=rb-d*t**3*co;s=sb-d*t**3*si;q=a+t*t/8;M=2*q-a
 mi=[mp.inf]*4
 for nsum in [-2,0,2]:
  for ssum in [-2,0,2]:
   dx=2*r+h*nsum*t**3*co-k*ssum*t**mp.mpf('1.5')*si
   dy=2*s+h*nsum*t**3*si+k*ssum*t**mp.mpf('1.5')*co
   gg=lambda zz:(zz-dx)*(zz-dy)-2*a*zz+1
   mi[0]=min(mi[0],gg(a)/t**3);mi[1]=min(mi[1],-gg(q)/t**3);mi[2]=min(mi[2],-gg(1/a)/t**3)
 for ns in BITS:
  for ss in BITS:
   cs=[]
   for (sx,sy),n,v in zip([(1,1),(-1,1),(-1,-1),(1,-1)],ns,ss):
    # Reflections of the NE normal/tangent frame.
    cs.append((sx*(r+h*n*t**3*co-k*v*t**mp.mpf('1.5')*si),sy*(s+h*n*t**3*si+k*v*t**mp.mpf('1.5')*co)))
   area=mp.mpf(0)
   for mask in range(1,16):
    sel=[cs[i] for i in range(4) if mask>>i&1]
    X=max(x for x,y in sel)-min(x for x,y in sel);Y=max(y for x,y in sel)-min(y for x,y in sel)
    area+=(-1)**(len(sel)+1)*(a-X)*(a-Y)
   mi[3]=min(mi[3],(area-4*M*(M-a)-1)/t**2)
 assert min(mi)>0
 U=mp.quad(lambda v:a*((1/a-a)-v)/(a-v),[0,1/a-a])
 return {'t':ts,'zeta':zs,'min_normalized_margins':[mp.nstr(x,55) for x in mi],
         'model_gap_over_t4':mp.nstr((2-t*t)*U/t**2,55)}
if __name__=='__main__':
 start=time.monotonic();ts=['.01','.001','.0001','.000001','.00000001'];zs=['-.125','-.0625','-.03125','0','.03125','.0625','.125'];jobs=list(product(ts,zs))
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 out={'python':platform.python_version(),'mpmath':mp.__version__,'precision':70,'workers':8,'parameter_pairs':len(rows),'opposite_pair_corner_cases':len(rows)*9,'full_union_corner_cases':len(rows)*256,'rows':rows,'seconds':time.monotonic()-start,'scope':'finite high precision regressions only; the uniform small-t theorem uses analytic compact Taylor bounds and continuous geometric proofs'}
 Path(__file__).resolve().parents[1].joinpath('results/arc_aggregation_python.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='rows'})
