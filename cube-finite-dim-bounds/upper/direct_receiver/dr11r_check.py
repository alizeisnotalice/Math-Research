"""Exactly six two-receiver-box LP certificates; finite geometry, not C_n bounds."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,combinations
from concurrent.futures import ProcessPoolExecutor
from math import prod
import json,hashlib,os,time

CONFIG_SHA='e5c1c21eb66f61765016d8351a204d0a9f5cc88e88e845c1972e23f7476ae6d3'
def overlap(y,l,u,r):return max(Q(0),min(u,y+r)-max(l,y-r))
def coeff(y,box,r):return prod(overlap(v,Q(l),Q(u),Q(r)) for v,(l,u) in zip(y,box))/Q((2*r)**len(box))
def volume(box):return prod(Q(u-l) for l,u in box)
def feasible(h,rows):return min(h)>=0 and all(a*h[0]+b*h[1]<=1 for a,b in rows)

def solve(config):
 start=time.time();n=config['n'];boxes=config['boxes'];rs=config['radii'];vols=list(map(volume,boxes))
 assert prod(max(Q(0),Q(min(b[1],c[1])-max(b[0],c[0]))) for b,c in zip(*boxes))==0
 axes=[]
 for j in range(n):
  axes.append(sorted({Q(x+sg*r) for box,r in zip(boxes,rs) for x in box[j] for sg in (-1,1)}))
 coeffs={}
 for y in product(*axes):
  aa=tuple(coeff(y,b,r) for b,r in zip(boxes,rs));coeffs.setdefault(aa,y)
 front=[];bestb=Q(-1)
 for aa in sorted(coeffs,reverse=True):
  if aa[1]>bestb:front.append(aa);bestb=aa[1]
 front=[aa for aa in front if any(aa)]
 candidates=[(Q(0),Q(0)),(min(1/a for a,b in front if a),Q(0)),(Q(0),min(1/b for a,b in front if b))]
 for (a,b),(c,d) in combinations(front,2):
  det=a*d-b*c
  if det:candidates.append(((d-b)/det,(a-c)/det))
 candidates=[h for h in candidates if feasible(h,front)]
 h=max(candidates,key=lambda h:(sum(x*y for x,y in zip(h,vols)),h[1],h[0]))
 mass=sum(x*y for x,y in zip(h,vols));total=sum(vols)
 fullmax=max(sum(x*y for x,y in zip(h,aa)) for aa in coeffs);assert fullmax==1
 active=[aa for aa in front if sum(x*y for x,y in zip(h,aa))==1]
 dual=None
 for (a,b),(c,d) in combinations(active,2):
  det=a*d-b*c
  if not det:continue
  lam=((vols[0]*d-vols[1]*c)/det,(vols[1]*a-vols[0]*b)/det)
  if min(lam)<0:continue
  assert lam[0]*a+lam[1]*c==vols[0] and lam[0]*b+lam[1]*d==vols[1]
  assert sum(lam)==mass
  dual=[{'input_point':list(map(str,coeffs[aa])),'load_vector':list(map(str,aa)),'price_exact':str(p)} for aa,p in zip(((a,b),(c,d)),lam)];break
 assert dual is not None
 um=max(a+b for a,b in coeffs);uh=1/um;umass=uh*total
 centers=[tuple(Q(l+u,2) for l,u in b) for b in boxes]
 centerloads=[sum(coeff(y,b,r) for b,r in zip(boxes,rs)) for y in centers]
 assert max(centerloads)==1 and um>1
 witness=max(coeffs,key=lambda aa:sum(aa));yp=coeffs[witness]
 if config['kind']=='volume_critical':assert total/mass==Q(5,4)
 elif config['kind']=='face_adjacent':assert total/mass==Q(8,7)
 else:assert total/mass==Q(2*4**n,2*4**n-1)
 constraint_hash=hashlib.sha256(json.dumps([list(map(str,aa)) for aa in sorted(coeffs)]).encode()).hexdigest()
 return {'id':config['id'],'n':n,'kind':config['kind'],'box_volumes':list(map(str,vols)),
         'breakpoints':[[str(x) for x in z] for z in axes],'grid_vertex_count':prod(map(len,axes)),
         'unique_load_vectors':len(coeffs),'load_vectors_sha256':constraint_hash,
         'Pareto_constraints':[{'vector':list(map(str,aa)),'witness':list(map(str,coeffs[aa]))} for aa in front],
         'optimal_h_exact':list(map(str,h)),'optimal_receiver_mass_exact':str(mass),'receiver_measure_exact':str(total),
         'U_receiver_constant_class_exact':str(total/mass),'U_receiver_constant_class_display':float(total/mass),
         'global_Th_max_exact':str(fullmax),'finite_LP_dual_atoms':dual,
         'uniform_1_weights_global_load_exact':str(um),'uniform_1_load_witness':list(map(str,yp)),
         'receiver_center_uniform_loads':list(map(str,centerloads)),
         'best_uniform_h_exact':str(uh),'best_uniform_mass_exact':str(umass),
         'nonuniform_mass_gain_exact':str(mass/umass),'elapsed_seconds':time.time()-start,
         'scope':'optimum only among h constant on each of these two boxes; no general E or all-scale C_n claim'}

def main():
 p=Path(__file__).parent;f=p/'dr11r_configs.json';assert hashlib.sha256(f.read_bytes()).hexdigest()==CONFIG_SHA
 configs=json.loads(f.read_text())['configs'];assert len(configs)==6
 start=time.time();workers=min(6,os.cpu_count() or 1)
 with ProcessPoolExecutor(max_workers=workers) as ex:rows=list(ex.map(solve,configs))
 out={'batch':'11R_direct_receiver','status':'all six exact continuous-capacity and finite-LP certificates PASS',
      'config_sha256':CONFIG_SHA,'worker_processes':workers,'total_grid_vertices':sum(z['grid_vertex_count'] for z in rows),
      'rows':rows,'elapsed_seconds':time.time()-start,'new_Cn_upper_bound':False,
      'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (p/'dr11r_results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
 for z in rows:print(z['id'],'h',z['optimal_h_exact'],'U',z['U_receiver_constant_class_exact'],'grid',z['grid_vertex_count'],flush=True)
 print('PASS',workers,'processes',out['elapsed_seconds'],'seconds',flush=True)
if __name__=='__main__':main()
