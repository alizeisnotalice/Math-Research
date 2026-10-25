#!/usr/bin/env python3
"""Bounded exact rational two-cover qualification; no five-dimensional union."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import hashlib,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def root(n):
 d=2**80;lo=0;hi=8*d
 while lo+1<hi:
  m=(lo+hi)//2
  if m**5<=n*d**5:lo=m
  else:hi=m
 return Q(lo,d),Q(hi,d)
def main():
 start=time.monotonic();assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
 reg=json.loads((HERE/'registration_K8_n2.json').read_text());assert sha(Path(__file__))==reg['source_sha256']
 for row in reg['inputs']:assert sha(PROJECT/row['path'])==row['sha256']
 runtime=lambda:{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
 print(json.dumps({'stage':'started_before_math','runtime':runtime()}),flush=True)
 (HERE/'running_handle.json').write_text(json.dumps(runtime(),indent=2)+'\n')
 # Registration permits a short handoff before the bounded computation.
 if os.environ.get('LOWER018_HANDOFF')=='1':
  deadline=time.monotonic()+20
  while not (HERE/'start_authorized.flag').exists() and time.monotonic()<deadline:time.sleep(.05)
 start=time.monotonic();nodes=0
 reps=[(U,V,(4*U+3*V)//25,(-3*U+4*V)//25) for U,V in product(range(25),repeat=2) if (4*U+3*V)%25==0 and (-3*U+4*V)%25==0]
 assert len(reps)==25
 atoms=[]
 for U,V,m,n in reps:
  for h,k in product((-1,0,1),repeat=2):
   mm,nn=m+4*h+3*k,n-3*h+4*k;x,y=Q(U,5)+5*h,Q(V,5)+5*k
   assert x==Q(4*mm-3*nn,5) and y==Q(3*mm+4*nn,5)
   atoms.append((mm,nn,x,y))
 assert len(atoms)==len(set((m,n) for m,n,x,y in atoms))==225
 rho=Q(4,5);windows=[]
 for i,(m,n,x,y) in enumerate(atoms):
  b=max(Q(0),x-rho),min(Q(5),x+rho),max(Q(0),y-rho),min(Q(5),y+rho)
  if b[0]<b[1] and b[2]<b[3]:windows.append((i,b))
 xs=sorted({Q(0),Q(5)}|{t for i,b in windows for t in b[:2]});ys=sorted({Q(0),Q(5)}|{t for i,b in windows for t in b[2:]})
 xm=[sum(1<<i for i,b in windows if b[0]<=xl and xh<=b[1]) for xl,xh in zip(xs,xs[1:])]
 ym=[sum(1<<i for i,b in windows if b[2]<=yl and yh<=b[3]) for yl,yh in zip(ys,ys[1:])]
 histogram={};mincount=225;cellcert=[];under=[]
 for xi,(xl,xh) in enumerate(zip(xs,xs[1:])):
  for yi,(yl,yh) in enumerate(zip(ys,ys[1:])):
   nodes+=1
   if nodes>10000 or time.monotonic()-start>30:
    (HERE/'K8_n2_results.json').write_text(json.dumps({'status':'CHECKPOINT_LIMIT','nodes':nodes,'runtime':runtime()},indent=2)+'\n');return
   mask=xm[xi]&ym[yi];count=mask.bit_count();histogram[count]=histogram.get(count,Q(0))+(xh-xl)*(yh-yl);mincount=min(mincount,count)
   if count<2:under.append({'cell':[str(v) for v in (xl,xh,yl,yh)],'count':count})
   else:
    bit1=mask&-mask;bit2=(mask-bit1)&-(mask-bit1);ids=[bit1.bit_length()-1,bit2.bit_length()-1];assert len(set(ids))==2
    for i in ids:
     m,n,x,y=atoms[i];assert x-rho<=xl and xh<=x+rho and y-rho<=yl and yh<=y+rho
    cellcert.append([[str(t) for t in (xl,xh,yl,yh)],ids])
 assert sum(histogram.values())==25
 ul,uh=root(2);sl,sh=(1+ul)/2,(1+uh)/2;rl,rh=root(16);rl/=2;rh/=2
 assert ul**5<=2<=uh**5 and (2*rl)**5<=16<=(2*rh)**5
 margin=rl-rho*sh;strict_cover=mincount>=2 and margin>0
 # The K=8 shapes include (2,2,2). Its distinct unrotated index set is {0,1}^3.
 vertices=list(product((0,1),repeat=3));assert len(vertices)==len(set(vertices))==8
 activation=(2*rl,2*rh)
 examples=[]
 if strict_cover:
  for ids in (cellcert[0][1],cellcert[len(cellcert)//2][1],cellcert[-1][1]):
   full=[(j1,atoms[i][0],atoms[i][1],j4,j5) for j1,j4,j5 in vertices for i in ids]
   assert len(full)==len(set(full))==16;examples.append({'planar_atom_ids':ids,'distinct_full_atom_count':len(full)})
 result={'status':'CERTIFIED_WHOLE_PLANE_TWO_CAPTURE' if strict_cover else 'RATIONAL_TEST_ONLY_not_resolved','K':8,'n':2,'selected_mass':16,'common_exact_radius':'16^(1/5)/2','tail_spacing':'s=(1+2^(1/5))/2','rotation':[['4/5','-3/5'],['3/5','4/5']],'normalized_period':['5','5'],'rational_normalized_radius':str(rho),'minimum_closed_period_capture':mincount,'whole_plane_proof':'closed windows cover every closed cell of the finite arrangement with two distinct integer lattice indices; 5Z^2 periods cover the plane','finite_atom_completeness':'25 representatives of RZ^2 modulo 5Z^2; translation h,k in {-1,0,1} covers all centers intersecting [0,5]^2 for rho<1; omitted h or k has coordinate distance >rho','normalized_lattice_atoms':[{'id':i,'m':m,'n':n,'x':str(x),'y':str(y)} for i,(m,n,x,y) in enumerate(atoms)],'cell_witnesses':cellcert,'undercovered_cells':under,'area_count_histogram':{str(k):str(v) for k,v in sorted(histogram.items())},'strict_physical_radius_margin_lower':str(margin),'s_interval':[str(sl),str(sh)],'r_interval':[str(rl),str(rh)],'plane_capture_set':'entire [0,5s]^2 and its translates' if strict_cover else 'not certified','plane_area_exact':'25*s^2' if strict_cover else None,'plane_area_interval':[str(25*sl*sl),str(25*sh*sh)] if strict_cover else None,'a1_activation_event':'a1=2r=16^(1/5) for shape (2,2,2); active for 3/2<=a1<2r, empty at and above 2r','a1_activation_interval':[str(x) for x in activation],'parameter_window':['3/2','2'],'planar_capture_independent_of_a1':True,'full_distinct_mass_examples':examples,'first_axis_witness_for_shape222':'[max(0,a1-r), min(a1/2,r)]','tail_axis_witness_for_shape222':'[max(0,s-r), min(s/2,r)] on each remaining unrotated tail','ten_endpoint_guard':'each of the five center intervals lies inside [atom_max-r,atom_min+r]; planar four inequalities checked per cell at rho*s<r; first/tail six follow directly from displayed extrema formulas','resolved_cases':[[8,2]] if strict_cover else [],'remaining_unresolved_cases':[x for x in json.loads((PROJECT/'lower_round017/append_last_two/checkpoint_10of10.json').read_text())['unresolved_capture_cases'] if not strict_cover or x!=[8,2]],'source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration_K8_n2.json'),'nodes':nodes,'runtime':runtime(),'limits':{'seconds':30,'candidates':10000},'scope':'only planar capture qualification for one current finite case; no union, event-cell objective, parameter optimum or C5 result'}
 (HERE/'K8_n2_results.json').write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('normalized_lattice_atoms','cell_witnesses','undercovered_cells')},indent=2),flush=True)
if __name__=='__main__':main()
