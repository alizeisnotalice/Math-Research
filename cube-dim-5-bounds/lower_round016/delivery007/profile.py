#!/usr/bin/env python3
"""One serial low-priority exact radius-profile checkpoint per invocation."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from collections import defaultdict
import hashlib,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;BASE=HERE.parent
DEN=2**80;GRID=2**40
THREAD_ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def root2():
 lo,hi=0,2*DEN
 while lo+1<hi:
  m=(lo+hi)//2
  if m**5<=2*DEN**5:lo=m
  else:hi=m
 return Q(lo,DEN),Q(hi,DEN)
def down(x):z=x*GRID;return z.numerator//z.denominator
def up(x):return -down(-x)
def iv(c,ab):return (c*ab[0],c*ab[1]) if c>=0 else (c*ab[1],c*ab[0])
def indices(ell):
 return tuple(range(-(ell-1)//2,(ell-1)//2+1)) if ell%2 else tuple(range(1-ell//2,1+ell//2))
def runtime(start):
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),
         'python_threads':threading.active_count(),
         'thread_environment':{k:os.environ.get(k) for k in THREAD_ENV},
         'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,
         'wall_seconds':time.monotonic()-start}

def geometry():
 ul,uh=root2();aa=(1+ul**3/2,1+uh**3/2);ss=((1+ul)/2,(1+uh)/2)
 ap=HERE/'rotation_atoms.json';atoms=json.loads(ap.read_text())
 assert len(atoms)==len({(a['m'],a['n']) for a in atoms})==225
 return aa,ss,atoms,ap

def overlap(side,ss,atoms):
 started=time.monotonic();nodes=0
 def check(node=False):
  nonlocal nodes
  if node:nodes+=1
  if nodes>200000:raise RuntimeError('overlap node checkpoint')
  if time.monotonic()-started>120:raise TimeoutError('overlap wall checkpoint')
 maxima=[];stats=[]
 for mode in (0,1):
  bs=[]
  for atom in atoms:
   box=[]
   for key in ('x_coefficient_of_s','y_coefficient_of_s'):
    low,high=iv(Q(atom[key]),ss)
    if mode==0:lo=up(max(Q(0),high-side[0]/2));hi=down(min(5*ss[0],low+side[0]/2))
    else:lo=down(max(Q(0),low-side[1]/2));hi=up(min(5*ss[1],high+side[1]/2))
    if hi<=lo:break
    box.extend((lo,hi))
   else:bs.append(tuple(box))
  xs=sorted({v for b in bs for v in b[:2]});ys=sorted({v for b in bs for v in b[2:]})
  xm=[sum(1<<i for i,b in enumerate(bs) if b[0]<=x and z<=b[1]) for x,z in zip(xs,xs[1:])]
  ym=[sum(1<<i for i,b in enumerate(bs) if b[2]<=y and z<=b[3]) for y,z in zip(ys,ys[1:])]
  maximum=0
  for xmask in xm:
   for ymask in ym:check(True);maximum=max(maximum,(xmask&ymask).bit_count())
  second=0
  for xl,xh in zip(xs,xs[1:]):
   events=defaultdict(int)
   for b in bs:
    if b[0]<=xl and xh<=b[1]:events[b[2]]+=1;events[b[3]]-=1
   count=0;last=None
   for y,diff in sorted(events.items()):
    check()
    if last is not None and y>last:second=max(second,count)
    count+=diff;last=y
   assert count==0
  assert maximum==second;maxima.append(maximum)
  stats.append({'mode':mode,'windows':len(bs),'segments':[len(xs)-1,len(ys)-1]})
 assert maxima[0]==maxima[1],('increase precision before certification',maxima)
 return {'exact_maximum':maxima[0],'inner_outer_maxima':maxima,'dual_cell_and_event_algorithms_agree':True,'nodes':nodes,'wall_seconds':time.monotonic()-started,'stats':stats}

def profile(i):
 start=time.monotonic();aa,ss,atoms,ap=geometry();sp=HERE/'profile_setup_v4.json';spec=json.loads(sp.read_text())
 assert sha(ap)==spec['atom_manifest_sha256'];left,right=spec['critical_points'][i:i+2]
 lb=tuple(Q(t) for t in left['interval']);rb=tuple(Q(t) for t in right['interval']);side=((lb[0]+rb[0])/2,(lb[1]+rb[1])/2)
 capture=overlap(side,ss,atoms)
 def maximum1D(ab):
  floor=0
  for k in range(1,10):
   if k*ab[1]<side[0]:floor=k
   elif k*ab[0]>side[1]:break
   else:raise ArithmeticError('unseparated1D event')
  return floor+1
 m1,mt=maximum1D(aa),maximum1D(ss);maximum=m1*mt*mt*capture['exact_maximum'];excluded=Q(maximum)<=lb[0]**5
 cover=[];extra=Q(0)
 if not excluded:
  for lens in product(range(1,m1+1),range(1,mt+1),range(1,mt+1)):
   K=lens[0]*lens[1]*lens[2]
   if K*capture['exact_maximum']<=lb[0]**5:continue
   box=[];v=Q(1)
   for ell,ab in zip(lens,(aa,ss,ss)):
    ii=indices(ell);mn=iv(Q(min(ii)),ab);mx=iv(Q(max(ii)),ab)
    lo=max(Q(0),mx[0]-rb[1]/2);hi=min(ab[1]/2,mn[1]+rb[1]/2)
    if hi<=lo:v=Q(0);break
    box.extend((lo,hi));v*=hi-lo
   if v:
    ratio=Q(8,25)*v*(5*ss[1])**2;extra+=ratio
    cover.append({'shape':lens,'box_outer':[str(t) for t in box],'extra_ratio_upper':str(ratio)})
 out={'status':'PASS_exact_radius_interval_profile','interval_index':i,'left':left,'right':right,'interior_side_interval':[str(t) for t in side],
      'capture':capture,'max_unrotated_counts':[m1,mt,mt],'max_total_mass':maximum,'excluded_for_all_sides_in_open_interval':excluded,
      'large_side_cover':cover,'large_side_extra_ratio_upper':str(extra),'source_sha256':sha(Path(__file__)),'setup_sha256':sha(sp),'runtime':runtime(start),
      'scope':'open side interval; strict critical-side witnesses are covered via right enlargement'}
 cache=HERE/'replay_checkpoints';cache.mkdir(exist_ok=True)
 (cache/f'profile_interval_{i:03}_v4.json').write_text(json.dumps(out,indent=2)+'\n');return out

if __name__=='__main__':
 for key in THREAD_ENV:assert os.environ.get(key)=='1',(key,'must be1')
 assert os.getpriority(os.PRIO_PROCESS,0)>=10,'run with nice+10 or lower priority'
 print(json.dumps({'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'checkpoint':sys.argv[1:]}),flush=True)
 out=profile(int(sys.argv[1]))
 print(json.dumps({k:v for k,v in out.items() if k not in ('critical_points','large_side_cover')},indent=2))
