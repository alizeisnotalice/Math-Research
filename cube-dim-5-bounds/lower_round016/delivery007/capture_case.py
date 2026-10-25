#!/usr/bin/env python3
"""One K=8,n=4 physical planar capture arrangement, exact enclosure."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
from decimal import Decimal,localcontext
import hashlib,json,time,sys,os,resource,threading
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
DEN=2**80;GRID=2**40;P,R,QC=4,3,5;K=int(sys.argv[1]);NEED=int(sys.argv[2]);SECONDS=120;NODES=200000
@lru_cache(None)
def roots(n):
 lo,hi=0,8*DEN
 while lo+1<hi:
  m=(lo+hi)//2
  if m**5<=n*DEN**5:lo=m
  else:hi=m
 if lo**5==n*DEN**5:lo-=1
 a,b=Q(lo,DEN),Q(lo+1,DEN);assert a**5<n<=b**5
 return a,b

def iv(c,ab):return (c*ab[0],c*ab[1]) if c>=0 else (c*ab[1],c*ab[0])
def down(x):
 z=x*GRID;return Q(z.numerator//z.denominator,GRID)
def up(x):return -down(-x)
def decimal(x):
 with localcontext() as c:c.prec=60;return str(Decimal(x.numerator)/Decimal(x.denominator))

def run():
 started=time.monotonic();nodes=0;stage='setup'
 def check(node=False):
  nonlocal nodes
  if node:nodes+=1
  if nodes>NODES:raise RuntimeError('node cap '+stage)
  if time.monotonic()-started>SECONDS:raise TimeoutError('wall cap '+stage)
 ul,uh=roots(2);ss=((1+ul)/2,(1+uh)/2);aa=(1+ul**3/2,1+uh**3/2)
 side_lo,side_hi=roots(K*NEED);rlo,rhi=side_lo/2,side_hi/2
 assert K in (1,2,3,4,6,8,9,12,16,18,24,27,32,36,48,64) and 1<=NEED<=16 and (2*rlo)**5<K*NEED
 reps=[]
 for U,V in product(range(QC*QC),repeat=2):
  if (P*U+R*V)%(QC*QC)==0 and (-R*U+P*V)%(QC*QC)==0:reps.append((U,V,(P*U+R*V)//25,(-R*U+P*V)//25))
 assert len(reps)==25
 atoms=[]
 for U,V,m,n in reps:
  for h,k in product((-1,0,1),repeat=2):
   atoms.append({'m':m+P*h+R*k,'n':n-R*h+P*k,'xcoef':Q(U,QC)+QC*h,'ycoef':Q(V,QC)+QC*k})
 assert len(atoms)==len({(p['m'],p['n']) for p in atoms})==225
 for p in atoms:
  assert p['xcoef']==Q(P*p['m']-R*p['n'],QC) and p['ycoef']==Q(R*p['m']+P*p['n'],QC)
 domains=(QC*ss[0],QC*ss[1]);windows=[[],[]];active_ids=[[],[]]
 def box_from_extrema(ex,mode):
  out=[]
  for mn,mx in ex:
   if mode==0:lo=up(max(Q(0),mx[1]-rlo));hi=down(min(domains[0],mn[0]+rlo))
   else:lo=down(max(Q(0),mx[0]-rhi));hi=up(min(domains[1],mn[1]+rhi))
   if hi<=lo:return None
   out.extend((lo,hi))
  return tuple(out)
 for i,p in enumerate(atoms):
  ex=[(iv(p[c],ss),iv(p[c],ss)) for c in ('xcoef','ycoef')]
  for mode in (0,1):
   b=box_from_extrema(ex,mode)
   if b:windows[mode].append(b);active_ids[mode].append(i)
 stage='unrotated_shape_screen'
 shapes=[lens for lens in product(range(1,5),repeat=3) if lens[0]*lens[1]*lens[2]==K]
 assert shapes
 def inds(ell):return tuple(range(-(ell-1)//2,(ell-1)//2+1)) if ell%2 else tuple(range(1-ell//2,1+ell//2))
 groups=[[],[]];excluded=[]
 for lens in shapes:
  bounds=[]
  for mode in (0,1):
   b=[]
   for ell,ab in zip(lens,(aa,ss,ss)):
    ii=inds(ell);mn=iv(Q(min(ii)),ab);mx=iv(Q(max(ii)),ab)
    if mode==0:lo=up(max(Q(0),mx[1]-rlo));hi=down(min(ab[0]/2,mn[0]+rlo))
    else:lo=down(max(Q(0),mx[0]-rhi));hi=up(min(ab[1]/2,mn[1]+rhi))
    b.extend((lo,hi))
   bounds.append(tuple(b))
  if any(bounds[1][2*q+1]<=bounds[1][2*q] for q in range(3)):excluded.append(lens);continue
  for mode in (0,1):
   if all(bounds[mode][2*q]<bounds[mode][2*q+1] for q in range(3)):groups[mode].append({'shape':list(lens),'box':[str(t) for t in bounds[mode]]})
 base={'K':K,'n':NEED,'selected_mass':K*NEED,'common_exact_radius':'(K*n)^(1/5)/2','unrotated_inner_outer_shapes_boxes':groups,'unrotated_shape_count':len(shapes),'excluded_by_unrotated_outer_empty':excluded,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'C_source_sha256':hashlib.sha256((HERE.parents[1]/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest(),'limits':{'seconds':SECONDS,'nodes':NODES}}
 if not groups[1]:
  return dict(base,status='CERTIFIED_EMPTY_all_unrotated_shapes',area_interval=['0','0'],nodes=nodes,wall_seconds=time.monotonic()-started)
 def analyze(bs,threshold,label,ids=None):
  nonlocal stage
  stage=label
  if not bs:return {'area':Q(0),'x_segments':0,'y_segments':0,'elementary_cells':0,'qualifying_cells':0,'witness_rectangles':[],'count_histogram':{}}
  xs=sorted({v for b in bs for v in b[:2]});ys=sorted({v for b in bs for v in b[2:]})
  xm=[sum(1<<i for i,b in enumerate(bs) if b[0]<=a and z<=b[1]) for a,z in zip(xs,xs[1:])]
  ym=[sum(1<<i for i,b in enumerate(bs) if b[2]<=a and z<=b[3]) for a,z in zip(ys,ys[1:])]
  area=Q(0);qualifying=0;strips=[];count_hist={}
  for xi,(xl,xh) in enumerate(zip(xs,xs[1:])):
   vertical=[];start=None;previous=None
   for yi,(yl,yh) in enumerate(zip(ys,ys[1:])):
    check(True);mask=xm[xi]&ym[yi];count=mask.bit_count();count_hist[count]=count_hist.get(count,0)+1
    if count>=threshold:
     qualifying+=1;area+=(xh-xl)*(yh-yl)
     chosen=[];bits=mask
     for _ in range(threshold):bit=bits&-bits;chosen.append(bit.bit_length()-1);bits-=bit
     witness=tuple(ids[i] for i in chosen) if ids is not None else tuple(chosen)
     if previous!=witness:
      if start is not None:vertical.append((start,yl,previous))
      start=yl;previous=witness
    else:
     if start is not None:vertical.append((start,yl,previous));start=None;previous=None
   if start is not None:vertical.append((start,ys[-1],previous))
   strips.append((xl,xh,vertical))
  # Independent integer y-event sweep per x slab; no covering bitmasks.
  second=Q(0)
  for xl,xh in zip(xs,xs[1:]):
   events={}
   for b in bs:
    if b[0]<=xl and xh<=b[1]:events[b[2]]=events.get(b[2],0)+1;events[b[3]]=events.get(b[3],0)-1
   count=0;last=None;length=Q(0)
   for y,diff in sorted(events.items()):
    check(True)
    if last is not None and count>=threshold:length+=y-last
    count+=diff;last=y
   assert count==0;second+=(xh-xl)*length
  assert area==second
  # Merge only rectangles carrying the SAME fixed atom witness tuple.
  merged=[];open_rects={}
  for xl,xh,vertical in strips:
   nxt={}
   for yl,yh,w in vertical:
    key=(yl,yh,w)
    if key in open_rects and open_rects[key][1]==xl:
     old=open_rects.pop(key);old[1]=xh;nxt[key]=old
    else:nxt[key]=[xl,xh,yl,yh,w]
   merged.extend(open_rects.values());open_rects=nxt
  merged.extend(open_rects.values())
  assert sum((b[1]-b[0])*(b[3]-b[2]) for b in merged)==area
  return {'area':area,'x_segments':len(xs)-1,'y_segments':len(ys)-1,'elementary_cells':(len(xs)-1)*(len(ys)-1),'qualifying_cells':qualifying,'witness_rectangles':merged,'count_histogram':count_hist}
 inner=analyze(windows[0],NEED,'inner_capture',active_ids[0]);outer=analyze(windows[1],NEED,'outer_capture',active_ids[1]);assert inner['area']<=outer['area']
 stage='direct_fixed_witness_guard';guards=0;witness_products=0;nonrect=0
 if outer['area']==0:
  assert inner['area']==0
  return dict(base,status='CERTIFIED_EMPTY_planar_capture',area_interval=['0','0'],maximum_outer_window_overlap=max(outer['count_histogram'],default=0),dual_area_engines_agree=True,nodes=nodes,wall_seconds=time.monotonic()-started)
 for shape_box in groups[0]:
  lens=tuple(shape_box['shape']);unrot=tuple(Q(t) for t in shape_box['box'])
  for b in inner['witness_rectangles']:
   xl,xh,yl,yh,ids=b;assert len(ids)==NEED==len(set(ids));chosen=[atoms[i] for i in ids]
   pts=[(j1,p['m'],p['n'],j4,j5) for j1,j4,j5 in product(*(inds(ell) for ell in lens)) for p in chosen];assert len(pts)==len(set(pts))==K*NEED
   full=unrot+(xl,xh,yl,yh);ex=[]
   for q,ab in zip((0,3,4),(aa,ss,ss)):ex.append((iv(Q(min(t[q] for t in pts)),ab),iv(Q(max(t[q] for t in pts)),ab)))
   for c in ('xcoef','ycoef'):ex.append((iv(min(p[c] for p in chosen),ss),iv(max(p[c] for p in chosen),ss)))
   for q,(mn,mx) in enumerate(ex):assert mx[1]-rlo<=full[2*q] and full[2*q+1]<=mn[0]+rlo;guards+=2
   if len({p['m'] for p in chosen})*len({p['n'] for p in chosen})!=NEED:nonrect+=1
   witness_products+=1;check(True)
 atom_data=[{'id':i,'m':a['m'],'n':a['n'],'x_coefficient_of_s':str(a['xcoef']),'y_coefficient_of_s':str(a['ycoef'])} for i,a in enumerate(atoms)]
 atom_text=json.dumps(atom_data,indent=2)+'\n';atom_path=HERE/'rotation_atoms.json'
 if atom_path.exists():assert atom_path.read_text()==atom_text
 else:atom_path.write_text(atom_text)
 packet={'atom_manifest_path':'rotation_atoms.json','atom_manifest_sha256':hashlib.sha256(atom_path.read_bytes()).hexdigest(),'inner_fixed_atom_rectangles':[{'rectangle':[str(t) for t in b[:4]],'atom_ids':list(b[4])} for b in inner['witness_rectangles']],'outer_count_rectangles':[{'rectangle':[str(t) for t in b[:4]],'covering_window_ids':list(b[4])} for b in outer['witness_rectangles']]}
 packet_path=HERE/f'K{K}_n{NEED}_partition.json';packet_path.write_text(json.dumps(packet,indent=2)+'\n')
 strip=lambda d:{k:v for k,v in d.items() if k not in ('witness_rectangles','area')}
 return dict(base,status='PASS_exact_all_shapes_physical_capture',strict_inner_radius=str(rlo),area_interval=[str(inner['area']),str(outer['area'])],area_decimal_display=[decimal(inner['area']),decimal(outer['area'])],inner_arrangement=strip(inner),outer_arrangement=strip(outer),compressed_rectangle_counts=[len(inner['witness_rectangles']),len(outer['witness_rectangles'])],direct_geometry_endpoint_inequalities=guards,checked_witness_products=witness_products,distinct_mass_per_witness=K*NEED,nonrectangular_witness_products=nonrect,dual_area_engines_agree=True,partition_sha256=hashlib.sha256(packet_path.read_bytes()).hexdigest(),atom_manifest_sha256=hashlib.sha256(atom_path.read_bytes()).hexdigest(),nodes=nodes,wall_seconds=time.monotonic()-started)
def main():
 assert os.getpriority(os.PRIO_PROCESS,0)>=10
 assert all(os.environ.get(k)=='1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'))
 started=time.monotonic()
 try:out=run()
 except (RuntimeError,TimeoutError) as e:out={'status':'COMPUTE_LIMIT','K':K,'n':NEED,'reason':str(e),'wall_seconds':time.monotonic()-started}
 out['portable_runtime']={'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if sys.platform=='darwin' else resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
 (HERE/f'K{K}_n{NEED}_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'K':K,'n':NEED,'status':out['status'],'nodes':out.get('nodes'),'wall_seconds':out['wall_seconds']}))
if __name__=='__main__':main()
