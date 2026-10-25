#!/usr/bin/env python3
"""Portable exact certificate replay; one serial low-priority task per call."""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
import hashlib,json,os,resource,sys,threading,time
HERE=Path(__file__).resolve().parent;CACHE=HERE/'replay_checkpoints'
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name):return json.loads((HERE/name).read_text())
def runtime(start):
 peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':peak if sys.platform=='darwin' else peak*1024,'wall_seconds':time.monotonic()-start}
def assets():
 if (HERE/'manifest.json').exists():
  manifest=load('manifest.json')
  for row in manifest['files']:assert sha(HERE/row['path'])==row['sha256'],row['path']
  ref=manifest['existing_reference_dependency'];assert sha(HERE/ref['path'])==ref['sha256']

def arithmetic(histograms,bundle,result):
 G=bundle['grid'];families=bundle['families'];answers=[];old2=[];old3=[]
 for mode in (0,1):
  ph=next(h['integer_units'] for h in histograms if h['kind']=='plane' and h['mode']==mode)
  bh=next(h['integer_units'] for h in histograms if h['kind']=='unrot' and h['mode']==mode);joint=defaultdict(int)
  for p,a in ph.items():
   for b,v in bh.items():
    mask=int(p)&int(b)
    if mask:joint[mask]+=a*v
  factor=Q(8,25*G**5)
  def query(Ks):
   active=sum(1<<i for i,f in enumerate(families) if f['K'] in Ks)
   return factor*sum(v for m,v in joint.items() if m&active)
  answers.append(factor*sum(joint.values()));old2.append(query({1,8}));old3.append(query({1,2,8}))
 assert answers==list(map(Q,result['ratio_interval']))
 assert old2==list(map(Q,result['prior_signature_ratios']['K1_K8']))
 assert old3==list(map(Q,result['prior_signature_ratios']['K1_K2_K8']))
 return answers

def profile_certificate():
 setup=load('profile_setup_v4.json');cert=load('radius_certificate.json')
 assert Q(setup['F8_rigorous_upper'])<1 and setup['side4_overlap']['exact_maximum']==15
 assert cert['setup_sha256']==sha(HERE/'profile_setup_v4.json') and len(cert['interval_results'])==setup['open_intervals']==22
 for i,d in enumerate(cert['interval_results']):
  assert d['interval_index']==i and d['left']==setup['critical_points'][i] and d['right']==setup['critical_points'][i+1]
  assert d['excluded_for_all_sides_in_open_interval'] and d['capture']['dual_cell_and_event_algorithms_agree']
  assert d['capture']['inner_outer_maxima'][0]==d['capture']['inner_outer_maxima'][1]
  assert Q(d['max_total_mass'])<=Q(d['left']['interval'][0])**5
 return True

def geometry(kind,mode,bundle):
 bs=[];which=[]
 for fno,f in enumerate(bundle['families']):
  if kind=='plane':group=f['inner_planar'] if mode==0 else f['outer_planar'];boxes=[b['box'] for b in group]
  else:boxes=f['inner_boxes'] if mode==0 else f['outer_boxes']
  for b in boxes:bs.append(tuple(b));which.append(fno)
 dim=2 if kind=='plane' else 3;axes=[sorted({b[2*q+k] for b in bs for k in (0,1)}) for q in range(dim)]
 return bs,which,axes

def block(kind,mode,begin,bundle):
 start=time.monotonic();nodes=0;bs,which,axes=geometry(kind,mode,bundle);width=36 if kind=='plane' else 16;end=min(begin+width,len(axes[0])-1);assert 0<=begin<end
 nf=len(bundle['families']);fm=[sum(1<<i for i,f in enumerate(which) if f==j) for j in range(nf)]
 masks=[[sum(1<<i for i,b in enumerate(bs) if b[2*q]<=lo and hi<=b[2*q+1]) for lo,hi in zip(a,a[1:])] for q,a in enumerate(axes)]
 hist=defaultdict(int);second=defaultdict(int)
 for xi in range(begin,end):
  xl,xh=axes[0][xi:xi+2];dx=xh-xl
  for yi,(yl,yh) in enumerate(zip(axes[1],axes[1][1:])):
   xy=masks[0][xi]&masks[1][yi];dy=yh-yl
   zs=[None] if kind=='plane' else list(zip(axes[2],axes[2][1:]))
   for zi,slab in enumerate(zs):
    nodes+=1;assert nodes<=200000 and time.monotonic()-start<120
    rmask=xy if slab is None else xy&masks[2][zi];mask=sum(1<<f for f,m in enumerate(fm) if rmask&m)
    if mask:hist[mask]+=dx*dy*(1 if slab is None else slab[1]-slab[0])
  yslabs=[None] if kind=='plane' else list(zip(axes[1],axes[1][1:]))
  for slab in yslabs:
   events=defaultdict(lambda:[0]*nf)
   for b,f in zip(bs,which):
    if b[0]<=xl and xh<=b[1] and (slab is None or (b[2]<=slab[0] and slab[1]<=b[3])):
     q=2 if kind=='plane' else 4;events[b[q]][f]+=1;events[b[q+1]][f]-=1
   counts=[0]*nf;last=None;factor=dx if slab is None else dx*(slab[1]-slab[0])
   for z,changes in sorted(events.items()):
    mask=sum(1<<f for f,c in enumerate(counts) if c)
    if last is not None and mask:second[mask]+=factor*(z-last)
    counts=[a+b for a,b in zip(counts,changes)];assert all(c>=0 for c in counts);last=z
   assert not any(counts)
 assert dict(hist)==dict(second)
 out={'status':'PASS_portable_exact_histogram_block','kind':kind,'mode':mode,'x_indices':[begin,end],'x_grid_bounds':[axes[0][begin],axes[0][end]],'total_segments':[len(a)-1 for a in axes],'integer_units':{str(k):v for k,v in sorted(hist.items())},'dual_algorithms_agree':True,'nodes':nodes,'runtime':runtime(start)}
 (CACHE/f'{kind}_m{mode}_x{begin:04}.json').write_text(json.dumps(out,indent=2)+'\n');return out

def from_blocks(bundle):
 out=[]
 for mode in (0,1):
  for kind in ('plane','unrot'):
   _,_,axes=geometry(kind,mode,bundle);width=36 if kind=='plane' else 16;hist=defaultdict(int);last=0
   for begin in range(0,len(axes[0])-1,width):
    pp=CACHE/f'{kind}_m{mode}_x{begin:04}.json'
    if not pp.exists():raise FileNotFoundError('Replay one block first: '+f'block {kind} {mode} {begin}')
    d=json.loads(pp.read_text());assert d['x_indices'][0]==last and d['dual_algorithms_agree'];last=d['x_indices'][1]
    for m,v in d['integer_units'].items():hist[m]+=v
   assert last==len(axes[0])-1
   out.append({'kind':kind,'mode':mode,'integer_units':dict(hist)})
 return out

def main():
 start=time.monotonic();assert all(os.environ.get(k)=='1' for k in ENV);assert os.getpriority(os.PRIO_PROCESS,0)>=10
 assets();CACHE.mkdir(exist_ok=True);bundle=load('geometry_bundle.json');result=load('result.json');cmd=sys.argv[1] if len(sys.argv)>1 else 'verify'
 if cmd=='block':
  kind=sys.argv[2];mode,begin=map(int,sys.argv[3:]);assert kind in ('plane','unrot') and mode in (0,1);d=block(kind,mode,begin,bundle)
  print(json.dumps({k:v for k,v in d.items() if k!='integer_units'},indent=2));return
 hist=from_blocks(bundle) if cmd=='rollup-blocks' else load('aggregate_histograms.json')['histograms']
 answer=arithmetic(hist,bundle,result);assert profile_certificate()
 reference=json.loads((HERE/result['reference_path']).read_text());assert answer[1]<Q(reference['ratio_interval'][0])
 d={'status':'PASS_portable_complete_fixed_input_certificate','ratio_interval':[str(x) for x in answer],'saved_profile_and_exact_prior_signatures_verified':True,'old_capture_sources_rerun':False,'full_geometry_union_rebuilt':cmd=='rollup-blocks','runtime':runtime(start)}
 (CACHE/'last_certificate_replay.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d,indent=2))
if __name__=='__main__':main()
