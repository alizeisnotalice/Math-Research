#!/usr/bin/env python3
"""One novel a1 finite five-dimensional union, via a joint planar stripe mask."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from collections import defaultdict
from decimal import Decimal,localcontext
import hashlib,json,os,sys,resource,threading,time
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parent;G=5*2**40;A=Q(151,100)
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def save(p,d):p.write_text(json.dumps(d,separators=(',',':'))+'\n')
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
@lru_cache(None)
def roots(n):
 D=2**80;lo=0;hi=8*D
 while lo+1<hi:
  mid=(lo+hi)//2
  if mid**5<=n*D**5:lo=mid
  else:hi=mid
 if lo**5==n*D**5:lo-=1
 assert lo**5<n*D**5<=hi**5
 return Q(lo,D),Q(hi,D)
def runtime(start):
 r=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return {'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'threads':threading.active_count(),'peak_rss_bytes':r*(1 if sys.platform=='darwin' else 1024),'wall_seconds':time.monotonic()-start,'thread_environment':{k:os.environ.get(k) for k in ENV}}
def dec(q):
 with localcontext() as c:c.prec=60;return str(Decimal(q.numerator)/Decimal(q.denominator))
def inds(ell):return tuple(range(-(ell-1)//2,(ell-1)//2+1)) if ell%2 else tuple(range(1-ell//2,1+ell//2))
def volume(b):
 v=1
 for i in range(0,len(b),2):v*=max(0,b[i+1]-b[i])
 return v
class Limit(Exception):pass
def main():
 start=time.monotonic();nodes=0;stage='inputs'
 def check(n=1):
  nonlocal nodes
  nodes+=n
  if nodes>200000 or time.monotonic()-start>120:raise Limit(f'{stage}:120seconds_or200000nodes')
 reg=load(HERE/'registration.json');assert sha(Path(__file__))==reg['source_sha256']
 for row in reg['inputs']:assert sha(PROJECT/row['path'])==row['sha256']
 print(json.dumps({'stage':'before_finite_union','runtime':runtime(start),'a1':str(A),'grid':G}),flush=True)
 if os.environ.get('LOWER021_HANDOFF')=='1':
  deadline=time.monotonic()+20
  while not (HERE/'start.flag').exists() and time.monotonic()<deadline:time.sleep(.05)
 start=time.monotonic()
 try:
  raw=load(PROJECT/'lower_round017/raw_endpoint_catalog.json');events=load(PROJECT/'lower_round017/append_last_two/complete_events_10of10.json');assert len(raw['raw_cases'])==1024 and len(events['ordered_events'])==1947
  interval=None
  for left,right in zip(events['ordered_events'],events['ordered_events'][1:]):
   if Q(left['interval'][1])<A<Q(right['interval'][0]):interval={'left_event':left,'right_event':right};break
  assert interval is not None
  ul,uh=roots(2);sl,sh=(1+ul)/2,(1+uh)/2
  old=load(PROJECT/'lower_round016/delivery007/case_evidence.json');oldmap={(r['K'],r['n']):r for r in old['cases']};new=load(PROJECT/'lower_round020/results.json');newmap={(r['K'],r['n']):r for r in new['all_case_exact_formulas']}
  n18={tuple(x) for x in load(PROJECT/'lower_round017/raw_endpoint_catalog.json')['not_previously_evaluated_planar_cases_needed']};assert len(n18)==18
  zero20={key for key,r in newmap.items() if r['normalized_exact_area_polynomial']==['0','0','0']}
  oldzero={key for key,r in oldmap.items() if r['status']!='CERTIFIED_EMPTY_all_unrotated_shapes' and not r['status'].startswith('SKIP_') and list(map(Q,r['area_interval']))==[0,0]}
  fam={};inventory=[];stage='all1024_raw_inventory'
  for r in raw['raw_cases']:
   check();K,n,N=r['K'],r['n'],r['N'];key=(K,n);lens=r['shape'];reason=None
   if r['uniform_geometric_exclusion']:reason='017_uniform_span_exclusion'
   elif key in zero20:reason='020_exact_zero_planar_area'
   elif key in oldzero:reason='old_certified_zero_planar_area'
   if reason:
    inventory.append({'shape':lens,'K':K,'n':n,'status':reason});continue
   assert r['planar_status']=='evaluated_nonzero' or key in n18,(key,r['planar_status'])
   nl,nh=roots(N);rl,rh=nl/2,nh/2
   alpha=floor(rl/A*G);beta=floor(rl/sh*G);alpha_out=ceil(rh/A*G);beta_out=ceil(rh/sl*G)
   inner=[];outer=[];indexsets=[]
   for ell,t,ot in zip(lens,(alpha,beta,beta),(alpha_out,beta_out,beta_out)):
    ii=inds(ell);assert len(ii)==ell;indexsets.append(ii)
    inner.extend((max(0,max(ii)*G-t),min(G//2,min(ii)*G+t)))
    outer.extend((max(0,max(ii)*G-ot),min(G//2,min(ii)*G+ot)))
   if not volume(outer):inventory.append({'shape':lens,'K':K,'n':n,'status':'candidate_outer_unrotated_empty'});continue
   # All six endpoint inequalities against the same strict common r_lower.
   for j,(ii,t,scale_upper) in enumerate(zip(indexsets,(alpha,beta,beta),(A,sh,sh))):
    if inner[2*j]<inner[2*j+1]:
     assert Q(max(ii))-Q(inner[2*j],G)<=Q(t,G) and Q(inner[2*j+1],G)-Q(min(ii))<=Q(t,G)
     assert Q(t,G)*scale_upper<=rl
   assert (2*rl)**5<N
   if key not in fam:fam[key]={'K':K,'n':n,'N':N,'alpha_grid':alpha,'rho_grid':beta,'rho_outer_interval_upper':str(rh/sl),'strict_common_radius':str(rl),'nominal_radius_interval':[str(rl),str(rh)],'inner_boxes':[],'outer_boxes':[],'shapes':[]}
   fam[key]['inner_boxes'].append(inner);fam[key]['outer_boxes'].append(outer);fam[key]['shapes'].append(lens)
   inventory.append({'shape':lens,'K':K,'n':n,'status':'included_in_finite_union','positive_inner_box':bool(volume(inner))})
  families=[fam[k] for k in sorted(fam)];fi={key:i for i,key in enumerate(sorted(fam))}
  for i,f in enumerate(families):f['family_id']=i
  all18=[]
  for key in sorted(n18):
   statuses=sorted({r['status'] for r in inventory if (r['K'],r['n'])==key});all18.append({'K':key[0],'n':key[1],'candidate_raw_statuses':statuses,'included_family_id':fi.get(key),'exact_plane_status_source':'017 analyticK8n1 /018 K8n2 /020 stratum formulas','not_ignored':True})
  asset={'a1':str(A),'s_interval':[str(sl),str(sh)],'grid':G,'raw_case_count':len(inventory),'raw_inventory':inventory,'new18_accounting':all18,'families':families,'event_component_containing_candidate':interval,'normalization':'8*a1*s^4 times stripe_x_unrot normalizedvolume;full plane period signature is25 times stripe signature'}
  save(HERE/'candidate_geometry.json',asset)
  # Joint planar family signatures on the exact lattice stripe[0,1/5]x[0,5].
  stage='joint_planar_stripe';atoms=load(PROJECT/'lower_round019/atoms.json');assert len(atoms)==len({(a['m'],a['n']) for a in atoms})==225
  masses=defaultdict(list)
  for f in families:masses[f['N']].append(f)
  ws={};xs={0,G//5};ys={0,5*G}
  for N,group in masses.items():
   rr=group[0]['rho_grid'];bs=[]
   for a in atoms:
    x=Q(a['x_grid'],100000)*G;y=Q(a['y_grid'],100000)*G;assert x.denominator==y.denominator==1;x,y=int(x),int(y)
    b=max(0,x-rr),min(G//5,x+rr),max(0,y-rr),min(5*G,y+rr)
    if b[0]<b[1] and b[2]<b[3]:bs.append((a['id'],b));xs.update(b[:2]);ys.update(b[2:])
   ws[N]=bs
  xs,ys=sorted(xs),sorted(ys);plane_cells=(len(xs)-1)*(len(ys)-1)
  print(json.dumps({'stage':'geometry_ready','families':len(families),'unique_masses':len(masses),'plane_segments':[len(xs)-1,len(ys)-1],'plane_cells':plane_cells,'nodes':nodes}),flush=True)
  def mask_axis(bs,q,axis):
   ev=defaultdict(list)
   for i,b in bs:ev[b[q]].append((i,True));ev[b[q+1]].append((i,False))
   bits=0;out=[]
   for lo in axis[:-1]:
    for i,opening in ev.get(lo,[]):
     if opening:bits|=1<<i
     else:bits&=~(1<<i)
    out.append(bits)
   return out
  xm={N:mask_axis(bs,0,xs) for N,bs in ws.items()};ym={N:mask_axis(bs,2,ys) for N,bs in ws.items()}
  masks_by_count={N:[sum(1<<f['family_id'] for f in group if f['n']<=c) for c in range(226)] for N,group in masses.items()}
  ph=defaultdict(int);ptrace=hashlib.sha256();Nlist=sorted(masses)
  for xi,(xl,xh) in enumerate(zip(xs,xs[1:])):
   for yi,(yl,yh) in enumerate(zip(ys,ys[1:])):
    check();mask=0
    for N in Nlist:mask|=masks_by_count[N][(xm[N][xi]&ym[N][yi]).bit_count()]
    ph[mask]+=(xh-xl)*(yh-yl);ptrace.update(f'{xi},{yi}:{mask:x}\n'.encode())
   save(HERE/'checkpoint.json',{'status':'plane_slab_saved','completed_x_slabs':xi+1,'plane_x_segments':len(xs)-1,'histogram':{str(k):v for k,v in ph.items()},'nodes':nodes,'runtime':runtime(start)})
  assert sum(ph.values())==G**2
  # Crosscheck the new symbolic masks through their exact marginal area identities.
  new_area_checks=[]
  for f in families:
   key=(f['K'],f['n']);bit=1<<f['family_id'];marginal=Q(25*sum(v for m,v in ph.items() if m&bit),G**2)
   if key in newmap:
    aa,bb,cc=map(Q,newmap[key]['normalized_exact_area_polynomial']);rho=Q(f['rho_grid'],G);expected=aa*rho*rho+bb*rho+cc;assert marginal==expected,(key,marginal,expected)
    new_area_checks.append({'K':key[0],'n':key[1],'rho_inner':str(rho),'normalized_marginal_area':str(marginal),'equals020symbolic_area':True})
   elif key in ((8,1),(8,2)):
    assert marginal==25;new_area_checks.append({'K':key[0],'n':key[1],'normalized_marginal_area':'25','equals017018wholeplane_certificate':True})
  plane={'status':'COMPLETE_joint_planar_stripe_signatures','stripe_grid':[0,G//5,0,5*G],'axes':[xs,ys],'family_count':len(families),'mass_window_boxes':{str(N):[[i,list(b)] for i,b in bs] for N,bs in ws.items()},'histogram':{str(k):v for k,v in sorted(ph.items())},'cell_trace_sha256':ptrace.hexdigest(),'cells':plane_cells,'new_symbolic_marginal_checks':new_area_checks,'period_factor_proof':'v=R(1,1)=(1/5,7/5) andw=R(3,4)=(0,5);jointcapture masks are invariant;25 translated stripes tile[0,5]^2 with y wrapping mod5','runtime':runtime(start)};save(HERE/'plane_histogram.json',plane)
  # Three-axis common-family signatures, retaining box identity until after intersection.
  stage='unrotated_3D_signatures';boxes=[];which=[]
  for f in families:
   for b in f['inner_boxes']:
    if volume(b):boxes.append(b);which.append(f['family_id'])
  axes=[sorted({0,G//2}|{b[2*q+k] for b in boxes for k in (0,1)}) for q in range(3)]
  am=[]
  for q,axis in enumerate(axes):am.append([sum(1<<i for i,b in enumerate(boxes) if b[2*q]<=a and z<=b[2*q+1]) for a,z in zip(axis,axis[1:])])
  boxfamily={};bh=defaultdict(int);btrace=hashlib.sha256()
  def family_from_boxmask(mask):
   if mask not in boxfamily:
    bits=mask;f=0
    while bits:
     bit=bits&-bits;i=bit.bit_length()-1;f|=1<<which[i];bits-=bit
    boxfamily[mask]=f
   return boxfamily[mask]
  for xi,(xl,xh) in enumerate(zip(axes[0],axes[0][1:])):
   for yi,(yl,yh) in enumerate(zip(axes[1],axes[1][1:])):
    xy=am[0][xi]&am[1][yi]
    for zi,(zl,zh) in enumerate(zip(axes[2],axes[2][1:])):
     check();mask=family_from_boxmask(xy&am[2][zi]);bh[mask]+=(xh-xl)*(yh-yl)*(zh-zl);btrace.update(f'{xi},{yi},{zi}:{mask:x}\n'.encode())
  assert sum(bh.values())==G**3//8
  save(HERE/'unrotated_histogram.json',{'status':'COMPLETE_3D_family_signatures','axes':axes,'boxes':boxes,'box_families':which,'histogram':{str(k):v for k,v in sorted(bh.items())},'cell_trace_sha256':btrace.hexdigest(),'cells':(len(axes[0])-1)*(len(axes[1])-1)*(len(axes[2])-1),'runtime':runtime(start)})
  # Intersection of signatures, not a sum of separate case areas.
  stage='five_dimensional_signature_AND';joint=defaultdict(int);total=0
  for pm,pa in ph.items():
   for bm,bv in bh.items():
    check();mask=pm&bm
    if mask:joint[mask]+=pa*bv;total+=pa*bv
  assert 0<=total<=G**5//8
  V=Q(total,G**5);shell=Q(0);shell_rows=[]
  for f in families:
   du=sum(Q(volume(ob)-volume(ib),G**3) for ib,ob in zip(f['inner_boxes'],f['outer_boxes']))
   vu=min(Q(1,8),sum((Q(volume(b),G**3) for b in f['outer_boxes']),Q(0)))
   lowrho=Q(f['rho_grid'],G);hirho=Q(f['rho_outer_interval_upper']);dp=min(Q(1),36*(hirho**2-lowrho**2));assert du>=0 and dp>=0
   err=du+vu*dp;shell+=err;shell_rows.append({'family_id':f['family_id'],'unrotated_shell_volume_upper':str(du),'outer_unrotated_union_volume_upper':str(vu),'planar_stripe_threshold_shell_area_upper':str(dp),'product_difference_upper':str(err)})
  ratio_low=8*A*sl**4*V;ratio_inner_hi=8*A*sh**4*V;ratio_hi=8*A*sh**4*min(Q(1,8),V+shell)
  C=load(PROJECT/'lower_round003/a2_equals_a3_results.json');cl,ch=map(Q,C['ratio_interval'])
  result={'status':'COMPLETE_one_new_finite_5D_union_candidate','a1':str(A),'tail_spacing':'s=(1+u)/2,u^5=2','raw_cases_accounted':len(inventory),'new18_cases_accounted':all18,'included_K_n_family_count':len(families),'normalized_inner_stripe_union_volume':str(V),'exact_inner_union_ratio_expression':f'8*(151/100)*s^4*({V})','inner_union_ratio_interval':[str(ratio_low),str(ratio_inner_hi)],'nominal_finite_family_ratio_interval':[str(ratio_low),str(ratio_hi)],'nominal_ratio_decimal_display':[dec(ratio_low),dec(ratio_hi)],'normalized_nominal_minus_inner_union_upper':str(shell),'below_C_cert_for_this_finite_candidate':ratio_hi<cl,'above_C_cert_certified':ratio_low>ch,'normalization':'periodmass25;threeunrotatedreflections give8;fullplanarmaskperiod=25*stripe;physicalJacobian=a1*s^4, so ratio=8*a1*s^4*stripevolume','all_1024_raw_families_included_before_exact_exclusions':True,'joint_union_rule':'only plane_family_mask AND unrotated_family_mask nonzero contributes, counted once regardless overlap','same_common_radius':'each family strict rationalr_lower has(2r_lower)^5<Kn;allfivecoordinates use this same radius;new18emptycases retainproof exclusions','distinct_mass':'n distinct planar lattice indices times the K canonical unrotated Cartesian indices gives Kn distinct unit atoms','boundary_guard_count':'six unrotated endpoint inequalities checked perbox;planar four peratomwindow from exact integercoordinates;signature intersection supplies onecommonfamily','shell_bound_proof':'threshold symmetric difference is contained in225 window annuli:wholeperiod area<=225*4*(rho_hi²-rho_inner²);period/stripe factor25 gives36*(...);product error<=unrotated_shell+outer_unrotated_volume*planar_shell, summed onlyaserrorbound,notunionvolume','shell_rows':shell_rows,'joint_histogram':{str(k):v for k,v in sorted(joint.items())},'event_component_indices':[interval['left_event']['event_index'],interval['right_event']['event_index']],'plane_histogram_sha256':sha(HERE/'plane_histogram.json'),'unrotated_histogram_sha256':sha(HERE/'unrotated_histogram.json'),'geometry_sha256':sha(HERE/'candidate_geometry.json'),'source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),'nodes':nodes,'runtime':runtime(start),'scope':'one genuinely new a1 finite-family union candidate only;no1946intervalscan,notallscalecomplete,notgeneralC5upper;mainownsL1/collartermination'}
  save(HERE/'results.json',result);save(HERE/'checkpoint.json',{'status':'COMPLETE_registered_candidate_only','nodes':nodes,'result_sha256':sha(HERE/'results.json'),'runtime':runtime(start),'no_next_candidate_started':True});print(json.dumps({k:v for k,v in result.items() if k not in ('joint_histogram','shell_rows','new18_cases_accounted')},indent=2),flush=True)
 except Limit as e:
  save(HERE/'limit_result.json',{'status':'CHECKPOINT_LIMIT_no_complete_ratio','reason':str(e),'nodes':nodes,'runtime':runtime(start),'no_unfinished_ratio_claim':True});print(json.dumps(load(HERE/'limit_result.json')),flush=True)
if __name__=='__main__':
 assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
 main()
