"""107 PREPARED ONLY: two fixed configurations, no run before parent release.

Latest user policy removes runtime/CPU caps. Inputs, geometry,1024-leaf/
2040-cell-bound/depth20 limits remain finite. One process; no old replay.
"""
import sys
sys.dont_write_bytecode=True
import datetime,hashlib,heapq,itertools,json,os,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_json(name,data):
 f=P/name;tmp=f.with_name(f.name+'.tmp')
 tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(f)
release_path=P/'parent_release.json'
if not release_path.exists():raise SystemExit('107 PREPARED_NOT_RELEASED: no mathematical execution')
release=json.loads(release_path.read_text())
if not(release.get('authorized')is True and release.get('round')=='107' and release.get('script_sha256')==sha(Path(__file__)) and release.get('registration_sha256')==sha(P/'registration.json')):
 raise SystemExit('107 release missing or mismatched')
if any((P/n).exists()for n in ['started.json','results.json','run_receipt.json','checkpoint.json','pair_geometry_certificate.json','cell_certificate.json']):
 raise SystemExit('107 one finite attempt only; existing runtime artifact')
started=time.monotonic();start_cpu=time.process_time();start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat();completed=0
write_json('started.json',{'round':'107','pid':os.getpid(),'started_utc':start_utc,'runtime_policy':'unlimited user policy effective2026-10-26T16:43:45.633722UTC','worker_processes':1,'children':0,'release_sha256':sha(release_path)})

def run_math():
 from flint import arb,fmpq,ctx
 ctx.prec=384
 reg=json.loads((P/'registration.json').read_text());up=json.loads((P/'upstream_constants.json').read_text())
 assert sha(P/'upstream_constants.json')==reg['upstream_constants_sha256']
 assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
 assert reg['configuration_cap']==2 and reg['geometry']['k']==1
 assert reg['benchmark_parameters']['old_c_squared']=='1024/3'
 assert reg['cell_limits']=={'max_leaves':1024,'max_cell_bounds':2040,'max_total_depth':20,'initial_leaves':8}
 def A(x):
  x=Q(x);return arb(fmpq(x.numerator,x.denominator))
 def endpoint(x,side):
  m,e=(x.lower()if side=='lower'else x.upper()).man_exp()
  return Q(int(m))*Q(2)**int(e)
 def enc(x):
  l,e=x.lower().man_exp();u,f=x.upper().man_exp();e=int(e);f=int(f);g=min(e,f)
  return {'lower_mantissa':str(int(l)*2**(e-g)),'upper_mantissa':str(int(u)*2**(f-g)),'common_binary_exponent':str(g)}
 def saved(x,side):return Q(int(x[side+'_mantissa']))*Q(2)**int(x['common_binary_exponent'])
 def poly(p,u):return p[0]+p[1]*u+p[2]*u*u
 def quad_range(p,a,b):
  points=[a,b]
  if p[2]:
   v=-p[1]/(2*p[2])
   if a<v<b:points.append(v)
  values=[(poly(p,u),u)for u in points]
  lo=min(z[0]for z in values);hi=max(z[0]for z in values)
  return lo,hi,[u for z,u in values if z==lo],[u for z,u in values if z==hi]
 def interval_strings(lo,hi):return {'lower':str(lo),'upper':str(hi)}
 configs=[];delta=Q(1,64);u_start=Q(1,4);u_end=Q(3,4);q=Q(6,5)

 # Configuration1: new exact algebraic(P,J2) aggregation; no kernel/peak.
 V=[(1,)+tuple(z)for z in itertools.product([-1,1],repeat=5)]
 assert len(V)==32 and V[0]==(1,-1,-1,-1,-1,-1)and V[31]==(1,1,1,1,1,1)
 grouped={};branch_data={};coverage={}
 for label,excluded in [('v',31),('w',0)]:
  table={};digest=hashlib.sha256();visits=0
  for i in range(32):
   if i==excluded:continue
   for j in range(32):
    if j==excluded or j==i:continue
    Bi=sum(V[i][1:]);Bj=sum(V[j][1:]);c=sum(a*b for a,b in zip(V[i],V[j]));D=36-c*c
    assert D in(20,32,36)
    p=(Q(1)-Q(12-2*c,D),-Q((12-2*c)*(Bi+Bj),D),Q(5)-Q(6*(Bi*Bi+Bj*Bj)-2*c*Bi*Bj,D))
    key=(p,D);table[key]=table.get(key,0)+1;visits+=1
    digest.update((f'{i},{j};{p[0]},{p[1]},{p[2]};{D}\n').encode())
  assert visits==930 and sum(table.values())==930 and len(table)<=105
  hist={D:sum(n for(p,d),n in table.items()if d==D)for D in(20,32,36)}
  assert hist=={20:180,32:450,36:300}
  rows=[];internal=[]
  for index,((p,D),n)in enumerate(sorted(table.items(),key=lambda x:(x[0][0],x[0][1]))):
   lo,hi,atlo,athi=quad_range(p,u_start,u_end)
   assert lo>delta*delta
   gap=A(lo).sqrt()-A(delta);assert gap.lower()>0
   name=f'{label}{index:03d}'
   rows.append({'class':name,'P_coefficients':[str(z)for z in p],'J2_squared':D,'ordered_multiplicity':n,'P_min_exact':str(lo),'P_max_exact':str(hi),'minimum_nodes':[str(z)for z in atlo],'maximum_nodes':[str(z)for z in athi],'normal_radius_gap_at_full_delta':enc(gap)})
   internal.append((name,p,D,n))
  grouped[label]=internal
  branch_data[label]={'excluded_unoriented_index':excluded,'pair_visits':visits,'ordered_count_by_J2_squared':hist,'algebraic_class_count':len(rows),'pair_visit_digest':digest.hexdigest(),'classes':rows}
 geometry={'status':'COMPLETE_TWO_BRANCH_K1_ALGEBRAIC_PAIR_AGGREGATION','frame_order':[list(v)for v in V],'first_face_coordinate':1,'k':1,'u_domain':['1/4','3/4'],'normal_radius':'1/64','grouping':'equal exact quadraticP andJ2 only; not asserted geometric orbits','branches':branch_data}
 write_json('pair_geometry_certificate.json',geometry)
 configs.append({'configuration':1,'kind':'k1 exact(P,J2,n) aggregation and every-class positivity','branches':{b:{k:v for k,v in row.items()if k!='classes'}for b,row in branch_data.items()},'geometry_certificate_sha256':sha(P/'pair_geometry_certificate.json'),'old_tuple_kernel_peak_replays':0})
 global completed
 completed=1
 write_json('checkpoint.json',{'round':'107','status':'CONFIGURATION1_COMPLETE','completed_configurations':1,'configuration1':configs[0],'elapsed_seconds':time.monotonic()-started})
 print(json.dumps({'stage':'geometry_complete','classes':{b:len(grouped[b])for b in grouped},'pair_visits_each':930}),flush=True)

 # Configuration2: coefficient-free full mixed q integral; both comparisons
 # use this ONE shape enclosure, never refresh the cells for a coefficient.
 pi=arb.pi();v4=pi*pi/2
 ap=1/(2*A(5).sqrt());aplo=endpoint(ap,'lower');aphi=endpoint(ap,'upper')
 mass=A(Q(10,7))*v4*A(5).sqrt()*A(delta)**A(Q(14,5))
 masslo=endpoint(mass,'lower');masshi=endpoint(mass,'upper')
 Kold=3*A(3).sqrt()/(64*pi**4);Kref=1/(8*pi**5)
 oldscale=Kold**A(q);refscale=Kref**A(q)
 scales={'old098_valid':(endpoint(oldscale,'lower'),endpoint(oldscale,'upper')),'ref108_conditional_shape_only':(endpoint(refscale,'lower'),endpoint(refscale,'upper'))}
 target=64*A(Q(1,65000)/6)**A(q);paid=Q(up['paid_CEH_I_upper'])
 remaining_lo=endpoint(target,'lower')-paid;remaining_hi=endpoint(target,'upper')-paid
 old104_lo=saved(up['old104_k1_weighted_I'],'lower');old104_hi=saved(up['old104_k1_weighted_I'],'upper')
 sqrtDelta={D:A(D).sqrt()for D in(20,32,36)};sqrt56=A(Q(5,6)).sqrt()
 Fv_inside=(Q(5,6),-Q(5,3),Q(5,6))
 Fv_boundary=(Q(3,2),-Q(5),Q(5))
 Fw_inside=(Q(5,6),Q(5,3),Q(5,6))
 Fw_boundary=(Q(147,128),-Q(5,8),Q(5))
 # Fw_boundary=(17/16)^2+5(u-1/16)^2=147/128-5u/8+5u^2.
 assert Fw_boundary[0]==Q(17,16)**2+5*Q(1,16)**2
 bounds_done=0
 def origin_range(label,a,b):
  if label=='v':
   if b<=Q(2,5):p=Fv_inside
   else:assert a>=Q(2,5);p=Fv_boundary
  else:
   if b<=Q(11,40):p=Fw_inside
   else:assert a>=Q(11,40);p=Fw_boundary
  lo,hi,_,_=quad_range(p,a,b);return lo,hi
 def branch_range(label,a,b,rlo,rhi):
  Wlo_ball=A(0);Whi_ball=A(0);gapmins=[]
  for name,p,D,n in grouped[label]:
   pmin,pmax,_,_=quad_range(p,a,b)
   denlo=A(pmin).sqrt()-A(rhi);denhi=A(pmax).sqrt()-A(rlo)
   assert denlo.lower()>0 and denhi.lower()>0
   gapmins.append(endpoint(denlo,'lower'))
   Whi_ball+=A(n)/(sqrtDelta[D]*denlo)
   Wlo_ball+=A(n)/(sqrtDelta[D]*denhi)
  Wlo=endpoint(Wlo_ball,'lower');Whi=endpoint(Whi_ball,'upper')
  fmin,fmax=origin_range(label,a,b)
  h0lo=endpoint(A(fmin).sqrt()-A(rhi),'lower');h0hi=endpoint(A(fmax).sqrt()-A(rlo),'upper')
  if label=='v':
   if b<=Q(1,2):mlo,mhi=a,b
   else:assert a>=Q(1,2);mlo,mhi=1-b,1-a
  else:mlo,mhi=Q(7,8)-b,Q(7,8)-a
  h1lo=endpoint(sqrt56*A(mlo)-A(rhi),'lower');h1hi=endpoint(sqrt56*A(mhi)-A(rlo),'upper')
  assert min(h0lo,h1lo)>0
  lower=endpoint(A(Wlo)/A(h0hi*h0hi*h1hi*h1hi),'lower')
  upper=endpoint(A(Whi)/A(h0lo*h0lo*h1lo*h1lo),'upper')
  return lower,upper,min(gapmins+[h0lo,h1lo])
 def cell(a,b,s0,s1,du,ds,cid,seed):
  nonlocal bounds_done
  assert bounds_done<reg['cell_limits']['max_cell_bounds'];bounds_done+=1
  rlo=delta*s0**5;rhi=delta*s1**5
  vl,vu,vg=branch_range('v',a,b,rlo,rhi);wl,wu,wg=branch_range('w',a,b,rlo,rhi)
  tl=aplo+rlo*(vl+wl);tu=aphi+rhi*(vu+wu)
  powerlo=endpoint(A(tl)**A(q),'lower');powerhi=endpoint(A(tu)**A(q),'upper')
  assert 0<powerlo and powerlo**5<=tl**6 and powerhi**5>=tu**6
  weight=(b-a)*(s1**14-s0**14)
  lo=masslo*weight*powerlo;hi=masshi*weight*powerhi
  assert 0<=lo<=hi
  return {'id':cid,'seed':seed,'a':a,'b':b,'s0':s0,'s1':s1,'u_depth':du,'s_depth':ds,'lower':lo,'upper':hi,'T_lower':tl,'T_upper':tu,'power_lower':powerlo,'power_upper':powerhi,'weight':weight,'positive_gap_lower':min(vg,wg)}
 def serial(c):
  return {'id':c['id'],'seed':c['seed'],'u':[str(c['a']),str(c['b'])],'s':[str(c['s0']),str(c['s1'])],'r':[str(delta*c['s0']**5),str(delta*c['s1']**5)],'depths':[c['u_depth'],c['s_depth']],'regular_shape_T':interval_strings(c['T_lower'],c['T_upper']),'regular_shape_q_power':interval_strings(c['power_lower'],c['power_upper']),'exact_weight':str(c['weight']),'weighted_shape_integral':interval_strings(c['lower'],c['upper']),'positive_denominator_lower':str(c['positive_gap_lower'])}
 def benchmark(L,U):
  out={}
  for name,(sl,su)in scales.items():
   lo=sl*L;hi=su*U
   status=('FIXED_EXPRESSION_OBSTRUCTION'if lo>remaining_hi else'FITS_K1_PARTIAL_ALLOCATION'if hi<remaining_lo else'UNSEPARATED')
   out[name]={'I':interval_strings(lo,hi),'S_plan_other_regions_zero':interval_strings(endpoint(target,'lower')-paid-hi,endpoint(target,'upper')-paid-lo),'allocation_status':status,'coefficient_validity':'accepted098 interface'if name=='old098_valid'else'108_NOT_ACCEPTED_CONDITIONAL_SHAPE_EXPRESSION_ONLY'}
  return out
 ubreak=[Q(1,4),Q(11,40),Q(2,5),Q(1,2),Q(3,4)];sbreak=[Q(0),Q(1,2),Q(1)]
 leaves={};heap=[];seeds=[];transitions=[];nextid=0
 for ua,ub in zip(ubreak,ubreak[1:]):
  for sa,sb in zip(sbreak,sbreak[1:]):
   cid=nextid;nextid+=1;c=cell(ua,ub,sa,sb,0,0,cid,cid)
   leaves[cid]=c;heapq.heappush(heap,(-(c['upper']-c['lower']),cid));seeds.append({'id':cid,'u':[str(ua),str(ub)],'s':[str(sa),str(sb)]})
 assert len(leaves)==8
 L=sum((c['lower']for c in leaves.values()),Q(0));U=sum((c['upper']for c in leaves.values()),Q(0));splits=0;outcome=None
 while True:
  benchmarks=benchmark(L,U)
  if all(x['allocation_status']=='FIXED_EXPRESSION_OBSTRUCTION'for x in benchmarks.values()):outcome='BOTH_FIXED_EXPRESSIONS_OBSTRUCTED_STOP';break
  if len(leaves)>=1024 or bounds_done+2>2040:outcome='FINITE_NODE_CAP_REACHED';break
  chosen=None
  while heap:
   _,cid=heapq.heappop(heap)
   if cid not in leaves:continue
   parent=leaves[cid]
   if parent['u_depth']+parent['s_depth']>=20:continue
   chosen=parent;break
  if chosen is None:outcome='FINITE_DEPTH_CAP_REACHED';break
  a,b,s0,s1=chosen['a'],chosen['b'],chosen['s0'],chosen['s1'];du,ds=chosen['u_depth'],chosen['s_depth']
  axis='u' if du<=ds else 's';newids=[nextid,nextid+1];nextid+=2
  if axis=='u':
   mid=(a+b)/2;children=[cell(a,mid,s0,s1,du+1,ds,newids[0],chosen['seed']),cell(mid,b,s0,s1,du+1,ds,newids[1],chosen['seed'])]
  else:
   mid=(s0+s1)/2;children=[cell(a,b,s0,mid,du,ds+1,newids[0],chosen['seed']),cell(a,b,mid,s1,du,ds+1,newids[1],chosen['seed'])]
  del leaves[chosen['id']]
  for cnew in children:leaves[cnew['id']]=cnew;heapq.heappush(heap,(-(cnew['upper']-cnew['lower']),cnew['id']))
  L+=sum((x['lower']for x in children),Q(0))-chosen['lower'];U+=sum((x['upper']for x in children),Q(0))-chosen['upper']
  transitions.append({'parent':chosen['id'],'axis':axis,'children':newids});splits+=1
  if splits%64==0:
   summary={'round':'107','status':'CONFIGURATION2_RUNNING_FINITE_BATCH','completed_configurations':1,'leaves':len(leaves),'cell_bounds':bounds_done,'shape_integral':interval_strings(L,U),'benchmarks':benchmark(L,U),'elapsed_seconds':time.monotonic()-started}
   write_json('checkpoint.json',summary);print(json.dumps({'stage':'cell_progress','leaves':len(leaves),'cell_bounds':bounds_done,'elapsed_seconds':summary['elapsed_seconds']}),flush=True)
 benchmarks=benchmark(L,U)
 assert L==sum((x['lower']for x in leaves.values()),Q(0))and U==sum((x['upper']for x in leaves.values()),Q(0))
 assert len(leaves)<=1024 and bounds_done<=2040
 certificate={'status':'COMPLETE_FINITE_K1_SHAPE_INTEGRAL_ENCLOSURE','k':1,'u_clamp_partition':[str(x)for x in ubreak],'substitution':'r=delta*s^5','delta':'1/64','q':'6/5','shape_definition':'J=20v4sqrt5delta^(14/5)integral du integral ds s^13[1/J_vw+r(Wv/Hv+Ww/Hw)]^(6/5)','mass_constant_after_s13_exact_integration':interval_strings(masslo,masshi),'large_shape_coefficient_1_over_Jvw':interval_strings(aplo,aphi),'initial_seeds':seeds,'refinements':transitions,'leaves':[serial(leaves[i])for i in sorted(leaves)],'shape_integral':interval_strings(L,U),'geometry_sha256':sha(P/'pair_geometry_certificate.json'),'parent_child_cover_not_merely_volume_sum':True}
 write_json('cell_certificate.json',certificate)
 oldI=benchmarks['old098_valid']['I'];strict_improvement=Q(oldI['upper'])<old104_lo
 configs.append({'configuration':2,'kind':'complete k1 mixed-q spatial budget with real pair geometry','stopping_outcome':outcome,'leaves':len(leaves),'cell_bounds':bounds_done,'max_total_depth_seen':max(x['u_depth']+x['s_depth']for x in leaves.values()),'shape_integral':interval_strings(L,U),'coefficient_independent_cells':True,'Kold':enc(Kold),'Kref_conditional_only':enc(Kref),'Kref_over_Kold':enc(Kref/Kold),'scale_old':enc(oldscale),'scale_ref_conditional_only':enc(refscale),'benchmarks':benchmarks,'saved104_k1_I_interval':interval_strings(old104_lo,old104_hi),'old098_strictly_improves_saved104_k1_upper_expression':strict_improvement,'old098_upper_over_saved104_lower_exact':str(Q(oldI['upper'])/old104_lo),'paid_CEH_I_upper_exact':str(paid),'conditional_remaining_I_interval':interval_strings(remaining_lo,remaining_hi),'global_g_I_S_status':'UNPAID: all other middle, far and globalg_L; core/k1 diagnostic not fullS','cell_certificate_sha256':sha(P/'cell_certificate.json'),'108_theorem_used':False,'extra_coefficient_or_cell_refreshes':0})
 completed=2
 return {'status':'COMPLETE_TWO_FINITE_CONFIGURATIONS_PROOF_REVIEW_SEPARATE','round':'107','configurations':configs,'precision_bits':384,'script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json'),'upstream_constants_sha256':sha(P/'upstream_constants.json'),'new_native_kernel_peak_old_grid_tuple_replays':0,'worker_processes':1,'children':0,'runtime_limit':None,'finite_parameter_scope':'k1 only, fullu/r, two registered benchmarks one shape enclosure','105_or108_unaccepted_bound_applied_as_theorem':False,'Cn_improved':False}

status='FIRST_FAILURE_STOP';failure=None
try:
 result=run_math();write_json('results.json',result);status='COMPLETED'
except KeyboardInterrupt:status='INTERRUPTED_NO_AUTOMATIC_RETRY';failure=traceback.format_exc()
except Exception:failure=traceback.format_exc()
finally:
 receipt={'status':status,'round':'107','pid':os.getpid(),'started_utc':start_utc,'receipt_constructed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_before_final_receipt_write':time.monotonic()-started,'CPU_process_total_before_final_receipt_write':time.process_time(),'completed_configurations':completed,'runtime_policy':'unlimited20261027; prior177/CPU60 superseded before first run','finite_nodes_not_refreshed':True,'worker_processes':1,'children':0,'retries':0,'failure':failure,'release_sha256':sha(release_path)}
 write_json('run_receipt.json',receipt);print(json.dumps({'status':status,'completed_configurations':completed,'pid':os.getpid(),'elapsed_before_receipt':receipt['elapsed_before_final_receipt_write']}),flush=True)
if status!='COMPLETED':raise SystemExit(1)
