#!/usr/bin/env python3
"""UNRUN one fixed full-G scalar point through audited exact finite part."""
import sys
sys.set_int_max_str_digits(0)
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter,defaultdict
import itertools,json,hashlib,math,time,resource,gzip,traceback
HERE=Path(__file__).resolve().parent
Z=(Q(0),)*4

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def req(x,s):
 if not x:raise ArithmeticError(s)
def atomic(n,d):
 p=HERE/n;t=p.with_name(p.name+'.own-tmp')
 with t.open('w')as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))
def det(a):
 a=[list(map(Q,row))for row in a];r=Q(1)
 for j in range(len(a)):
  k=next((k for k in range(j,len(a))if a[k][j]),None)
  if k is None:return Q(0)
  if k!=j:a[k],a[j]=a[j],a[k];r=-r
  p=a[j][j];r*=p
  for k in range(j+1,len(a)):
   z=a[k][j]/p
   for l in range(j+1,len(a)):a[k][l]-=z*a[j][l]
 return r
def jet(x,k):
 if x:return (Q(1 if x>0 else-1),Q(0),Q(0),Q(0))if k==0 else Z
 if k==1:return(Q(0),Q(1),Q(0),Q(0))
 if k==3:return(Q(0),Q(0),Q(-1,3),Q(0))
 if k==5:return(Q(0),Q(0),Q(0),Q(1,5))
 return Z
def dd(nodes):
 nodes=sorted(nodes);tab=[jet(x,0)for x in nodes]
 for k in range(1,6):
  nxt=[]
  for j in range(6-k):
   nxt.append(jet(nodes[j],k)if nodes[j]==nodes[j+k]else tuple((v-u)/(nodes[j+k]-nodes[j])for u,v in zip(tab[j],tab[j+1])))
  tab=nxt
 return tab[0]
def ends(a):
 req(a.is_finite(),'finite Arb normalization');out=[]
 for z in(a.lower(),a.upper()):m,e=z.man_exp();out.append(Q(int(m))*Q(2)**int(e))
 return out

def main():
 start=time.monotonic();cpu=time.process_time();inp=json.loads((HERE/'inputs.json').read_text());reg=json.loads((HERE/'registration.json').read_text());rel=json.loads((HERE/'parent_release.json').read_text())
 req(rel.get('authorized')is True and str(rel.get('round'))=='135','one primary scalar-point release')
 req(rel['script_sha256']==sha(__file__)and rel['registration_sha256']==sha(HERE/'registration.json'),'released exact source/registration')
 for n,pin in reg['prepared_pins'].items():req(sha(HERE/n)==pin,'prepared '+n)
 for row in json.loads((HERE/'source_manifest.json').read_text())['sources']:req(sha(row['path'])==row['sha256'],'source '+row['role'])
 req(os.getpriority(os.PRIO_PROCESS,0)>=10,'nice priority at least10')
 with(HERE/'started.json').open('x')as f:json.dump({'pid':os.getpid(),'release_sha256':sha(HERE/'parent_release.json'),'unix_time':time.time()},f);f.write('\n');f.flush();os.fsync(f.fileno())
 status='FAILED_NO_RETRY';phase='true_pair_qualification';failure=None
 try:
  point=list(map(Q,inp['point']));req(point==[Q(1),Q(1),Q(1),Q(3,8),Q(1,4),Q(1,8)],'only registered point')
  req(max(map(abs,point))==1 and point==sorted(point,reverse=True),'true ordered cube boundary point')
  V=[(1,)+z for z in itertools.product((-1,1),repeat=5)];basis=[(1,)*6]+[tuple(-1 if j==i else 1 for j in range(6))for i in range(1,6)]
  req(abs(det(basis))==32,'rank6 original-body basis witness')
  pairrows=[];X=dot(point,point)
  for i in range(32):
   for j in range(i+1,32):
    gii=dot(V[i],V[i]);gjj=dot(V[j],V[j]);gij=dot(V[i],V[j]);D=gii*gjj-gij*gij;ui=dot(V[i],point);uj=dot(V[j],point)
    req(D>=0,'honest pair Gram');rho=X-(gjj*ui*ui-2*gij*ui*uj+gii*uj*uj)/D if D else X-ui*ui/gii
    req(rho>=0,'nonnegative exact pair residual');pairrows.append({'i':i,'j':j,'Gram_determinant':str(D),'distance_squared':str(rho)})
  req(len(pairrows)==496,'all actual unordered pairs')
  rho=min(Q(r['distance_squared'])for r in pairrows)
  atomic('pair_qualification.json',{'point':inp['point'],'body_columns':V,'basis_witness':basis,'pair_count':496,'minimum_squared_distance':str(rho),'pairs':pairrows,'off_true_S2':rho>0})
  req(rho>0,'S2_REJECT_INFINITE_GREEN_no_finite_part_value')
  phase='frozen_complete_fan_binding';fan=json.loads(Path(inp['fan_path']).read_text());audit=json.loads(Path(inp['fan_audit_path']).read_text())
  req(audit['status']=='PASS'and audit['input_sha256']==sha(inp['fan_path']),'independently audited exact same fan')
  cells=fan['cells'];req(len(cells)==21 and audit['full_cells']==21,'exact original21 fan')
  req(Q(audit['simplex_section_volume_sum'])==Q(1,120),'complete chamber volume from audited source')
  req(Q(audit['full_exp_minus_symbol_integral'])==Q(fan['full_exp_minus_symbol_integral'])==Q(inp['R0']),'same original normalization')
  cones=[]
  for cell in cells:
   rays=[tuple(map(Q,r))for r in cell['rays_a0_equals_one']];rates=list(map(Q,cell['rates']));D=Q(cell['determinant'])
   req(len(rays)==6 and all(len(r)==6 for r in rays)and abs(det(rays))==D>0,'bound source simplex determinant')
   req(all(q==sum(abs(dot(v,r))for v in V)for q,r in zip(rates,rays)),'original psi rates')
   req(all(min(dot(v,r)for r in rays)>=0 or max(dot(v,r)for r in rays)<=0 for v in V),'all body signs linear on simplex')
   cones.append((rays,rates,D/math.prod(rates)))
  req(46080*sum(c[2]for c in cones)==Q(inp['R0']),'precise full-group multiplicity normalization')
  phase='one_point_exact_orbit_and_Laurent';ints=inp['point_integer_numerators'];den=inp['point_denominator'];req([Q(v,den)for v in ints]==point,'exact integer point binding')
  orbit=Counter()
  for p in itertools.permutations(range(6)):
   for mask in range(64):orbit[tuple((-1 if(mask>>k)&1 else 1)*ints[p[k]]for k in range(6))]+=1
  req(sum(orbit.values())==46080 and len(orbit)==7680 and set(orbit.values())=={6},'complete actual orbit and stabilizer')
  merged=defaultdict(Q)
  for transformed,multiplicity in sorted(orbit.items()):
   y=[Q(v,den)for v in transformed]
   for rays,rates,alpha in cones:
    nodes=tuple(sorted(dot(y,r)/q for r,q in zip(rays,rates)));merged[nodes]+=multiplicity*alpha
  req(len(merged)<=161280 and sum(merged.values())==Q(inp['R0']),'complete exact node aggregate')
  total=list(Z);saved=[];zeros=repeated=0
  for index,(nodes,weight)in enumerate(sorted(merged.items())):
   ch=dd(nodes);total=[a+weight*b for a,b in zip(total,ch)];zeros+=int(0 in nodes);repeated+=int(len(set(nodes))<6)
   saved.append({'nodes':[str(v)for v in nodes],'weight':str(weight),'channels':[str(v)for v in ch]})
   if(index+1)%10000==0:atomic('progress.json',{'phase':phase,'merged_tuples_done':index+1,'merged_tuples_total':len(merged)})
  with(HERE/'exact_tuple_certificate.json.gz').open('xb')as raw:
   with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as f:f.write(json.dumps({'channel_order':['Gamma','epsilon^-1','epsilon^-3','epsilon^-5'],'tuples':saved},separators=(',',':')).encode())
  req(all(z==0 for z in total[1:]),'complete exact negative-channel cancellation')
  req(total[0]>0,'positive complete Green coefficient')
  from flint import arb,fmpq,ctx
  ctx.prec=512;q=total[0]/128;value=arb(fmpq(q.numerator,q.denominator))/arb.pi()**5;lo,hi=ends(value);S=1<<384
  lo=Q((lo*S).numerator//(lo*S).denominator,S);hi=Q(-((-hi*S).numerator//(-hi*S).denominator),S)
  req(0<lo<=hi,'strict complete point enclosure');target=Q(inp['target'])
  saved128=json.loads(Path(inp['round128_certificate_path']).read_text());M6=Q(saved128['M_power6_upper']);A=Q(8,9)**6*24576**5;gap=M6-A*hi**6
  certificate_budget_at_GU_power5=Q(9,8)**6*M6/(384**5*hi**6)
  req((gap>=0)==(certificate_budget_at_GU_power5>=64**5),'exact saved-upper certificate budget equivalence')
  status='COMPLETE_ONE_FULL_GREEN_POINT'
  atomic('results.json',{'status':status,'point':inp['point'],'off_true_S2':True,'minimum_pair_squared_distance':str(rho),'orbit_size':7680,'stabilizer':6,'uncompressed_group_cone_terms':967680,'compressed_orbit_cone_terms':161280,'merged_exact_node_tuples':len(merged),'zero_node_tuples':zeros,'repeated_node_tuples':repeated,'complete_channels':[str(v)for v in total],'Q_over_pi5':str(q),'G_interval':[str(lo),str(hi)],'G_display':[format(float(lo),'.17g'),format(float(hi),'.17g')],'target':str(target),'G_upper_below_floor_target':hi<target,'G_lower_at_least_target':lo>=target,'essential_floor_target_excluded_by_offS2_continuity':hi<target,'saved128_M_power6_upper':str(M6),'saved128_M6_minus_threshold_at_GpointU':str(gap),'saved128_certificate_budget_at_GpointU_power5':str(certificate_budget_at_GU_power5),'saved128_upper_certificate_cannot_reach_below64_for_any_valid_g':gap>=0,'actual_M_or_model_cost_or_C6_lower_bound':False,'target_minus_G_upper_exact':str(target-hi),'G_lower_minus_target_exact':str(lo-target),'full_positive_time_G':True,'numerical_time_truncation':False,'unpaid_time_tail':False,'global_G_floor_proved':False,'new_C6_bound':False,'tuple_certificate_sha256':sha(HERE/'exact_tuple_certificate.json.gz')})
  print(json.dumps({'status':status,'G_interval_display':[format(float(lo),'.17g'),format(float(hi),'.17g')],'G_upper_below_floor_target':hi<target,'merged_tuples':len(merged)}),flush=True)
 except BaseException as e:
  failure=type(e).__name__+': '+str(e);atomic('failure.json',{'phase':phase,'failure':failure,'traceback':traceback.format_exc()});raise
 finally:atomic('run_receipt.json',{'status':status,'phase':phase,'failure':failure,'pid':os.getpid(),'wall_before_receipt':time.monotonic()-start,'CPU_before_receipt':time.process_time()-cpu,'science_cores':1,'priority':os.getpriority(os.PRIO_PROCESS,0),'children':0,'time_or_CPU_ceiling':None,'automatic_retry':False,'Arb_bits':512,'new_points':1,'MC_native_quad_MATLAB_calls':0})
if __name__=='__main__':main()
