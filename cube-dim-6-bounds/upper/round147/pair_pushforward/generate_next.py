#!/usr/bin/env python3
"""PREPARED ONLY: exact m3/m6 density compiler; requires main exact-hash release."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import hashlib,json,os,sys,time,traceback
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(v,m):
 if not v:raise ArithmeticError(m)
def write(n,v):
 with(P/n).open('x')as f:json.dump(v,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def det(a,b):return a[0]*b[1]-a[1]*b[0]
def direction(v):return tuple(x/sum(v)for x in v)
def release_requirements():
 return {'round':'147_NEXT_M3_M6','authorized':True,'scope':'TWO_FIXED_REAL_PAIR_EXACT_DENSITY_GENERATIONS_NO_K_INTEGRAL','workers':1,'nice_minimum':10,'time_limit':None,'max_real_cases':2,'max_density_stage_checks':7,'max_nonzero_terms_per_stage':10000,'max_raw_endpoint_contributions':100000,'new_k_integral_kernel_native_calls_permitted':0,'run_sha256':sha(__file__),'case_packet_sha256':sha(P/'next_cases.json'),'proof_sha256':sha(P/'PROOF.md'),'registration_sha256':sha(P/'next_registration.json')}
def preflight():
 reg=read(P/'next_registration.json');packet=read(P/'next_cases.json')
 need(sha(__file__)==reg['run_sha256']and sha(P/'PROOF.md')==reg['proof_sha256']and sha(P/'next_cases.json')==reg['case_packet_sha256'],'FIXED_NEXT_PACKET')
 for p,h in packet['source_hashes'].items():need(sha(p)==h,'SOURCE_HASH')
 cones=read(packet['cone_table_path'])['cones'];need(len(packet['cases'])==2,'TWO_CASES_ONLY')
 for row,m in zip(packet['cases'],[3,6]):
  cone=cones[row['cone_id']];i,j=row['body_pair_zero_based'];cols=[list(v)for v in zip(cone['nodes_by_body'][i],cone['nodes_by_body'][j])]
  need(cols==row['all_six_columns']and sum(any(F(z)for z in v)for v in cols)==m,'EXACT_REAL_PATTERN_IDENTITY')
 return packet,reg

def run():
 rel=read(P/'next_parent_release.json')
 for k,v in release_requirements().items():need(rel.get(k)==v,'RELEASE:'+k)
 need(rel.get('main_approved')is True and rel.get('independent_review_approved')is True,'DUAL_GATE_REQUIRED')
 need(not(P/'next_started.json').exists(),'NO_AUTOMATIC_RESTART');packet,reg=preflight()
 os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
 need(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']),'ONE_THREAD_ENV')
 write('next_started.json',{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'release_sha256':sha(P/'next_parent_release.json'),'time_limit':None})
 start=time.monotonic();cpu=time.process_time();counters={'density_stage_checks':0,'raw_endpoint_contributions':0,'new_k_integral_kernel_native_calls':0,'old_6D_grid_replays':0}
 try:
  import sympy as sp
  need(sp.__version__==reg['sympy_version'],'PINNED_SYMPY')
  X,Y,R,U,V,S,T=sp.symbols('X Y R U V S T')
  def q(x):return sp.Rational(x.numerator,x.denominator)
  def linear(v):return q(v[0])*X+q(v[1])*Y
  def polynomial(v):
   p=sp.Poly(sp.expand(v),X,Y,domain=sp.QQ);return p.as_expr()
  def interval_for_sector(left,right,v,sample):
   lowers=[(F(0),F(0))];uppers=[]
   for L in [(-left[1],left[0]),(right[1],-right[0])]:
    d=dot(L,v);value=dot(L,sample)
    if not d:
     if value<0:return None
    elif d>0:uppers.append(tuple(z/d for z in L))
    else:lowers.append(tuple(z/d for z in L))
   need(bool(uppers),'POINTED_NONNEGATIVE_SUPPORT_GIVES_FINITE_UPPER')
   lo=max(lowers,key=lambda z:dot(z,sample));hi=min(uppers,key=lambda z:dot(z,sample))
   return None if dot(hi,sample)<=dot(lo,sample)else(lo,hi)
  def verify_stage(cells,active,stage):
   counters['density_stage_checks']+=1;need(counters['density_stage_checks']<=7,'STAGE_CAP')
   need(len(cells)<=len(active)-1,'AT_MOST_M_MINUS_ONE_SECTORS')
   count=sum(len(terms)for _,_,terms in cells);need(count<=10000,'TERM_CAP');total=sp.S(0);encoded=[]
   for left,right,terms in cells:
    need(det(left,right)>0,'ORIENTED_SECTOR');D=det(left,right);out=[]
    for lam,poly in sorted(terms.items()):
     need(sp.Poly(poly,X,Y).total_degree()<=len(active)-2,'DEGREE_AT_MOST_M_MINUS_TWO')
     ab=sp.Poly(sp.expand(poly.subs({X:q(left[0])*U+q(right[0])*V,Y:q(left[1])*U+q(right[1])*V},simultaneous=True)),U,V,domain=sp.QQ)
     L=(q(lam[0])+S)*q(left[0])+(q(lam[1])+T)*q(left[1]);W=(q(lam[0])+S)*q(right[0])+(q(lam[1])+T)*q(right[1])
     for (i,j),coeff in ab.terms():total+=q(D)*coeff*factorial(i)*factorial(j)/(L**(i+1)*W**(j+1))
     out.append({'lambda':list(map(str,lam)),'polynomial':str(poly)})
    encoded.append({'left_ray':list(map(str,left)),'right_ray':list(map(str,right)),'terms':out})
   target=sp.prod(1/(1+q(v[0])*S+q(v[1])*T)for v in active)
   need(sp.cancel(total-target)==0,'EXACT_BIVARIATE_LAPLACE_NUMERATOR_ZERO')
   normal=sp.cancel(total);need(sp.simplify(normal.subs({S:0,T:0}))==1,'WHOLE_RATIONAL_NORMALIZATION_ONE')
   return {'stage':stage,'active_columns':len(active),'open_sector_count':len(cells),'nonzero_term_count':count,'density_cells':encoded,'joint_Laplace_product':str(target),'exact_Laplace_identity':True,'normalization':1,'positivity_basis':'positive exponential-ray convolution theorem in PROOF.md; signed terms not separately positive'}
  outputs=[]
  for case in packet['cases']:
   allcols=[tuple(map(F,v))for v in case['all_six_columns']];cols=[v for v in allcols if any(v)];need(len(cols)in (3,6),'REGISTERED_ACTIVE_COUNTS')
   pivot=next((i,j)for i in range(len(cols))for j in range(i+1,len(cols))if det(cols[i],cols[j]))
   u,w=cols[pivot[0]],cols[pivot[1]]
   if det(u,w)<0:u,w=w,u
   D=det(u,w);lam=((w[1]-u[1])/D,(u[0]-w[0])/D);left,right=direction(u),direction(w)
   need(dot(lam,u)==dot(lam,w)==1,'INITIAL_COVECTOR')
   fan=[left,right];cells=[(left,right,{lam:q(1/D)})];active=[u,w];stages=[verify_stage(cells,active,2)]
   remaining=[cols[j]for j in range(len(cols))if j not in pivot]
   for v in remaining:
    newfan=sorted(set(fan+[direction(v)]),key=lambda z:z[1]);newcells=[]
    for a,b in zip(newfan,newfan[1:]):
     sample=(a[0]+b[0],a[1]+b[1]);newterms={}
     for oldleft,oldright,terms in cells:
      bounds=interval_for_sector(oldleft,oldright,v,sample)
      if bounds is None:continue
      lo,hi=bounds
      for lam,poly in terms.items():
       eta=dot(lam,v)-1;pr=sp.Poly(sp.expand(poly.subs({X:X-q(v[0])*R,Y:Y-q(v[1])*R},simultaneous=True)),R)
       if not eta:
        primitive=sum(coef*R**(deg[0]+1)/sp.Integer(deg[0]+1)for deg,coef in pr.terms());value=polynomial(primitive.subs(R,linear(hi))-primitive.subs(R,linear(lo)))
        newterms[lam]=polynomial(newterms.get(lam,0)+value);counters['raw_endpoint_contributions']+=1
       else:
        expr=pr.as_expr();Q=sum((-1)**k*sp.diff(expr,R,k)/q(eta)**(k+1)for k in range(pr.degree()+1))
        for bound,sign in [(hi,1),(lo,-1)]:
         newlam=tuple(lam[j]-eta*bound[j]for j in range(2));value=polynomial(sign*Q.subs(R,linear(bound)))
         newterms[newlam]=polynomial(newterms.get(newlam,0)+value);counters['raw_endpoint_contributions']+=1
       need(counters['raw_endpoint_contributions']<=100000,'RAW_ENDPOINT_CONTRIBUTION_CAP')
     newterms={lam:poly for lam,poly in newterms.items()if poly!=0}
     if newterms:newcells.append((a,b,newterms))
    active.append(v);fan=newfan;cells=newcells;stages.append(verify_stage(cells,active,len(active)))
   outputs.append({'case_id':case['case_id'],'cone_id':case['cone_id'],'body_pair_zero_based':case['body_pair_zero_based'],'all_six_columns':case['all_six_columns'],'zero_columns_removed':6-len(cols),'stages':stages})
  need(counters['density_stage_checks']==7,'EXACT_SEVEN_STAGE_CHECKS')
  for p,h in packet['source_hashes'].items():need(sha(p)==h,'SOURCE_HASH_AT_COMPLETION')
  need(sha(__file__)==reg['run_sha256']and sha(P/'PROOF.md')==reg['proof_sha256']and sha(P/'next_cases.json')==reg['case_packet_sha256'],'FROZEN_PACKET_AT_COMPLETION')
  out={'status':'COMPLETE_TWO_FIXED_EXACT_DENSITY_GENERATIONS','cases':outputs,'actual_counts':counters,'run_sha256':sha(__file__),'case_packet_sha256':sha(P/'next_cases.json'),'release_sha256':sha(P/'next_parent_release.json'),'new_k_integral_evaluations':0,'new_moment_upper':False,'new_C6_bound':False,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
  write('next_result.json',out);write('next_exit_receipt.json',{'status':out['status'],'result_sha256':sha(P/'next_result.json'),'actual_counts':counters,'wall_seconds':out['wall_seconds'],'CPU_seconds':out['CPU_seconds'],'time_limit':None});print(json.dumps({'status':out['status'],'actual_counts':counters}))
 except BaseException:
  write('next_exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'actual_counts':counters,'new_k_integral_evaluations':0});raise
if __name__=='__main__':
 if sys.argv[1:]==['--saved-only-preflight']:
  packet,_=preflight();print(json.dumps({'status':'SAVED_ONLY_NEXT_PREFLIGHT_PASS','fixed_cases':len(packet['cases']),'new_density_generations':0,'new_k_integral_evaluations':0}))
 elif not sys.argv[1:]:run()
 else:raise SystemExit('usage: generate_next.py [--saved-only-preflight]')
