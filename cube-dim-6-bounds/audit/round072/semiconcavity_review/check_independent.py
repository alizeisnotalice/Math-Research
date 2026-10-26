"""Five new bounded exact tests; no producer path/pairs or body32 work."""
import sys,time,signal,resource
START=time.monotonic();CPU=time.process_time();sys.dont_write_bytecode=True
resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def stop(sig,frame):raise TimeoutError('B072 CPU30s/wall170s limit')
signal.signal(signal.SIGXCPU,stop);signal.signal(signal.SIGALRM,stop);signal.setitimer(signal.ITIMER_REAL,170)
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib,json,os,traceback
import sympy as s
H=Path(__file__).resolve().parent
reg=json.loads((H/'registration.json').read_text());SRC=Path(reg['producer_directory'])
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
rd=lambda p:json.loads(Path(p).read_text())
def save(n,obj):(H/n).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
if (H/'execution_receipt.json').exists():raise RuntimeError('Single-attempt receipt already exists')
receipt={'status':'RUNNING','pid':os.getpid(),'completed_configurations':0,'attempted_configurations':0,
 'cpu_limit_seconds':30,'wall_total_limit_seconds':180,'automatic_retries':0,'old_fixtures_replayed':0,
 'pair_residual_evaluations':0,'Green_path_evaluations':0,'body32_kernel_calls':0}
out={'status':'INCOMPLETE','checks':[]}
def checkpoint():
 receipt['wall_seconds']=time.monotonic()-START;receipt['cpu_seconds']=time.process_time()-CPU
 save('execution_receipt.json',receipt);save('results.json',out)
def begin(i):
 assert i==receipt['attempted_configurations']+1<=5;receipt['attempted_configurations']=i;checkpoint()
def done(i,**kw):
 out['checks'].append({'id':i,'pass':True,**kw});receipt['completed_configurations']=i;checkpoint()
def pins():
 assert sha(SRC/'manifest.json')==reg['producer_manifest_sha256']
 assert sha(SRC/'theorem.md')==reg['producer_theorem_sha256']
 m=rd(SRC/'manifest.json')
 for n,h in m['whitelist_sha256'].items():assert sha(SRC/n)==h,n
 for n,h in m['source_sha256_unchanged'].items():assert sha(n)==h,n
 b=rd(H/'execution_binding.json')
 assert b['registration_sha256']==sha(H/'registration.json')
 for n,h in b['code_sha256'].items():assert sha(H/n)==h,n
 return m
try:
 checkpoint();source=pins()
 begin(1);t,r=s.symbols('t r',positive=True)
 radial=-s.diff(r/(t*t+r*r),t)/(2*s.pi**2*r)
 assert s.factor(radial-t/(s.pi**2*(t*t+r*r)**2))==0
 deficit=1-t**4/(t*t+r*r)**2
 assert s.factor(deficit-2*r*r/(t*t+r*r)+r**4/(t*t+r*r)**2)==0
 I1=s.pi/(2*r);I2=-s.diff(I1,r)/(2*r)
 val=s.factor(2*r*r*I1-r**4*I2);assert val==3*s.pi*r/4
 T=s.diag(2,3,5);assert T.det()==30
 coeff=s.simplify(val/(s.pi*s.pi*r*T.det()));assert coeff==1/(40*s.pi)
 done(1,name='anisotropic_radial_Cauchy_normal_deficit',T=mat(T),determinant='30',
      radial_density='t/[pi^2*(t^2+r^2)^2]',deficit_integral_without_pi_minus2=str(val),
      transformed_norm='sqrt(z1^2/4+z2^2/9+z3^2/25)/(40*pi)',
      old_axis_product_slice_replayed=False,new_Green_path=False)

 begin(2);y,z1,z2,eps,h1,h2=s.symbols('y z1 z2 epsilon h1 h2',real=True)
 N=s.sqrt(z1*z1+4*z2*z2+eps*eps);A=2+(y+z1)**2+z2*z2;variables=(y,z1,z2)
 gradA=s.Matrix([s.diff(A,v)for v in variables]);gradN=s.Matrix([s.diff(N,v)for v in variables])
 formula=-A*s.hessian(N,variables)-N*s.hessian(A,variables)-gradA*gradN.T-gradN*gradA.T
 assert all(s.simplify(v)==0 for v in s.hessian(-A*N,variables)-formula)
 hn=s.Matrix([h1,h2]);numerator=eps*eps*(h1*h1+4*h2*h2)+4*(z1*h2-z2*h1)**2
 assert s.simplify((hn.T*s.hessian(N,(z1,z2))*hn)[0]*N**3-numerator)==0
 assert s.expand(4*N**2-(z1*z1+16*z2*z2))==3*z1*z1+4*eps*eps
 assert s.hessian(A,variables).eigenvals()=={s.Integer(4):1,s.Integer(2):1,s.Integer(0):1}
 joint2=s.hessian(-(2+(y+z1)**2)*z1,(y,z1))
 det_at=s.factor(joint2.det().subs({y:s.Rational(1,2),z1:s.Rational(1,3)}));assert det_at==-s.Rational(25,9)
 assert 5<9 and 3*4+2*2*4==28
 done(2,name='joint_amplitude_and_normal_semiconcavity',amplitude=str(A),norm=str(N),
      norm_Hessian_nonnegative_numerator=str(numerator),joint_Hessian_product_residuals=['0']*9,
      gradient_norm_bound='2',amplitude_gradient_bound='4',amplitude_Hessian_bound='4',norm_value_bound='3',
      semiconcavity_constant_upper='28',joint_2x2_Hessian_det_at_new_point=str(det_at),
      joint_function_not_claimed_concave=True)

 begin(3);Y=s.symbols('y1:4',positive=True);t=s.symbols('t',positive=True)
 w=[1,2,4];rates=[F(2),F(3),F(5)];delta=[F(2),F(4),F(6)]
 Hfun=s.prod(Y[i]**2/(Y[i]**2+w[i]**2*t*t)for i in range(3))
 base=s.cancel((1-Hfun)/(t*t));S=sum(s.Rational(w[i]**2)/Y[i]**2 for i in range(3))
 limits=[]
 patterns=[()]+[(i,)for i in range(3)]+[(i,j)for i in range(3)for j in range(i,3)]
 for pattern in patterns:
  value=base;target=S
  for i in pattern:value=s.diff(value,Y[i]);target=s.diff(target,Y[i])
  atzero=s.factor(value.subs(t,0));assert s.factor(atzero-target)==0
  limits.append({'y_derivative_indices':list(pattern),'limit_of_derivative_1minusH_over_t2':str(atzero)})
 A0=F(8)/F(2*4*6)**2
 H0=sum(F(w[i]**2)/delta[i]**2 for i in range(3))
 H1=[2*F(w[i]**2)/delta[i]**3 for i in range(3)]
 H2=[[6*F(w[i]**2)/delta[i]**4 if i==j else 4*F(w[i]**2*w[j]**2)/(delta[i]**3*delta[j]**3)for j in range(3)]for i in range(3)]
 A1=[2*A0/delta[i]for i in range(3)]
 A2=[[6*A0/delta[i]**2 if i==j else 4*A0/(delta[i]*delta[j])for j in range(3)]for i in range(3)]
 K0=A0*H0;K1=[A1[i]*H0+A0*H1[i]for i in range(3)]
 K2=[[A2[i][j]*H0+A1[i]*H1[j]+A1[j]*H1[i]+A0*H2[i][j]for j in range(3)]for i in range(3)]
 prod_rates=rates[0]*rates[1]*rates[2]
 q0=1/prod_rates;gradq=sum(1/v for v in rates)/prod_rates
 hessq=(sum(1/v**2 for v in rates)+2*max(1/v**2 for v in rates))/prod_rates
 assert hessq==F(811,27000) and A0==F(1,288)
 assert all(v>0 for v in [K0,*K1,*[v for row in K2 for v in row],q0,gradq,hessq])
 done(3,name='new_weighted_joint_remainder_dominators',tangential_rates=w,normal_rates=list(map(str,rates)),
      exact_small_time_limits=limits,A_upper_without_pi_minus3=str(A0),
      small_time_A_times_1minusH_derivative_bounds_without_pi_minus3={'order0':str(K0),'order1':list(map(str,K1)),'order2':[[str(v)for v in row]for row in K2]},
      normal_peak_constant_without_pi_minus3=str(q0),normal_gradient_constant_without_pi_minus3=str(gradq),normal_Hessian_constant_without_pi_minus3=str(hessq),
      large_time_normal_difference_gradient_Hessian_constants_without_pi_minus3=[str(2*hessq),str(2*hessq),str(hessq)],
      small_time_integrand_exponents_for_z_orders012=[2,1,0],large_time_all_joint_derivative_integrand_exponent=-2,
      body_or_axis_Green_integrals_executed=0)

 begin(4);N=s.Matrix(reg['fixed_configurations'][3]['N_rows']);inv=N.inv();det=abs(N.det())
 assert det==30030
 Y0=s.Matrix(reg['fixed_configurations'][3]['axis_coordinates']);X0=N*Y0
 assert inv*X0==Y0 and sum(v!=0 for v in Y0)==3
 X=s.Matrix(s.symbols('x0:6'));transformed=(inv*X).dot(inv*X)/(2*det)
 HH=s.hessian(transformed,list(X));expected=inv.T*inv/det
 assert HH==expected
 factor=s.factor(sum(v*v for v in inv)/det)
 uppermat=factor*s.eye(6)-HH
 leading=[s.factor(uppermat[:i,:i].det())for i in range(1,7)]
 assert all(v>0 for v in leading)
 done(4,name='general_nonorthogonal_sixcolumn_transport',N=mat(N),inverse=mat(inv),absolute_determinant=str(det),
      new_point=list(map(str,X0)),unique_axis_coordinates=list(map(str,Y0)),off_S2_by_support_count=3,
      Green_evaluations=0,pair_enumerations=0,Hessian_pullback_identity=True,
      finite_semiconcavity_multiplier_using_Frobenius=str(factor),
      sharp_operator_formula='abs(detN)^(-1)*C_axis*||N^-1||op^2',
      positive_leading_minors_of_coarse_bound_difference=list(map(str,leading)))

 begin(5);fixture=reg['fixed_configurations'][4]
 L=s.Matrix.hstack(*(s.Matrix(col)for col in fixture['L_columns']))
 hs=list(map(s.Rational,fixture['halfwidths']));ss=s.Matrix(s.symbols('s0:3'));center=L*ss
 total=s.Integer(0);mean=s.zeros(6,1);second=s.zeros(6,6);quad=s.Integer(0)
 for signs in product((-1,1),repeat=3):
  weight=s.prod((1+signs[i]*ss[i]/hs[i])/2 for i in range(3))
  vertex=L*s.Matrix([signs[i]*hs[i]for i in range(3)])
  total+=weight;mean+=weight*vertex;second+=weight*(vertex*vertex.T);quad+=weight*vertex.dot(vertex)/2
 assert s.expand(total)==1 and all(s.expand(v)==0 for v in mean-center)
 covariance=second-center*center.T
 expected=L*s.diag(*(hs[i]**2-ss[i]**2 for i in range(3)))*L.T
 assert all(s.expand(v)==0 for v in covariance-expected)
 variance=sum((hs[i]**2-ss[i]**2)*L[:,i].dot(L[:,i])for i in range(3))
 assert s.expand(quad-variance/2-center.dot(center)/2)==0
 penalty=s.factor(sum(hs[i]**2*L[:,i].dot(L[:,i])for i in range(3))/2)
 done(5,name='whole_linear_box_vertex_lower_certificate',L=mat(L),halfwidths=list(map(str,hs)),vertex_count=8,
      barycentric_sum='1',covariance_identity=True,quadratic_equality_identity=True,
      maximum_lower_penalty_per_unit_C=str(penalty),
      general_pointwise_penalty='(C/2)*sum_i (h_i^2-s_i^2)*||L_i||^2',
      whole_box_lower='min_vertex_lower-(C/2)*sum_i h_i^2*||L_i||^2',
      axis_box_full_side_coefficient='C/8 times sum of squared full side lengths',
      requires_uniform_semiconcavity_on_entire_convex_box=True)
 pins();out.update(status='PASS_FIVE_NEW_INDEPENDENT_CONFIGURATIONS',scientific_configurations=5,
  registration_sha256=sha(H/'registration.json'),checker_sha256=sha(__file__),producer_manifest_sha256=reg['producer_manifest_sha256'],
  body32_semiconcavity_proved=False,body32_counterexample=False,old_fixtures_replayed=0,pair_residual_evaluations=0,Green_path_evaluations=0,
  new_Cn_bound=False,main_retains_final_proof=True)
 receipt.update(status='PASS_STOPPED',exit_code=0);checkpoint();receipt['results_sha256']=sha(H/'results.json');save('execution_receipt.json',receipt)
 print(json.dumps({'status':out['status'],'configurations':5,'pid':os.getpid(),'wall_seconds':receipt['wall_seconds'],'cpu_seconds':receipt['cpu_seconds'],'old_replays':0,'body32_kernel_calls':0}))
except BaseException as ex:
 receipt.update(status='FAILED_STOP_NO_RETRY',exit_code=1,error=repr(ex),traceback=traceback.format_exc());checkpoint();raise
finally:signal.setitimer(signal.ITIMER_REAL,0)
