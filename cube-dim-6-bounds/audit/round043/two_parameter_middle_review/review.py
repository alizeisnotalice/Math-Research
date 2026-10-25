"""F043: one independent, six-configuration audit; no imported author code."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import hashlib,json,time,os,platform
import sympy as S
from flint import arb,fmpq,ctx
P=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def B(q):
 q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def ends(x):return Q(str(x.lower().fmpq())),Q(str(x.upper().fmpq()))
def pack(x):
 lo,hi=ends(x);return {'lower_exact':str(lo),'upper_exact':str(hi),'display':str(x)}
def sb(v):return Q(int(v['lower_mantissa']))*Q(2)**int(v['lower_binary_exponent']),Q(int(v['upper_mantissa']))*Q(2)**int(v['upper_binary_exponent'])
def same(x,v):
 a,b=ends(x);c,d=sb(v);assert max(a,c)<=min(b,d)

def main():
 if (P/'execution_receipt.json').exists():raise SystemExit('one audit only')
 clock=json.loads((P/'input_manifest.json').read_text());source=Path(clock['source_directory']);tic=time.monotonic()
 receipt={'status':'RUNNING','pid':os.getpid(),'started_unix':time.time(),'configurations_completed':0,'kernel_calls':0,'author_writer_calls':0}
 (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 try:
  assert time.time()-clock['started_unix']<900
  for f,h in clock['input_sha256'].items():assert sha(f)==h,f
  man=json.loads((source/'manifest.json').read_text())
  for f,h in man['whitelist_sha256'].items():assert sha(source/f)==h,f
  for f,h in man['external_readonly_sha256'].items():assert sha(f)==h,f
  saved=json.loads((source/'results.json').read_text());old={r['id']:r for r in saved['configurations']}
  rows=[]
  # 1. Independent residue derivation of the full Cauchy-resolvent integral.
  L,s,g=S.symbols('L s gamma',positive=True);a=S.symbols('alpha',real=True)
  residues_sum=1/(L**2*((S.I*s-a)**2+g*g))+s/(g*L**2*((a+S.I*g)**2+s*s))
  formula=(s+g)/(g*L**2*(a*a+(s+g)**2))
  residual=S.cancel(residues_sum-formula);assert residual==0
  rows.append({'id':1,'name':'Cauchy_average_by_two_upper_half_plane_residues','residual':str(residual),'coincident_poles':'continuous positive parameter extension proved in proof_review','b_zero':'direct resolvent bound','fiber_distance_zero_allowed':True})
  # 2. Binomial multiplicities, not the author's 32-vector sum.
  psi=Q(sum(comb(6,k)*abs(6-2*k) for k in range(7)),2)
  assert psi==60 and psi/6==10
  loss=(1/(1-Q(1,10)))**3;assert loss==Q(1000,729)
  assert old[1]['psi_ones']==psi and Q(old[1]['dimension3_peak_loss'])==loss
  rows.append({'id':2,'name':'symbol_lower_bound_by_binomial_multiplicities','psi_ones_exact':str(psi),'l1_coercivity_exact':'10','missing_symbol_fraction':'9/10','dimension3_peak_loss_exact':str(loss),'full_frame_or_fan_enumeration':False})
  # 3. One fixed geometry, analytic all-pair proof plus a single exact zero witness.
  v=S.ones(6,1);w=S.Matrix([1,1,1,-1,-1,-1]);T=S.Matrix.hstack(v,w)
  b=S.Matrix([1,1,1,1,1,-1]);c=S.Matrix([1,1,1,-1,-1,1]);R=S.Matrix.hstack(b,c)
  x=S.Matrix([1,1,1,S.Rational(1,4),S.Rational(1,4),S.Rational(1,4)])
  proj=S.eye(6)-R*(R.T*R).inv()*R.T
  assert x-v/4==S.Rational(3,8)*(b+c) and proj*(x-v/4)==S.zeros(6,1)
  assert (proj*x).dot(proj*x)==S.Rational(1,6)
  assert (T.T*T).det()==36
  JA=S.Matrix([[-1,-1],[2,-1],[-1,2]])/3
  JB=S.Matrix([[1,1,0],[1,0,1],[1,-1,-1]])
  assert (JA.T*JA).det()==S.Rational(1,3) and abs(JB.det())==3
  ctx.prec=768;pi=arb.pi();rho=Q(1,64);eps=B(Q(1,1000));eta=B(3).sqrt()*B(rho)
  hv0=B(3).sqrt()/2-eta;hv1=B(Q(3,32)).sqrt()-eta
  hw0=15*B(3).sqrt()/16-eta;hw1=B(Q(3,128)).sqrt()-eta
  gap=1/(4*B(2).sqrt())-eta;factor=1+B(Q(6,928)).sqrt()
  assert all(ends(y)[0]>0 for y in [hv0,hv1,hw0,hw1,gap])
  for n,y in [('eta',eta),('hv0',hv0),('hv1',hv1),('hw0',hw0),('hw1',hw1),('g',gap),('cauchy_average_factor',factor)]:same(y,old[3][n])
  volume=pi*pi*B(rho)**4/4
  rows.append({'id':3,'name':'fixed_middle_geometry_and_true_measure','rho':'1/64','u_interval':['1/4','3/4'],'pair_Gram_determinant':36,'coarea_J2':6,'fiber_zero_exact':True,'original_distance_squared':'1/6','A_residual_Jacobian_Gram_determinant':'1/3','B_coordinate_Jacobian_absolute':3,'actual_face_volume_formula':'pi^2*rho^4/4','volume':pack(volume),'eta':pack(eta),'hv0':pack(hv0),'hv1':pack(hv1),'hw0':pack(hw0),'hw1':pack(hw1),'other_plane_gap':pack(gap),'all_pairs_covered_by':'analytic signed-partition argument; no pair enumeration'})
  # 4. Frozen projection R2, exact coarea coefficient, alternate sixth-power budget.
  peakpath=source.parents[1]/'round033/exact_projection_peaks/results.json'
  peaks={r['id']:r for r in json.loads(peakpath.read_text())['cases']}
  r2=Q(peaks['r2_k3']['R_exact']);assert r2==Q(2051,47038464)
  C2=B(64*r2)/pi**6
  good=B(Q(2,3))*(C2*eps**9/pi**5)**B(Q(1,6))
  same(C2,old[4]['C2']);same(good,old[4]['big_branch_short_norm_upper'])
  rows.append({'id':4,'name':'two_large_parameters_short_budget','R2_exact':str(r2),'C2_times_pi6_exact':str(64*r2),'C2':pack(C2),'short_good':pack(good),'coarea_Q4_over_J2':'R2/(2*pi)^4'})
  # 5. Reconstruct restricted ordered-triple coefficient from five saved R3 values.
  ids=['r3_4110','r3_3210','r3_2220','r3_3111','r3_2211'];counts=[480,1920,480,640,1440]
  restricted=[Q(29*n,32) for n in counts];assert restricted==[435,1740,435,580,1305]
  assert sum(restricted)==comb(31,3) and 6*sum(restricted)==26970
  coefficient=24576*sum(n*Q(peaks[k]['R_exact']) for n,k in zip(restricted,ids))
  assert coefficient==Q(44743079,2970)
  crem=B(coefficient)/pi**6;common=B(loss)*crem*factor/gap**2
  mv=common/(hv0*hv1)**2;mw=common/(hw0*hw1)**2;M=mv+mw
  for n,y in [('integrated_distance_coefficient',common),('Mv',mv),('Mw',mw),('M_combined',M)]:same(y,old[5][n])
  assert old[3]['true_time_coefficients']==[32*31,32*30,32*29]
  rows.append({'id':5,'name':'all_remainder_coefficients_from_frozen_R3','restricted_unordered_counts':[int(x) for x in restricted],'ordered_count':26970,'Crem_times_pi6_exact':str(coefficient),'missing_peak_loss':str(loss),'time_coefficients':[992,960,928],'common':pack(common),'Mv':pack(mv),'Mw':pack(mw),'M':pack(M),'actual_B_le_p':True,'enlarged_sum_le_p_assumed':False})
  # 6. One actual disjoint complementary L1, integrated with alternate power order.
  bad=B(Q(6,11))*(M*eps**11/pi**5)**B(Q(1,6));total=good+bad
  same(bad,old[6]['combined_actual_complement_short_norm_upper']);same(total,old[6]['total_short_norm_upper'])
  target=Q(615635661087,125000000000000000)
  assert target==Q(man['rational_short_middle_norm_upper'])
  summed=ends(good)[1]+ends(bad)[1]
  assert ends(total)[1]<target and summed<target
  rows.append({'id':6,'name':'single_actual_complement_short_norm_and_sum','bad_short':pack(bad),'total_short':pack(total),'sum_of_independent_upper_endpoints_exact':str(summed),'strict_rational_bound':str(target),'strict_positive_margin_exact':str(target-summed),'all_saved_comparisons':'exact endpoint overlap; final independent upper test separately passed','full_time_paid':False,'whole_I_certified':False})
  assert len(rows)==6
  for f,h in clock['input_sha256'].items():assert sha(f)==h,('source drift',f)
  elapsed=time.time()-clock['started_unix'];assert elapsed<900
  out={'status':'PASS_ONE_MIDDLE_SHORT_BUDGET_ONLY','configurations':rows,'precision_bits':768,'python':platform.python_version(),'sympy':S.__version__,'pid':os.getpid(),'runtime_seconds':time.monotonic()-tic,'audit_elapsed_seconds':elapsed,'source_sha256':sha(__file__),'author_manifest_sha256':sha(source/'manifest.json'),'kernel_calls':0,'native_calls':0,'grid_calls':0,'matlab_calls':0,'old_full_enumerations':0,'author_writer_calls':0,'new_C6_bound':False}
  (P/'results.json').write_text(json.dumps(out,indent=2)+'\n')
  receipt.update(status='PASS',configurations_completed=6,elapsed_seconds=time.monotonic()-tic,finished_unix=time.time(),audit_elapsed_seconds=elapsed)
  print(json.dumps({'status':out['status'],'pid':os.getpid(),'runtime_seconds':receipt['elapsed_seconds'],'audit_elapsed_seconds':elapsed,'good':str(good),'bad':str(bad),'total':str(total),'strict_target':str(target)}))
 except BaseException as exc:
  receipt.update(status='FAILED_STOP_NO_RETRY',error=repr(exc),finished_unix=time.time(),elapsed_seconds=time.monotonic()-tic)
  raise
 finally:
  (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':main()
