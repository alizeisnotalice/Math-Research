"""Exact rational ledger for sufficient inequalities in the explicit cubic proof.
The ledger covers scalar bounds on the entire t interval; it is not a search
certificate for arbitrary measurable inputs. That implication is proved in TeX.
"""
from fractions import Fraction as F
from pathlib import Path
import json,platform,time,sys
ROOT=Path(__file__).resolve().parents[1]/'results'

def polynomial_range(coefs,tau):
 lo=hi=coefs[0]
 for k,c in enumerate(coefs[1:],1):
  lo+=min(c,0)*tau**k;hi+=max(c,0)*tau**k
 return lo,hi

def verify(e_eta,e_tau,e_kappa):
 st=time.monotonic();tau=F(1,2**e_tau);eta=F(1,2**e_eta);kappa=F(1,2**e_kappa);sqrt_k=F(1,2**(e_kappa//2));sqrt_tau=F(1,2**(e_tau//2));s0=F(1,2**(e_eta//2+4));a=F(11,20);b=F(9,20);d=a-b
 rows=[]
 def positive(name,expr,reason):
  assert expr>0,(name,expr)
  rows.append(dict(name=name,positive_margin=str(expr),reason=reason))
 positive('prior_curvature_range',F(1,100)-tau,'prior cap and curvature lemmas require t<=.01')
 positive('chart_sampling_range',F(1,100)-s0*F(17,16),'s0 +/- h stay in the chart interval')
 positive('pair_coordinate_sign',2*b-4*eta,'coordinate signs of opposite pair are stable')
 # F(t+c*t^2)/t^2; derivative F_r<0 on 0<r<2t.
 poly=[]
 for name,c in [('lower',F(9,10)),('upper',F(2))]:
  coefs=[1+4*a*b-2*c,2*c+8*a*b*c,4*a*b*c*c];lo,hi=polynomial_range(coefs,tau)
  assert lo>0 if name=='lower' else hi<0
  poly.append(dict(name=name,c=str(c),coefficients=list(map(str,coefs)),interval=list(map(str,[lo,hi]))))
 positive('root_monotonicity',1-2*tau-16*a*b*tau,'F_r <= -2+2tau+16ab*tau < -1')
 positive('normal_length_squared',2*(1-3*tau)**2-1,'each unnormalized normal coordinate >=1-3tau')
 positive('tangent_coordinate_sum',F(3,10)-2*d*(1+2*tau),'N>=1 gives |v1+v2|/t <=2d(1+2tau)')
 positive('pair_sum_perturbation',4-2*(F(3,2)+F(3,10)),'sqrt2<3/2, two points')
 positive('pair_product_perturbation',10-(8*(1+2*tau)+16*eta),'coordinate perturbations <=4eta*t')
 positive('endpoint_order',F(18,10)-1/(1-tau)-4*eta-F(3,4),'S-W >= (9/5-1/alpha-4eta)t^2 >3t^2/4')
 positive('pair_budget',16-(4+10),'G(alpha)<=14eta*t^2 <=16eta*t^2')
 positive('edge_uncertainty',20-(2+16),'point coordinate 2eta*t plus individual half-scale error16eta*t')
 positive('mixed_low_side',d-4*b*tau-40*eta-F(9,100),'short side <=1-.09t')
 positive('mixed_low_mass',F(9,100)-32*eta*tau,'lost mass budget <=32eta*t^2')
 high_gain=2*d-2*(3*tau+40*eta)-2*d*tau*(3*tau+40*eta)-3*d*d*tau
 positive('mixed_high_union',high_gain-F(1,10),'cross-rectangle union >=1+.1t')
 positive('mixed_high_mass',F(1,10)-32*eta*tau,'missing mass <=32eta*t^2')
 positive('corner_width',1+d-40*eta-1,'two bad corners width >=t')
 positive('corner_height',1-d*(1+2*tau)-40*eta-F(4,5),'height >=.8t')
 positive('corner_contradiction',F(8,5)-32*eta,'corner area >=1.6t^2 > pair budgets')
 angle=d*(1+2*tau)/(1-2*tau-2*tau*tau)
 positive('normal_window',F(11,100)-angle,'arctan(x)<=x gives offset <=.11t')
 positive('chart_window',F(1,8)-F(11,100)-F(1,100),'homogeneous chart |s|<=.01 stays in |theta-theta_k|<t/8')
 positive('chart_curvature_lower',F(9,4)-(1+F(1,10000)*tau*tau)**3,'rho>=1/2 and denominator <=3/2 gives f_T second derivative >=1/3')
 positive('local_cap_condition',F(1,2**22)-kappa,'delta<=kappa*t^3 implies prior local cap hypothesis')
 E=32*F(1,2**((2*e_kappa)//3))+2*sqrt_k
 positive('support_chart_error',s0*s0/4096-E,'kappa=2^-m implies kappa^(2/3)<=2^-floor(2m/3); exact rational majorant')
 positive('positive_subgradient',F(1,3)-F(1,32)-F(1,128)-F(1,4),'h=s0/16, E<=s0^2/4096')
 positive('bounded_subgradient',2-(1+F(1,32)+F(1,128)),'upper derivative <=2s0')
 positive('positive_chart_value',F(1,6)-F(1,4096),'f_A>=s0^2/6-E>0')
 positive('exposed_normal',eta/64-3*s0*s0,'|f_A-s*p|<=3s0^2<eta/64')
 positive('origin_interior',F(1,4)-2*sqrt_k*sqrt_tau,'h_A/t>=1/2-2sqrt(kappa*t)>1/4')
 positive('triangle_height',F(1,2)-eta*tau/64-F(1,3),'core to chord normal separation >=t/3')
 positive('triangle_core_weight',1-(F(1,64)+F(5,8)),'weight <=3eta*t')
 positive('triangle_width',(1-3*eta*tau)*s0/4-6*eta*tau-eta/8,'all target depths have left/right endpoints outside +/-eta*t/8')
 positive('inner_box_geometry',F(3,4)*eta-F(25,32)*eta*eta*tau*tau-eta*eta/32,'F(p+z)/t^2 >=2Nu/t^2-2|z|^2/t^2')
 positive('inner_box_positive_coordinates',b-eta,'p coordinates >=b*t, displacement <=eta*t')
 positive('inner_box_branch',F(1,2)*(1-tau)-(a*(1+2*tau)+eta)*tau,'coordinates stay below alpha/2')
 positive('hole_area_budget',eta*eta/16-kappa,'one missed P has area eta^2*t^3/16 >delta')
 out=dict(python=platform.python_version(),method='exact Fraction monotone interval bounds',eta=str(eta),t_max=str(tau),kappa=str(kappa),root_polynomials=poly,inequalities=rows,count=len(rows),seconds=time.monotonic()-st,scope='all scalar sufficient conditions on 0<t<=t_max; analytic geometry and mass implications are proved separately')
 return out
if __name__=='__main__':
 profiles=[verify(*powers) for powers in [(14,8,72),(12,14,68)]]
 if '--ledger-only' in sys.argv:
  print({'exact_scalar_inequalities':sum(p['count'] for p in profiles),'profiles':len(profiles),'dependencies':'Python standard library only'});raise SystemExit
 import mpmath as mp  # pinned in requirements.txt; optional for --ledger-only
 mp.mp.dps=80;k1=mp.mpf(2)**-72;k2=mp.mpf(2)**-68;a1=1-mp.mpf(2)**-8;a2=1-mp.mpf(2)**-14
 def exact_tail(p,a):return sum((-1)**j*[1,3,3,1][j]*p*a**(p+j)/(p+j) for j in range(4))
 moments=[]
 for pp in [3,10,100,1000,10000,100000,1000000]:
  p=mp.mpf(pp);full=6*k2/((p+1)*(p+2)*(p+3));tail=k1*exact_tail(p,a1)+(k2-k1)*exact_tail(p,a2);D=full-tail
  independent=k1*p*mp.betainc(4,p,0,1-a1)+(k2-k1)*p*mp.betainc(4,p,0,1-a2)
  assert D>0 and abs(D-independent)<mp.mpf('1e-65')*full
  assert 0<=tail<=k1*a1**p+(k2-k1)*a2**p
  moments.append(dict(p=pp,de=mp.nstr(D,55),normalized=mp.nstr(p**3*D,55),limit=mp.nstr(6*k2,55),independent_beta_error=mp.nstr(abs(D-independent),10)))
 (ROOT/'explicit_cubic_python.json').write_text(json.dumps({'profiles':profiles,'mpmath':mp.__version__,'moment_precision':80,'moment_regressions':moments},indent=2)+'\n')
 print([dict(count=p['count'],eta=p['eta'],t_max=p['t_max'],kappa=p['kappa']) for p in profiles])
