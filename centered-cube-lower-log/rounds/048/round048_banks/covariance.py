from pathlib import Path
from fractions import Fraction as F
import json,math,time
import numpy as np
from ranks import frac,mul,power
HERE=Path(__file__).resolve().parent

def ramp(x,w):
 x=frac(x)
 return min(F(1),max(F(0),(x-F(1,2)+w/2)/w))-x

def J(s,w):
 """Exact piecewise integration of a periodic ramp difference squared."""
 s=frac(s);knots=[F(0),F(1),F(1,2)-w/2,F(1,2)+w/2]
 cuts=sorted(set(knots+[frac(x-s) for x in knots]));total=F(0)
 for left,right in zip(cuts,cuts[1:]):
  a=ramp(left+s,w)-ramp(left,w);b=ramp(right+s,w)-ramp(right,w)
  total+=(right-left)*(a*a+a*b+b*b)/3
 return total,len(cuts)-1

def enc(x):return str(x)

def outward(x,places=9,up=False):
 scale=10**places;t=x*scale;q=t.numerator//t.denominator
 if up and t.denominator!=1:q+=1
 sign='-'if q<0 else '';q=abs(q)
 return f'{sign}{q//scale}.{q%scale:0{places}d}'

def source_value(labels,j,theta,off,eps,B,K):
 def end(t):
  k=t.numerator//t.denominator;r=t-k;u=F(2*labels[(j+k)%K]+1,2*B)
  return min(F(1),max(F(0),(r-u+eps/2)/eps))-r
 return end(theta+off)-end(theta-off)

def small_direct_check():
 d=8;p=285;K=255;q=256;B=4;bases=[[1,6],[254,134]];L=F(10,9);eps=F(1,128);off=K*L/2
 states=[];x=1
 for _ in range(K):states.append(x);x=mul(x,2,d,p)
 def tr(x):
  z=0
  for _ in range(d):z^=x;x=mul(x,x,d,p)
  assert z in [0,1]
  return z
 labels=[[sum(tr(mul(beta,x,d,p))<<i for i,beta in enumerate(basis))for x in states]for basis in bases]
 cuts=sorted({F(0),F(1),frac(off),frac(-off)});cross=F(0);diag=F(0)
 for j in range(K):
  for left,right in zip(cuts,cuts[1:]):
   mid=(left+right)/2;local={left,right}
   for labs in labels:
    for sign in [1,-1]:
     t=sign*off;k=(mid+t).numerator//(mid+t).denominator;u=F(2*labs[(j+k)%K]+1,2*B)
     for e in [-1,1]:
      z=k+u+e*eps/2-t
      if left<z<right:local.add(z)
   points=sorted(local)
   for a,z in zip(points,points[1:]):
    fa=source_value(labels[0],j,a,off,eps,B,K);fz=source_value(labels[0],j,z,off,eps,B,K);ga=source_value(labels[1],j,a,off,eps,B,K);gz=source_value(labels[1],j,z,off,eps,B,K)
    cross+=(z-a)*(2*fa*ga+fa*gz+fz*ga+2*fz*gz)/6
    diag+=(z-a)*(fa*fa+fa*fz+fz*fz)/3
 cross/=K;diag/=K;jmean,_=J(B*K*L,B*eps);jmean/=B*B;jzero,_=J(K*L,eps);var=F(B*B-1,6*B*B)-eps*F(B-1,6*B)
 predicted_cross=F(q,K)*jmean-jzero/K;predicted_diag=F(q,K)*(2*var+jmean)-jzero/K
 assert cross==predicted_cross and diag==predicted_diag
 return {'field':{'d':d,'modulus':p,'generator':2},'bases':bases,'L':str(L),'epsilon':str(eps),'exact_direct_H_cross':str(cross),'exact_rank_formula_H_cross':str(predicted_cross),'exact_direct_H_diagonal':str(diag),'exact_rank_formula_H_diagonal':str(predicted_diag),'passed':True,'scope':'Exact full small-field state sum and piecewise-theta integration; small test bases need no Sidon for this identity.'}

def main():
 start=time.perf_counter();m=json.loads((HERE/'MANIFEST.json').read_text());plan=json.loads((HERE/'FROZEN_PLAN.json').read_text());ranks=json.loads((HERE/'RANKS.json').read_text());K=m['K'];q=m['q'];B=m['B'];b=m['b'];p=F(plan['p']);n=int(plan['n']);w=plan['w'];M=plan['M'];v=F(9)*p*p*n/(29*K*K);eps=1/(w*v);delta=eps/K
 # First: independent exact small finite orbit and piecewise integration.
 small=small_direct_check();(HERE/'SMALL_COVARIANCE_CHECK.json').write_text(json.dumps(small,indent=2)+'\n')
 var=F(B*B-1,6*B*B)-eps*F(B-1,6*B)
 harmonic=sum(F(1,j)for j in range(1,M+1));lam_up=F(math.isqrt(q)+1,K);U=2*(b+1);Vglobal=F(40,9);tail=F(48,29*M**3);gamma=F(27,29)
 mm=np.arange(1,M+1);alpha=F(729,29);bm=4*float(alpha)*((np.sin(np.pi*mm*10/27)/(np.pi*mm))**4+3*(np.sin(np.pi*mm/9)/(np.pi*mm))**4);cm=np.sqrt(2*bm/3);rows=[]
 for ls in plan['scales']:
  L=F(ls);cases=[r for r in ranks['rank_cases']if r['L']==str(L)];assert cases
  assert all(x==22 for r in cases for t in r['pair_banks']for row in t['cross_endpoint_pair_ranks']for x in row)
  near=L-round(L);sin=np.array([((-1)**(int(round(L))*int(j)))*math.sin(math.pi*int(j)*float(near))for j in mm]);BB=float(np.sum(bm*sin*sin)/3);VL=float(np.sum(cm*np.abs(sin)))
  if L.denominator==1:
   rows.append({'L':str(L),'integer_scale':True,'exact_H_diagonal':'0','exact_H_cross':'0','exact_BB':'0','correlation':None,'reason':'All Y_b=1; variance0 so Pearson correlation undefined.'});continue
  assert all(r['same_bank_two_endpoint_ranks']==[22]*4 for r in cases)
  jmean,stepsmean=J(B*K*L,B*eps);jmean/=B*B;jzero,stepszero=J(K*L,eps)
  HHcross=F(q,K)*jmean-jzero/K;HHdiag=F(q,K)*(2*var+jmean)-jzero/K
  dist=abs(near);Vbound=min(Vglobal,F(440,63)*dist*harmonic);HBbound=lam_up*U*Vbound
  if F(10,9)<=L<=F(11,9):
   plateau=gamma*L*L-1;BBlo=max(F(0),(plateau-tail)/3);BBhi=plateau/3;BBscope='Accepted F_infty=1+sum b_m sin²=gamma L²; F_M includes the constant1; 0<=F_infty-F_M<=48/(29*M³).'
  else:BBlo=F(0);BBhi=Vbound*Vbound/2;BBscope='Orthogonal Fourier modes: BB=1/2 sum c_m²sin² <= V_L²/2; sin bound pi*m*dist(L,Z), pi<22/7.'
  assert float(BBlo)-1e-15<=BB<=float(BBhi)+1e-15
  Dlo=HHdiag+BBlo-2*HBbound;Dhi=HHdiag+BBhi+2*HBbound;Clo=HHcross+BBlo-2*HBbound;Chi=HHcross+BBhi+2*HBbound;assert Dlo>0
  corrlo=max(F(-1),Clo/(Dhi if Clo>=0 else Dlo));corrhi=min(F(1),Chi/(Dlo if Chi>=0 else Dhi))
  center=(float(HHcross)+BB)/(float(HHdiag)+BB);share=BB/(float(HHdiag)+BB);na2=float(F(n)*p*p/(K*K*L*L))
  rows.append({'L':str(L),'integer_scale':False,'exact_H_diagonal':enc(HHdiag),'exact_H_cross_all_six_pairs':enc(HHcross),'exact_mean_H_squared_integral':enc(jmean),'exact_zero_state_H_squared_integral':enc(jzero),'theta_segment_counts':{'compressed_uniform_mean_ramp':stepsmean,'zero_state_ramp':stepszero},'BB_numerical':BB,'BB_rational_bounds':[enc(BBlo),enc(BBhi)],'BB_bound_basis':BBscope,'HB_abs_rational_upper_bound_each_bank':enc(HBbound),'HB_bound_decimal_reference':float(HBbound),'V_L_numerical':VL,'V_L_rational_upper':enc(Vbound),'total_variance_rational_bounds':[enc(Dlo),enc(Dhi)],'cross_bank_covariance_rational_bounds':[enc(Clo),enc(Chi)],'Pearson_correlation_rational_bounds':[enc(corrlo),enc(corrhi)],'Pearson_correlation_outward_decimal_9_places':[outward(corrlo),outward(corrhi,up=True)],'baseline_BB_variance_fraction_numeric':share,'baseline_correlation_center_numeric':center,'n_scaled_linear_d_variance_baseline_numeric':na2*(float(HHdiag)+BB),'n_scaled_shared_B_variance_numeric':na2*BB,'note':'Strict rational bounds use accepted spectral identity/tail and rational Gauss-Walsh envelope. Numerical BB/centers are references, not interval endpoints. This is covariance of H+B or linear d, not log-score covariance.'})
 out={'parameters':{'q':q,'K':K,'B':B,'p':str(p),'n':str(n),'w':w,'M':M,'v':str(v),'v_numeric':float(v),'epsilon_Kdelta':str(eps),'delta':str(delta),'epsilon_less_than_1_over_B':eps<F(1,B),'lambda_sqrtq_over_K_rational_upper':str(lam_up),'source_density_lower_bound':str(1-p-F(20,3)*p*harmonic/K)},'small_exact_test':small,'scales':rows,'seconds':time.perf_counter()-start,'claim_scope':'Exact HH second moments and analytic/rational HB,BB,corr enclosures for these measured-rank scales. No Fourier independence, no Gaussian transfer, no n-dimensional Monte Carlo.'}
 assert F(out['parameters']['source_density_lower_bound'])>0
 (HERE/'COVARIANCE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'small_exact_passed':small['passed'],'seconds':out['seconds'],'rows':[{k:r[k] for k in ['L','BB_numerical','baseline_BB_variance_fraction_numeric','Pearson_correlation_outward_decimal_9_places']}if not r['integer_scale']else{'L':r['L'],'correlation':None}for r in rows]},indent=2))
if __name__=='__main__':main()
