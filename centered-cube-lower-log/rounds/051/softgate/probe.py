"""Discrete positive soft sources, actual box averages, stable posterior diagnostics.
No Gaussian gate is substituted for the source. Uses only installed NumPy.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import math,json,hashlib,time
import numpy as np
HERE=Path(__file__).resolve().parent

def metrics(Y,c,pi):
 d=Y-1;weights=Y/Y.sum();shift=c*d/Y;logY=np.log1p(d);ms=float(weights@shift);ml=float(weights@logY);ds=shift-ms;dl=logY-ml;vp=float(weights@(ds*ds));vl=float(weights@(dl*dl));cov=float(weights@(ds*dl))
 return {'r_maximal_binary_query_channel':math.sqrt(vp/(pi*(1-pi))),'Corr_p_logY_under_Y_anchor':cov/math.sqrt(vp*vl) if vp*vl else None,'posterior_std':math.sqrt(vp),'query_log_std':math.sqrt(vl),'posterior_shift_abs_max':float(np.max(np.abs(shift))),'mean_Y':float(Y.mean())}

def main():
 start=time.perf_counter();rng=np.random.default_rng(20261029051);B=16;rows=[]
 for K in [31,127,511]:
  Q=K*B;positions=rng.integers(0,B,size=K);peaks=np.zeros(Q,dtype=np.int64);peaks[np.arange(K)*B+positions]=1;prefix=np.r_[0,np.cumsum(np.tile(peaks,2))];idx=np.arange(Q)
  for minrho in [.1,.01,1e-8,1e-14]:
   rho=minrho+(1-minrho)*B*peaks;m=float(rho.min());M=float(rho.max())
   for pi in [.5,.01]:
    cm=-min((1-pi)*m/(1-m),pi*M/(M-1));cp=min(pi*m/(1-m),(1-pi)*M/(M-1));c=.95*cp;rho1=pi*rho+c*(rho-1);rho0=rho-rho1;assert rho1.min()>=-1e-25 and rho0.min()>=-1e-25
    rinput=math.sqrt(float(np.mean((c*(rho-1))**2/rho))/(pi*(1-pi)))
    for L in [F(10,9),F(7,6),F(11,9)]:
     S=(L*Q).numerator//(L*Q).denominator
     # Actual cyclic box of S bins, including its complete-period contribution.
     full,rem=divmod(S,Q);count=full*K+prefix[idx+rem]-prefix[idx];d=(1-minrho)*(B*count-S)/S;Y=1+d
     met=metrics(Y,c,pi)
     # Exponential latent-count tilt: exact affine product family at any t.
     tilt=[]
     for t in [-3.,0.,3.]:
      e=math.exp(t);Z=1-pi+pi*e;alpha=1+c*(e-1)/Z;rt=(rho0+e*rho1)/Z;expected=1+alpha*(rho-1);assert np.max(np.abs(rt-expected))<1e-12
      tilt.append({'t':t,'alpha':alpha,'density_mean':float(rt.mean()),'source_identity_max_abs_error':float(np.max(np.abs(rt-expected)))})
     rows.append({'K_coarse_cells':K,'B_bins_per_cell':B,'Q_query_nodes':Q,'source_minrho':m,'source_maxrho':M,'pi':pi,'legal_c_interval':[cm,cp],'chosen_c':c,'input_binary_channel_r':rinput,'L_target':str(L),'actual_box_bins':S,'actual_L':str(F(S,Q)),'max_abs_Y_minus1':float(np.max(np.abs(d))),**met,'exponential_tilt_identity_checks':tilt,'classification':'Actual finite discrete positive source/box averages; numerical evidence, not a continuous HL transfer.'})
 # Exact rational real-source envelope from the accepted positive peak/Fourier construction.
 K=2**47-1;q=K+1;n=10**32;p=F(99,100);w=64;M=40;L=F(7,6);pi=F(1,2);harm=sum(F(1,j)for j in range(1,M+1));v=9*p*p*n/(29*K*K);low=F(20,3)*p*harm/K
 mlower=1-p-low;mupper=1-p+low;Mlower=1-p+p*w*v-low;Mupper=1-p+p*w*v+low
 assert 0<mlower<=mupper<1<Mlower<=Mupper
 cplus=min(pi*mlower/(1-mlower),(1-pi)*Mupper/(Mupper-1));c=F(19,20)*cplus;h=p/(K*L)*F(58,9);assert h<1
 r2=c*c*h*h/(pi*(1-pi)*(1-h));a_up=mupper/(1-mupper);b_up=Mlower/(Mlower-1);r2_all_pi=a_up*b_up*h*h/(1-h)
 with localcontext()as ctx:
  ctx.prec=75
  def decimal(x):return D(x.numerator)/D(x.denominator)
  real={'n':str(n),'K':K,'p':str(p),'L':str(L),'v':str(v),'source_min_lower':str(mlower),'source_min_upper_from_peak_gap':str(mupper),'source_max_lower_at_peak':str(Mlower),'source_max_upper':str(Mupper),'safe_chosen_c_pi_half':str(c),'h_query_uniform_upper':str(h),'r2_upper_chosen_c':str(r2),'r_upper_chosen_c_decimal':str(decimal(r2).sqrt()),'r2_upper_every_legal_pi_c':str(r2_all_pi),'r_upper_every_legal_pi_c_decimal':str(decimal(r2_all_pi).sqrt()),'Gaussian_reference_u':str(decimal(v).sqrt()),'Gaussian_reference_shift2_upper_every_legal_pi_c':str(r2_all_pi*(v+2+1/v)),'Gaussian_reference_shift_upper_decimal':str(decimal(r2_all_pi*(v+2+1/v)).sqrt()),'classification':'Exact rational squared envelopes; Decimal square roots are display only. Gaussian shift formula is a calibration conditional on a Gaussian reference, not actual rare-event transfer.'}
 out={'discrete_cases':rows,'real_source_envelope':real,'elapsed_seconds':time.perf_counter()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Exact general formula is proved separately. High Corr(p,logY) does not mean a strong Bernoulli/query channel. Count conditioning beyond exponential tilt is not reduced to a product here.'}
 (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'cases':len(rows),'selected_extreme_rows':[{k:r[k]for k in ['K_coarse_cells','source_minrho','pi','L_target','chosen_c','r_maximal_binary_query_channel','Corr_p_logY_under_Y_anchor']}for r in rows if r['source_minrho']==1e-14 and r['pi']==.5 and r['L_target']=='7/6'],'real':real,'seconds':out['elapsed_seconds']},indent=2))
if __name__=='__main__':main()
