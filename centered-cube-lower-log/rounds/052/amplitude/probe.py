"""Real continuous cosine source / centered interval kernels, same queries and priors.
Monte Carlo calibrations only, followed by a separate exact finite-atom identity check.
"""
from pathlib import Path
from fractions import Fraction as F
import numpy as np,math,json,hashlib,time
HERE=Path(__file__).resolve().parent

def weak(values):
 z=np.sort(values)[::-1];scores=z*np.arange(1,len(z)+1)/len(z);k=int(np.argmax(scores));return float(scores[k]),float(z[k])

def main():
 start=time.perf_counter();rng=np.random.default_rng(20261029052);N=32768;batch=2048;amp=.9;scales=np.array([10/9,7/6,11/9]);alphas=np.array([0.,.2,.35,.5,.65,.8,.95,1.1]);priors={'uniform_positive':np.r_[0.,np.ones(7)/7],'two_endpoints':np.array([0.,.5,0,0,0,0,0,.5]),'high_heavy':np.array([0.,.05,.05,.05,.05,.05,.05,.7]),'constant_plus_high':np.array([.5,0,0,0,0,0,0,.5])};rows=[]
 assert 1-amp*max(alphas)>0
 for freq in [1,3]:
  coeff=amp*np.sin(np.pi*freq*scales)/(np.pi*freq*scales)
  for n in [16,64,128]:
   logs=np.empty((N,len(alphas),3));Sarr=np.empty((N,3));Varr=np.empty((N,3))
   for left in range(0,N,batch):
    right=min(N,left+batch);cos=np.cos(2*np.pi*freq*rng.random((right-left,n)));d=cos[:,:,None]*coeff[None,None,:];Sarr[left:right]=d.sum(axis=1);Varr[left:right]=(d*d).sum(axis=1)
    for j,a in enumerate(alphas):logs[left:right,j,:]=np.log1p(a*d).sum(axis=1)
   T=np.exp(logs);Malpha=T.max(axis=2);individual=[weak(Malpha[:,j])[0]for j in range(len(alphas))];winners=np.argmax(logs[:,1:,:],axis=2);changed=float(np.mean(np.any(winners!=winners[:,:1],axis=1)));linear_winner=np.argmax(Sarr,axis=1);agreement=[float(np.mean(winners[:,j]==linear_winner))for j in range(len(alphas)-1)]
   for name,prior in priors.items():
    mixture=np.einsum('a,qas->qs',prior,T);W,thr=weak(mixture.max(axis=1));avg=float(prior@individual);best=max(individual[j]for j in range(len(alphas))if prior[j]>0)
    rows.append({'frequency':freq,'n':n,'queries':N,'prior':name,'source_min_density':1-amp*float(max(alphas)),'q_L_exact':(n*coeff*coeff/2).tolist(),'q_L_spread':float(n*(coeff.max()**2-coeff.min()**2)/2)if False else float(np.ptp(n*coeff*coeff/2)),'max_query_variance_V_minus_q_observed':float(np.max(np.abs(Varr-n*coeff*coeff/2))),'alpha_winner_scale_changes_observed_fraction':changed,'agreement_with_max_S_winner_by_positive_alpha':agreement,'individual_product_max_weak_MC':individual,'mixture_max_weak_MC':W,'prior_average_product_weak_MC':avg,'best_supported_product_weak_MC':best,'ratio_to_prior_average_MC':W/avg,'ratio_to_best_MC':W/best,'empirical_optimizing_threshold':thr,'classification':'Same-query Monte Carlo of actual continuous cosine density and exact centered kernels. No probability certificate or accepted-source transfer.'})
 # Exact n2 rational example: no single alpha product pointwise dominates the nondegenerate mixture.
 delta=F(1,10);u=F(1,5);v=F(1);mean=(u+v)/2;variance=(v-u)**2/4
 atom=[]
 for signs in [(1,1),(1,-1),(-1,1),(-1,-1)]:
  def product(a):return (1+a*delta*signs[0])*(1+a*delta*signs[1])
  mix=(product(u)+product(v))/2;single=product(mean);atom.append({'signs':signs,'mixture':str(mix),'product_at_mean_alpha':str(single),'difference':str(mix-single)})
 assert sum(F(a['mixture'])for a in atom)/4==1
 assert F(atom[0]['difference'])==delta*delta*variance>0
 out={'source':'rho_alpha(x)=1+alpha*.9*cos(2*pi*m*x); legal amplitudes alpha in [0,1/.9]','scales':['10/9','7/6','11/9'],'amplitudes':alphas.tolist(),'cases':rows,'rational_n2_nondomination_example':{'delta':str(delta),'prior_atoms':[str(u),str(v)],'prior_weights':['1/2','1/2'],'mean_alpha':str(mean),'atoms':atom,'scope':'Exact pointwise nondomination, not a weaknorm counterexample.'},'elapsed_seconds':time.perf_counter()-start,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Numerical exploration precedes proof. Nonflat q_L cosine examples do not satisfy the accepted flat-variance hypothesis.'}
 (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'cases':len(rows),'max_ratio_to_average_MC':max(r['ratio_to_prior_average_MC']for r in rows),'max_ratio_to_best_MC':max(r['ratio_to_best_MC']for r in rows),'winning_scale_change_fractions':sorted(set(r['alpha_winner_scale_changes_observed_fraction']for r in rows)),'q_spreads':sorted(set(r['q_L_spread']for r in rows)),'seconds':out['elapsed_seconds']},indent=2))
if __name__=='__main__':main()
