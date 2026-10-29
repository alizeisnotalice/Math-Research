"""General-R true-source quantile pilot, streaming all 233 common lengths.

Counts and N*4*R*J sufficient statistics only: never forms 4*R*J*Q features.
Ordinary FP, not a certified finite C_n lower bound. No asymptotic order claim.
"""
import os
for key in ['VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']:os.environ.setdefault(key,'4')
from pathlib import Path
from statistics import NormalDist
from fractions import Fraction as F
import argparse,gc,hashlib,json,math,resource,sys,time
import numpy as np
from helpers import lse,cdf_enclosure,log_prob,log_mixture,tail_enclosures,summarize

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'round040_source'
sys.path.insert(0,str(SOURCE))
from source_tables import Source
from build_cache import exact_grid
SOURCE_HASH='3d49b03dab553e1d828638248a724ee5da9d5201b9eba57a7a8643fa76f3c616'
N=256
RVALUES=[4,8,16,32]
H=np.arange(0,10.00001,.5)
MODELS=['fixed0','joint','prior']
ALPHA=.05/(2*len(RVALUES)*len(MODELS)*len(H))
ND=NormalDist()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rss():return {'raw':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'unit':'bytes'if sys.platform=='darwin'else'KiB','process_lifetime_peak':True}

def make_cuts(src):
    mu=src.n*src.pchi;var=src.n*src.pchi*(1-src.pchi)
    cuts=[int(round(mu+math.sqrt(var)*ND.inv_cdf(b/src.R)))for b in range(1,src.R)]
    assert 0<cuts[0] and cuts[-1]<src.n and cuts==sorted(set(cuts))
    lo,up=cdf_enclosure(np.array(cuts),np.full(len(cuts),mu),np.full(len(cuts),var),
                         np.full(len(cuts),src.n*1e-13),1/math.sqrt(src.n))
    pi_l=np.maximum(0,np.r_[lo,1]-np.r_[0,up])
    pi_u=np.minimum(1,np.r_[up,1]-np.r_[0,lo])
    assert pi_l.min()>0 and np.all(pi_l<=pi_u)
    return cuts,pi_l,pi_u

def feature(src,length,b):
    e=src.evaluate(float(length),b)
    total=e['mean_plus']+e['mean_minus']
    z=np.log(total);p=e['mean_plus']/total
    columns=np.column_stack([z,p,p*(1-p),np.log1p(e['relative_variance_class_max'])])
    nominal=np.exp(z)
    # Proposal is defined by precisely this FP nominal vector. In ideal real
    # arithmetic its normalizer is the displayed mean, not forcibly 1.
    norm=float(nominal.mean())
    error=float(np.max(e['mean_error_abs']));min_s=float(total.min())
    t=src.n*(-math.log1p(-2*error/min_s))
    pe=src.n*error/min_s
    return columns,nominal,norm,t,pe,e

def bounds(stats,src,cuts,pi_l,pi_u,infl,perror):
    z,mu,var,lr=stats
    log_r=log_prob(np.expm1(lr))
    al=z-infl[None,:,:];au=z+infl[None,:,:]
    low_cuts=np.array([-1]+cuts)[None,:,None]
    high_cuts=np.array(cuts+[src.n])[None,:,None]
    ll,lu=cdf_enclosure(low_cuts,mu,var,perror[None,:,:],1/math.sqrt(src.n))
    hl,hu=cdf_enclosure(high_cuts,mu,var,perror[None,:,:],1/math.sqrt(src.n))
    ll[:,0,:]=0;lu[:,0,:]=0;hl[:,-1,:]=1;hu[:,-1,:]=1
    pl=np.maximum(0,hl-lu);pu=np.minimum(1,hu-ll)
    joint=log_mixture(al+log_prob(pl),au+log_prob(pu),log_r)
    prior=log_mixture(al+log_prob(pi_l)[None,:,None],au+log_prob(pi_u)[None,:,None],log_r)
    fixed=(al[:,0,:],au[:,0,:],2*au[:,0,:]+log_r[:,0,:])
    out={'fixed0':fixed,'joint':joint,'prior':prior}
    diag={'minimum_segment_probability_lower':float(pl.min()),'mean_segment_enclosure_width':float((pu-pl).mean()),
          'r_quantiles':np.quantile(np.exp(log_r),[0,.5,.9,.99,1]).tolist()}
    return out,diag

def run_case(R,plan_hash):
    started=time.time();src=Source(R=R);grid=exact_grid(src);J=len(grid)
    cuts,pi_l,pi_u=make_cuts(src)
    seed=202610290400+R;rng=np.random.default_rng(seed)
    counts=np.empty((N,src.Q),np.int64);tags=[];anchor_hash=hashlib.sha256()
    sample_start=time.time()
    for i in range(N):
        if rng.random()<.5:
            probability=np.full(src.Q,1/src.Q);tags.append([-1,-1])
        else:
            b=int(rng.integers(R));j=int(rng.integers(J))
            _,nominal,norm,_,_,_=feature(src,grid[j],b)
            probability=nominal/nominal.sum();tags.append([b,j])
            anchor_hash.update(np.array([b,j],dtype=np.int64).tobytes());anchor_hash.update(nominal.tobytes())
        counts[i]=rng.multinomial(src.n,probability)
    count_hash=hashlib.sha256(counts.tobytes()).hexdigest()
    counts=counts.astype(float)
    draw_seconds=time.time()-sample_start
    stats=np.empty((4,N,R,J),dtype=float)
    norms=np.empty((R,J));infl=np.empty((R,J));perror=np.empty((R,J))
    max_analytic=max_guard=0.
    stream_hash=hashlib.sha256()
    print('BEGIN',R,'N',N,'J',J,'draw_seconds',draw_seconds,flush=True)
    stream_start=time.time()
    for b in range(R):
        for j,length in enumerate(grid):
            cols,_,norm,t,pe,ev=feature(src,length,b)
            values=counts@cols
            stats[:,:,b,j]=values.T
            stream_hash.update(values.tobytes())
            norms[b,j]=norm;infl[b,j]=t;perror[b,j]=pe
            max_analytic=max(max_analytic,float(np.max(ev['analytic_mean_error_abs'])))
            max_guard=max(max_guard,float(np.max(ev['floating_point_guard_abs'])))
        if b%4==3:
            peak=rss();peak_bytes=peak['raw']*(1 if peak['unit']=='bytes'else 1024)
            if peak_bytes>2_000_000_000:raise MemoryError('Registered RSS gate exceeded')
            print('STREAM',R,b+1,'/',R,'seconds',time.time()-stream_start,flush=True)
    stream_seconds=time.time()-stream_start
    del counts;gc.collect()
    C_Z=math.exp(max(0.,src.n*float(np.log(norms).max())))
    C_E=math.exp(float(infl.max()));cap_all=2*R*J*C_Z*C_E
    component_sum=lse((stats[0]-src.n*np.log(norms)[None,:,:]).reshape(N,-1),axis=1)
    logD=np.logaddexp(math.log(.5),math.log(.5)-math.log(R*J)+component_sum)
    weights=np.exp(-logD)
    assert weights.max()<=2+1e-12
    out,diag=bounds(stats,src,cuts,pi_l,pi_u,infl,perror)
    del stats;gc.collect()
    lower=np.empty((3,N,len(H)));upper=np.empty_like(lower)
    models={}
    for mi,model in enumerate(MODELS):
        rows=[]
        for hi,h in enumerate(H):
            B,U=tail_enclosures(*out[model],float(h))
            x=np.exp(h-logD)*B;y=np.exp(h-logD)*U
            cap=min(2*math.exp(h),cap_all)
            assert max(x.max(),y.max())<=cap*(1+1e-8)
            lower[mi,:,hi]=x;upper[mi,:,hi]=y
            rows.append({'h':float(h),'lambda':math.exp(h),'cap':cap,
                         'lower':summarize(x,cap,ALPHA,False),'upper':summarize(y,cap,ALPHA,True)})
        best=max(rows,key=lambda row:row['lower']['mean'])
        ci=max(rows,key=lambda row:row['lower']['CI'])
        models[model]={'profile':rows,'best_lower_MC':best,
                       'best_lower_CI_with_L1_anchor':max(1.,ci['lower']['CI']),
                       'anchor_scale_is_extra_L1':True}
    paired=[]
    for hi,h in enumerate(H):
        delta=lower[1,:,hi]-upper[2,:,hi]
        paired.append({'h':float(h),'joint_lower_minus_prior_upper_MC':float(delta.mean()),
                       'paired_SE_descriptive':float(delta.std(ddof=1)/math.sqrt(N)),
                       'simultaneous_CI_endpoint_gap':models['joint']['profile'][hi]['lower']['CI']-models['prior']['profile'][hi]['upper']['CI']})
    payload=HERE/f'R{R}_paired.npz'
    assert not payload.exists()
    np.savez_compressed(payload,h=H,logD=logD,component_log_sum_normalized=component_sum,
                        lower=lower,upper=upper,tags=np.array(tags,dtype=np.int64))
    result={'plan_SHA256':plan_hash,'R':R,'N':N,'seed':seed,'n':src.n,'K':src.K,'Q':src.Q,'J':J,
            'cuts':cuts,'cuts_definition':'fixed rounded Gaussian-quantile integer cuts, not exact Bin quantiles',
            'actual_eta_nominal':src.pchi,'eta_FP_guard':1e-13,'prior_lower':pi_l.tolist(),'prior_upper':pi_u.tolist(),
            'models':models,'paired_comparisons':paired,'cap_all':cap_all,'C_Z':C_Z,'C_E':C_E,
            'max_logmean_inflation':float(infl.max()),'max_p_sum_error':float(perror.max()),
            'max_analytic_mean_error_abs':max_analytic,'empirical_FP_guard_abs':max_guard,
            'normalizer_max_abs_residual':float(np.abs(norms-1).max()),
            'normalizer_min':float(norms.min()),'normalizer_max':float(norms.max()),
            'stored_statistic_shape':[4,N,R,J],'stored_statistic_bytes':4*N*R*J*8,
            'mass_weight_mean':float(weights.mean()),'mass_weight_SE_descriptive':float(weights.std(ddof=1)/math.sqrt(N)),
            'weight_max':float(weights.max()),'uniform_rows':sum(b<0 for b,j in tags),
            'count_stream_SHA256':count_hash,'sampled_anchor_nominal_stream_SHA256':anchor_hash.hexdigest(),
            'stream_statistics_SHA256':stream_hash.hexdigest(),'paired_payload':payload.name,'paired_payload_SHA256':sha(payload),
            'diagnostics':diag,'draw_seconds':draw_seconds,'stream_seconds':stream_seconds,'total_seconds':time.time()-started,
            'rss':rss(),'full_feature_cache_allocated':False,
            'asymptotic_qualification':{'B':4,'B_ge_v':False,'v_ge_16logn':False,'original_matching_family_claimed':False},
            'scope':'Finite genuine cube conditional tail envelopes only; ordinary FP not interval certified; no C_n log gain claim'}
    (HERE/f'R{R}_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print('DONE',R,'mass',result['mass_weight_mean'],'joint_lower',models['joint']['best_lower_MC']['lower']['mean'],
          'prior_lower',models['prior']['best_lower_MC']['lower']['mean'],'seconds',result['total_seconds'],'rss',result['rss'],flush=True)
    return result

def micro():
    assert F(2803,5000)**2*F(29,50)**2/F(21,100)<1
    src=Source(R=4);cuts,pl,pu=make_cuts(src)
    assert len(cuts)==3 and pl.sum()<=1<=pu.sum()
    # Array cut broadcasting agrees with the scalar legacy interface.
    mu=np.full((2,4,3),src.n*src.pchi);v=np.full_like(mu,src.n*.25);e=np.full((1,4,3),.01)
    c=np.array([-1]+cuts)[None,:,None]
    lo,up=cdf_enclosure(c,mu,v,e,1/math.sqrt(src.n))
    for b in range(4):
        a,z=cdf_enclosure(int(c[0,b,0]),mu[:,b,:],v[:,b,:],e[:,b,:],1/math.sqrt(src.n))
        assert np.array_equal(a,lo[:,b,:])and np.array_equal(z,up[:,b,:])
    return {'BE_rational_constant':True,'vectorized_cut_CDF_matches_scalar':True}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    assert sha(SOURCE/'source_tables.py')==SOURCE_HASH
    checks=micro()
    if args.check:print(json.dumps(checks));return
    if (HERE/'PLAN.json').exists()or any((HERE/f'R{r}_RESULTS.json').exists()for r in RVALUES):raise RuntimeError('Refusing repeated pilot')
    sources=[Source(R=r)for r in RVALUES]
    plan={'R_values':RVALUES,'N_per_R':N,'seeds':{str(r):202610290400+r for r in RVALUES},'models':MODELS,
          'log_thresholds':H.tolist(),'CI_family_count':2*len(RVALUES)*len(MODELS)*len(H),'alpha_total':.05,'alpha_each':ALPHA,
          'grid':list(map(str,exact_grid(sources[0]))),'grid_J':233,
          'source_SHA256':SOURCE_HASH,'builder_SHA256':sha(SOURCE/'build_cache.py'),
          'script_SHA256':sha(Path(__file__)),'helper_SHA256':sha(HERE/'helpers.py'),
          'source_eta_nominal':sources[0].pchi,
          'cuts_and_prior_enclosures':{str(s.R):{'cuts':make_cuts(s)[0],'lower':make_cuts(s)[1].tolist(),'upper':make_cuts(s)[2].tolist()}for s in sources},
          'proposal':'half uniform, half uniform over all R*J nominal exp(log total) product count tilts; exact mixture RN',
          'streaming':'counts N*Q and statistics 4*N*R*J only; no full 4*R*J*Q cache, no raw counts saved',
          'source_only_caps':'after streaming: C_Z=max(1,max norm**n), C_E=exp(max loginflate); cap=min(2lambda,2RJ C_Z C_E); rules frozen before sampling',
          'resource_policy':{'RSS_gate_bytes':2_000_000_000,'gate_check':'each four streamed tables',
                             'max_count_bytes':N*sources[0].Q*8,'max_statistic_bytes':4*N*max(RVALUES)*233*8},
          'paired_statistics':'saved core lower/upper contributions plus logD/tags; paired SE descriptive only, endpoint gaps use the 504 simultaneous intervals',
          'FP_qualification':'analytic mean bounds plus empirical guard; complete FP not outward certified',
          'python':sys.version,'numpy':np.__version__,'micro_checks':checks,
          'status':'all sources,cuts,N,thresholds,cap rules frozen before any tail sampling'}
    (HERE/'PLAN.json').write_text(json.dumps(plan,indent=2)+'\n');plan_hash=sha(HERE/'PLAN.json')
    started=time.time();summaries=[]
    for R in RVALUES:
        d=run_case(R,plan_hash)
        summaries.append({'R':R,'total_seconds':d['total_seconds'],'rss':d['rss'],'results_SHA256':sha(HERE/f'R{R}_RESULTS.json')})
        del d;gc.collect()
    (HERE/'RUN_RECEIPT.json').write_text(json.dumps({'plan_SHA256':plan_hash,'cases':summaries,'total_seconds':time.time()-started},indent=2)+'\n')
    print('ALL DONE',time.time()-started,flush=True)
if __name__=='__main__':main()
