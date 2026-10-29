"""Deterministic variance/3 lower sharpening of the fixed round40 samples."""
from pathlib import Path
import hashlib,json,math,sys,time
import numpy as np

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round040_stream'
sys.path.insert(0,str(OLD))
from helpers import summarize
MODELS=['fixed0','joint','prior']
RVALUES=[4,8,16,32]
FAMILY=1008 # 756 used now; reserve252 for a separately proved new upper.
ALPHA=.05/FAMILY
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=HERE/'RESULTS.json'
    if target.exists()or(HERE/'PLAN.json').exists():raise RuntimeError('Refusing replacement')
    oldplan=json.loads((OLD/'PLAN.json').read_text())
    inputs={f'R{r}_{name}':sha(OLD/f'R{r}_{name}')for r in RVALUES for name in ['RESULTS.json','paired.npz']}
    plan={'old_plan_SHA256':sha(OLD/'PLAN.json'),'input_SHA256':inputs,'script_SHA256':sha(Path(__file__)),
          'models':MODELS,'R':RVALUES,'thresholds':oldplan['log_thresholds'],
          'map':'recover B=x*exp(logD-h), preserve zero, sharpen conditional B by3B/(1+2B), multiply back by exp(h-logD)',
          'variance_factor':'1/3 fixed analytically from uniform-bin Lipschitz proof; no c optimization',
          'CI_family_reserved':FAMILY,'functions_used_now':756,'alpha_each':ALPHA,'alpha_total_budget':.05,
          'reserved_upper_functions':252,'upper_status':'unchanged round40 upper; a new upper needs separate proof and old-sample per-scale reconstruction',
          'posthoc_disclosure':'Round40 summaries and a descriptive hypothetical-c1/3 preview had been seen; not a blind or independent confirmation.',
          'sample_policy':'No new samples, no seed changes, no source/kernel reconstruction',
          'FP_qualification':'Ordinary FP and original mean/BE premises, no outward certification',
          'python':sys.version,'numpy':np.__version__}
    (HERE/'PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
    started=time.time();cases=[]
    for R in RVALUES:
        d=json.loads((OLD/f'R{R}_RESULTS.json').read_text())
        with np.load(OLD/f'R{R}_paired.npz',allow_pickle=False)as p:
            h=p['h'];logD=p['logD'];lower=p['lower'];upper=p['upper']
        assert list(h)==plan['thresholds']
        B=lower*np.exp(logD[None,:,None]-h[None,None,:])
        assert B.min()>=0 and B.max()<=1+1e-10
        B=np.clip(B,0,1)
        new=3*B/(1+2*B)*np.exp(h[None,None,:]-logD[None,:,None])
        assert np.min(new-lower)>-1e-10 and np.max(new-upper)<1e-10
        models={}
        for mi,name in enumerate(MODELS):
            rows=[]
            for hi,hr in enumerate(h):
                cap=min(2*math.exp(hr),d['cap_all'])
                assert new[mi,:,hi].max()<=cap+1e-8
                rows.append({'h':float(hr),'lambda':math.exp(hr),'cap':cap,
                    'old_lower':summarize(lower[mi,:,hi],cap,ALPHA,False),
                    'sharp_lower':summarize(new[mi,:,hi],cap,ALPHA,False),
                    'old_upper':summarize(upper[mi,:,hi],cap,ALPHA,True)})
            models[name]={'profile':rows,'best_old_lower':max(rows,key=lambda x:x['old_lower']['mean']),
                          'best_sharp_lower':max(rows,key=lambda x:x['sharp_lower']['mean']),
                          'best_sharp_lower_CI_with_extra_L1':max(1.,max(x['sharp_lower']['CI']for x in rows))}
        pairs=[]
        for hi,hr in enumerate(h):
            delta=new[1,:,hi]-upper[2,:,hi]
            pairs.append({'h':float(hr),'sharp_joint_lower_minus_old_prior_upper':float(delta.mean()),
                          'paired_SE_descriptive':float(delta.std(ddof=1)/math.sqrt(d['N'])),
                          'simultaneous_endpoint_gap':models['joint']['profile'][hi]['sharp_lower']['CI']-models['prior']['profile'][hi]['old_upper']['CI']})
        cases.append({'R':R,'N':d['N'],'models':models,'pairs':pairs,
                      'positive_same_h_MC':sum(p['sharp_joint_lower_minus_old_prior_upper']>0 for p in pairs),
                      'positive_endpoint_gaps':sum(p['simultaneous_endpoint_gap']>0 for p in pairs),
                      'raw_B_min':float(B.min()),'raw_B_max':float(B.max()),
                      'max_pointwise_new_lower_minus_old_upper':float((new-upper).max())})
    out={'plan_SHA256':sha(HERE/'PLAN.json'),'cases':cases,'seconds':time.time()-started,
         'functions_recomputed':756,'CI_reserved_family':FAMILY,
         'new_tail_samples':0,'independent_confirmation':False,'Cn_log_gain':False}
    target.write_text(json.dumps(out,indent=2)+'\n')
    for c in cases:
        best=c['models']['joint']['best_sharp_lower']
        print(c['R'],'joint sharp',best['sharp_lower']['mean'],'h',best['h'],'positive pairs',c['positive_same_h_MC'],flush=True)
    print('seconds',out['seconds'])
if __name__=='__main__':main()
