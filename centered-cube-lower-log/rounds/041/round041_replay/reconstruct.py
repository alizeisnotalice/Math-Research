"""Non-mutating observer replay of round40, retaining minimal upper statistics.

Same seeds/calls, never a new independent MC. Count and source-statistic stream
hashes must match round40 before the reconstructed upper is accepted.
"""
import os
for key in ['VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']:os.environ.setdefault(key,'4')
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,shutil,sys,time
import numpy as np
from upper_math import conditional_uppers,micro_checks,C_LOG

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round040_stream'
POST=HERE.parent/'round041_post'
MCD=HERE.parent/'round041_mcd'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    micro=micro_checks()
    if args.check:print(json.dumps(micro));return
    if (HERE/'PLAN.json').exists()or(HERE/'base_replay').exists():raise RuntimeError('Refusing repeated reconstruction')
    assert (MCD/'PROOF.md').exists(),'Mathematical audit required before replay'
    oldplan=json.loads((OLD/'PLAN.json').read_text());postplan=json.loads((POST/'PLAN.json').read_text())
    assert postplan['CI_family_reserved']==1008
    plan={'old_plan_SHA256':sha(OLD/'PLAN.json'),'postprocessing_plan_SHA256':sha(POST/'PLAN.json'),
          'base_driver_SHA256':sha(OLD/'probe.py'),'observer_SHA256':sha(Path(__file__)),
          'upper_math_SHA256':sha(HERE/'upper_math.py'),'MCD_proof_SHA256':sha(MCD/'PROOF.md'),
          'variance_factor':'1/3','source_r0':'1/1000000','C_log_bound':C_LOG,
          'new_upper':'per-scale min of old Markov/Cantelli, variance/3 Cantelli, component-union McDiarmid; then same233-scale union',
          'targets':{'R':oldplan['R_values'],'models':oldplan['models'],'h':oldplan['log_thresholds']},
          'budget':'1008 reserved one-sided functions across old/sharp lower and old/new upper; source-only old caps retained',
          'observer':'Only feature max-r checks and bounds side outputs; no writes/mutations to stats or returned old bounds; no added RNG calls',
          'acceptance':'Exact count, sampled-anchor and source-statistic stream hashes, plus old contribution replay checks',
          'posthoc_disclosure':'Same previously observed round40 sample; diagnostic reconstruction, no independent confirmation',
          'FP_qualification':'Original mean/BE/FP qualifications persist; this is not a finite C_n certificate',
          'micro_checks':micro,'python':sys.version,'numpy':np.__version__}
    (HERE/'PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
    sys.path.insert(0,str(OLD))
    spec=importlib.util.spec_from_file_location('frozen40base',OLD/'probe.py')
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    dest=HERE/'base_replay';dest.mkdir();shutil.copy2(OLD/'helpers.py',dest/'helpers.py')
    base.HERE=dest
    original_feature=base.feature;original_bounds=base.bounds
    max_r={R:0. for R in oldplan['R_values']};side={}
    def observe_feature(src,*args):
        value=original_feature(src,*args)
        r=float(np.max(value[-1]['relative_variance_class_max']))
        assert r<=1e-6
        max_r[src.R]=max(max_r[src.R],r)
        return value
    def observe_bounds(stats,src,*args):
        out,diag=original_bounds(stats,src,*args)
        logmu=np.stack([out[name][1]for name in base.MODELS])
        logv=np.stack([out[name][2]for name in base.MODELS])
        maxlog=stats[3].max(axis=1).copy();firstlog=stats[3][:,0,:].copy()
        oldresult=json.loads((OLD/f'R{src.R}_RESULTS.json').read_text())
        with np.load(OLD/f'R{src.R}_paired.npz',allow_pickle=False)as p:
            h=p['h'];logD=p['logD'];oldlower=p['lower'];oldupper=p['upper']
        new=[];variance=[];max_old_diff=0.
        for j,hr in enumerate(h):
            oldU,varU,newU=conditional_uppers(logmu,logv,maxlog,firstlog,float(hr),src.R)
            multiplier=np.exp(hr-logD)[None,:]
            max_old_diff=max(max_old_diff,float(abs(oldU*multiplier-oldupper[:,:,j]).max()))
            assert max_old_diff<1e-9
            new.append(newU*multiplier);variance.append(varU*multiplier)
        new=np.stack(new,axis=2);variance=np.stack(variance,axis=2)
        B=oldlower*np.exp(logD[None,:,None]-h[None,None,:])
        sharpened=3*B/(1+2*B)*np.exp(h[None,None,:]-logD[None,:,None])
        assert np.max(sharpened-new)<1e-9
        assert np.max(new-oldupper)<1e-9 and np.max(new-variance)<1e-9
        # These are the minimal per-scale means/variance and maximum table
        # logproducts needed to independently recompute every new upper.
        target=HERE/f'R{src.R}_upper_inputs.npz'
        np.savez_compressed(target,logmean_upper=logmu,logvariance_old=logv,
                            max_table_logproduct=maxlog,table0_logproduct=firstlog,
                            h=h,logD=logD,new_upper=new,variance_only_upper=variance)
        side[src.R]={'upper_inputs_file':target.name,'upper_inputs_SHA256':sha(target),
                     'max_old_upper_reconstruction_difference':max_old_diff}
        return out,diag
    base.feature=observe_feature;base.bounds=observe_bounds
    started=time.time();base.main()
    acceptance=[]
    for R in oldplan['R_values']:
        old=json.loads((OLD/f'R{R}_RESULTS.json').read_text())
        replay=json.loads((dest/f'R{R}_RESULTS.json').read_text())
        for key in ['count_stream_SHA256','sampled_anchor_nominal_stream_SHA256','stream_statistics_SHA256']:
            assert replay[key]==old[key],(R,key)
        with np.load(OLD/f'R{R}_paired.npz',allow_pickle=False)as a,np.load(dest/f'R{R}_paired.npz',allow_pickle=False)as b:
            for key in a.files:assert np.array_equal(a[key],b[key]),(R,key)
        acceptance.append({'R':R,**side[R],'count_stream_matches':True,'sampled_anchor_stream_matches':True,
                           'source_statistics_stream_matches':True,'all_old_payload_arrays_identical':True,
                           'max_checked_source_r':max_r[R],'replay_seconds':replay['total_seconds'],'replay_RSS':replay['rss']})
    result={'plan_SHA256':sha(HERE/'PLAN.json'),'cases':acceptance,'total_seconds':time.time()-started,
            'new_independent_samples':0,'same_protocol_replay':True,'all_acceptance_checks_passed':True}
    (HERE/'RECONSTRUCTION_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
    print('RECONSTRUCTION ACCEPTED',result['total_seconds'],flush=True)
if __name__=='__main__':main()
