"""Final fixed-target lower/upper analysis on exactly the round40 sample."""
from pathlib import Path
import hashlib,json,math,sys,time
import numpy as np
from uniform_upper import uniform_upper
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round040_stream';REPLAY=HERE.parent/'round041_replay';MCD=HERE.parent/'round041_mcd';AUDIT=HERE.parent/'round041_audit'
sys.path.insert(0,str(OLD))
from helpers import summarize
MODELS=['fixed0','joint','prior'];ALPHA=.05/1008
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if (HERE/'FINAL_PLAN.json').exists()or(HERE/'FINAL_RESULTS.json').exists():raise RuntimeError('Refusing replacement')
    receipt=json.loads((REPLAY/'RECONSTRUCTION_RECEIPT.json').read_text())
    assert receipt['all_acceptance_checks_passed']
    plan={'script_SHA256':sha(Path(__file__)),'uniform_upper_SHA256':sha(HERE/'uniform_upper.py'),
          'uniform_proof_SHA256':sha(MCD/'uniform_PROOF.md'),'method_uniform_CI_proof_SHA256':sha(AUDIT/'METHOD_UNIFORM_CI.md'),
          'reconstruction_receipt_SHA256':sha(REPLAY/'RECONSTRUCTION_RECEIPT.json'),
          'previous_lower_plan_SHA256':sha(HERE/'PLAN.json'),'CI_alpha_each':ALPHA,
          'reported_one_sided_endpoints':1008,'fixed_target_directions':504,
          'upper_method':'per-scale min of Markov/Cantelli with variance/3 and continuous-uniform Lipschitz tail exp(-6 gap^2/S); then233-scale union',
          'exp_floor':'upper exponent floored upward at-700; zero S handled separately; full FP not certified',
          'posthoc':'Same previously observed data. Analytic pointwise envelopes fixed by geometry; no target,threshold,source orcap change. Not independent confirmation.',
          'diagnostic_only_method':'factor2 McD and variance-only upper means reported without extra endpoint families',
          'inputs':{rec['upper_inputs_file']:sha(REPLAY/rec['upper_inputs_file'])for rec in receipt['cases']}}
    (HERE/'FINAL_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
    cases=[];started=time.time();floor_total=0
    for rec in receipt['cases']:
        R=rec['R'];meta=json.loads((OLD/f'R{R}_RESULTS.json').read_text())
        with np.load(REPLAY/rec['upper_inputs_file'],allow_pickle=False)as inp:
            lm=inp['logmean_upper'];lv=inp['logvariance_old'];sm=inp['max_table_logproduct'];s0=inp['table0_logproduct']
            h=inp['h'];logD=inp['logD'];mcd_upper=inp['new_upper'];var_upper=inp['variance_only_upper']
        with np.load(OLD/f'R{R}_paired.npz',allow_pickle=False)as p:
            oldlower=p['lower'];oldupper=p['upper']
            assert np.array_equal(p['h'],h)and np.array_equal(p['logD'],logD)
        B=np.clip(oldlower*np.exp(logD[None,:,None]-h[None,None,:]),0,1)
        newlower=3*B/(1+2*B)*np.exp(h[None,None,:]-logD[None,:,None])
        newupper=[];floors=0
        for hi,hr in enumerate(h):
            u,fc=uniform_upper(lm,lv,sm,s0,float(hr),R);floors+=fc
            newupper.append(u*np.exp(hr-logD)[None,:])
        newupper=np.stack(newupper,axis=2)
        assert np.max(newlower-newupper)<1e-8
        assert np.max(newupper-mcd_upper)<1e-8 and np.max(newupper-oldupper)<1e-8
        models={}
        for mi,name in enumerate(MODELS):
            rows=[]
            for hi,hr in enumerate(h):
                cap=min(2*math.exp(hr),meta['cap_all'])
                assert np.max(newupper[mi,:,hi])<=cap+1e-8
                rows.append({'h':float(hr),'lambda':math.exp(hr),'cap':cap,
                    'old_lower':summarize(oldlower[mi,:,hi],cap,ALPHA,False),
                    'sharp_lower':summarize(newlower[mi,:,hi],cap,ALPHA,False),
                    'old_upper':summarize(oldupper[mi,:,hi],cap,ALPHA,True),
                    'sharp_upper':summarize(newupper[mi,:,hi],cap,ALPHA,True),
                    'mcd2_upper_mean_descriptive':float(mcd_upper[mi,:,hi].mean()),
                    'variance_only_upper_mean_descriptive':float(var_upper[mi,:,hi].mean())})
            models[name]={'profile':rows,'best_sharp_lower':max(rows,key=lambda r:r['sharp_lower']['mean']),
                          'best_sharp_lower_CI_with_extra_L1':max(1.,max(r['sharp_lower']['CI']for r in rows))}
        pairs=[]
        for hi,hr in enumerate(h):
            delta=newlower[1,:,hi]-newupper[2,:,hi]
            pairs.append({'h':float(hr),'sharp_joint_lower_minus_sharp_prior_upper':float(delta.mean()),
                          'paired_SE_descriptive':float(delta.std(ddof=1)/math.sqrt(meta['N'])),
                          'simultaneous_endpoint_gap':models['joint']['profile'][hi]['sharp_lower']['CI']-models['prior']['profile'][hi]['sharp_upper']['CI']})
        cases.append({'R':R,'N':meta['N'],'models':models,'pairs':pairs,
                      'positive_same_h_MC':sum(r['sharp_joint_lower_minus_sharp_prior_upper']>0 for r in pairs),
                      'positive_same_h_endpoint_gaps':sum(r['simultaneous_endpoint_gap']>0 for r in pairs),
                      'upper_exp_floor_applications':floors,
                      'pointwise_lower_minus_upper_max':float((newlower-newupper).max()),
                      'pointwise_new_upper_minus_mcd_max':float((newupper-mcd_upper).max())})
        floor_total+=floors
    out={'plan_SHA256':sha(HERE/'FINAL_PLAN.json'),'cases':cases,'total_reported_endpoints':1008,
         'underlying_target_direction_events':504,'alpha_each':ALPHA,'new_independent_samples':0,
         'upper_exp_floor_count':floor_total,'seconds':time.time()-started,
         'full_FP_certificate':False,'Cn_log_gain_proved':False}
    (HERE/'FINAL_RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
    for c in cases:
        row=c['models']['joint']['best_sharp_lower'];p=next(r for r in c['models']['prior']['profile']if r['h']==row['h'])
        print(c['R'],'h',row['h'],'lower',row['sharp_lower']['mean'],'prior_upper',p['sharp_upper']['mean'],
              'positive_pairs',c['positive_same_h_MC'],'CIpositive',c['positive_same_h_endpoint_gaps'])
    print('seconds',out['seconds'],'floors',floor_total)
if __name__=='__main__':main()
