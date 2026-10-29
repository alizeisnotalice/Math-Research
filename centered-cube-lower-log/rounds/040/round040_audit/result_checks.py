"""Independent saved payload/RN/cap/paired checks. No MC, no source cache rebuild."""
from pathlib import Path
import json,math,hashlib,sys,argparse
import numpy as np
HERE=Path(__file__).resolve().parent;OTHER=HERE.parent/'round040_stream'
parser=argparse.ArgumentParser()
parser.add_argument('--replace-own-metadata-review',action='store_true')
args=parser.parse_args()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mean(x):return math.fsum(map(float,x))/len(x)
def se(x):
    m=mean(x);return math.sqrt(math.fsum((float(v)-m)**2 for v in x)/(len(x)-1)/len(x))
def ess(x):
    m=max(x)
    if m<=0:return 0.
    vals=[float(v)/m for v in x];return math.fsum(vals)**2/math.fsum(v*v for v in vals)
plan=json.loads((OTHER/'PLAN.json').read_text());receipt=json.loads((OTHER/'RUN_RECEIPT.json').read_text())
planhash=sha(OTHER/'PLAN.json')
assert receipt['plan_SHA256']==planhash
for name,key,folder in [('probe.py','script_SHA256',OTHER),('helpers.py','helper_SHA256',OTHER),
                        ('source_tables.py','source_SHA256',HERE.parent/'round040_source'),
                        ('build_cache.py','builder_SHA256',HERE.parent/'round040_source')]:
    assert sha(folder/name)==plan[key]
prior=json.loads((HERE/'pilot_checks.json').read_text())
assert all(prior['reviewed_SHA256'][name]==sha(OTHER/name)for name in ['probe.py','helpers.py'])
assert plan['CI_family_count']==504 and plan['N_per_R']==256 and plan['grid_J']==233
rows=[];maxmean=maxse=maxess=maxpair=0.
for R,rec in zip([4,8,16,32],receipt['cases']):
    d=json.loads((OTHER/f'R{R}_RESULTS.json').read_text());p=OTHER/d['paired_payload']
    assert d['plan_SHA256']==planhash and sha(p)==d['paired_payload_SHA256']
    assert rec['R']==R and rec['results_SHA256']==sha(OTHER/f'R{R}_RESULTS.json')
    assert d['R']==R and d['N']==256 and d['J']==233 and d['n']==6704800
    h=np.array(plan['log_thresholds']);N=d['N'];J=d['J']
    with np.load(p,allow_pickle=False)as data:
        assert np.array_equal(data['h'],h)
        lo=data['lower'];up=data['upper'];logD=data['logD'];cs=data['component_log_sum_normalized'];tags=data['tags']
        assert lo.shape==up.shape==(3,N,21)and tags.shape==(N,2)
        expected=np.logaddexp(math.log(.5),math.log(.5)-math.log(R*J)+cs)
        derr=float(abs(expected-logD).max());assert derr==0.
        assert np.all(np.isfinite(lo))and np.all(np.isfinite(up))and np.all(lo>=0)and np.all(up>=lo-1e-12)
        assert np.all((tags[:,0]==-1)==(tags[:,1]==-1))
        tilted=tags[:,0]>=0
        assert np.all(tags[tilted,0]<R)and np.all(tags[tilted,1]<J)
        assert int((~tilted).sum())==d['uniform_rows']
        w=np.exp(-logD);assert w.max()<=2+1e-12
        assert abs(mean(w)-d['mass_weight_mean'])<1e-13
        assert abs(se(w)-d['mass_weight_SE_descriptive'])<1e-13
        ce=math.exp(d['max_logmean_inflation']);assert ce==d['C_E']
        norm_max=d['normalizer_max'];norm_min=d['normalizer_min']
        assert max(abs(norm_min-1),abs(norm_max-1))==d['normalizer_max_abs_residual']
        cz=math.exp(max(0.,d['n']*math.log(norm_max)))
        assert abs(cz-d['C_Z'])<1e-15
        fullcap=2*R*J*d['C_Z']*d['C_E'];assert fullcap==d['cap_all']
        maxratio=0.
        for mi,model in enumerate(plan['models']):
            for hi,hr in enumerate(h):
                r=d['models'][model]['profile'][hi]
                cap=min(2*math.exp(hr),fullcap)
                assert r['h']==hr and r['cap']==cap
                for array,key in [(lo,'lower'),(up,'upper')]:
                    values=array[mi,:,hi]
                    maxratio=max(maxratio,float(values.max()/cap))
                    assert values.max()<=cap*(1+1e-8)
                    maxmean=max(maxmean,abs(mean(values)-r[key]['mean']))
                    maxse=max(maxse,abs(se(values)-r[key]['SE_descriptive']))
                    maxess=max(maxess,abs(ess(values)-r[key]['ESS']))
                    assert abs(mean(values)-r[key]['mean'])<1e-10
                    assert abs(ess(values)-r[key]['ESS'])<1e-8
        differences=[];endpoint=[]
        for hi,hr in enumerate(h):
            delta=lo[1,:,hi]-up[2,:,hi]
            record=d['paired_comparisons'][hi]
            gap=d['models']['joint']['profile'][hi]['lower']['CI']-d['models']['prior']['profile'][hi]['upper']['CI']
            maxpair=max(maxpair,abs(mean(delta)-record['joint_lower_minus_prior_upper_MC']),
                        abs(se(delta)-record['paired_SE_descriptive']))
            assert abs(gap-record['simultaneous_CI_endpoint_gap'])<1e-13
            assert mean(delta)<0 and gap<0
            differences.append((float(hr),mean(delta),se(delta)))
            endpoint.append((float(hr),gap))
        rows.append({'R':R,'N':N,'J':J,'uniform_rows':d['uniform_rows'],
                     'same_h_MC_positive_count':0,'same_h_endpoint_positive_count':0,
                     'largest_same_h_MC_difference':max(differences,key=lambda x:x[1]),
                     'largest_endpoint_difference':max(endpoint,key=lambda x:x[1]),
                     'logD_reconstruction_max_abs_error':derr,
                     'max_contribution_over_cap':maxratio,'C_Z':cz,'C_E':ce,
                     'saved_abs_normalizer_residual':d['normalizer_max_abs_residual'],
                     'saved_normalizer_min':norm_min,'saved_normalizer_max':norm_max,
                     'joint_best_lower':d['models']['joint']['best_lower_MC']['lower']['mean'],
                     'prior_best_lower':d['models']['prior']['best_lower_MC']['lower']['mean'],
                     'mass_weight_mean':d['mass_weight_mean'],'payload_SHA256':sha(p)})
out={'cases':rows,'CI_family_count':504,'same_h_comparisons':84,
     'max_mean_fsum_difference':maxmean,'max_SE_fsum_difference':maxse,
     'max_ESS_scaled_fsum_difference':maxess,'max_paired_stat_fsum_difference':maxpair,
     'all_frozen_hashes_match':True,'source_CZ_consistency_method':'CZ recomputed directly as exp(max(0,n log(saved normalizer_max))); min/max reproduce saved max absolute residual',
     'own_check_revision':'Initial metadata-only assertion treated R16 max absolute residual as positive signed residual and stopped before output. Final review reads the saved signed min/max directly; zero new source or MC samples.',
     'runtime_seconds_receipt':receipt['total_seconds'],
     'max_RSS_bytes':max(r['rss']['raw']for r in receipt['cases']),
     'scope':'Zero new samples/source recomputation; independent payload arithmetic; no full CI inversion (separate reviewer); FP certification not claimed'}
target=HERE/'result_checks.json'
if target.exists():
    if not args.replace_own_metadata_review:raise RuntimeError('Refusing replacement without explicit own-review revision flag')
    out['superseded_own_review_SHA256']=sha(target)
target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
