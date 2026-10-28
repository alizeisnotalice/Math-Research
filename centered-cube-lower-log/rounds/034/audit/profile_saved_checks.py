"""Independent arithmetic audit of saved summaries; no source or MC rerun."""
from pathlib import Path
import hashlib,json,math

ROOT=Path(__file__).resolve().parent.parent
SRC=ROOT/'block'
reg=json.loads((SRC/'profile_registration.json').read_text())
saved=json.loads((SRC/'profile_results.json').read_text())
checks={p:hashlib.sha256((SRC/p).read_bytes()).hexdigest()==h
        for p,h in reg['preserved_pilot_files_SHA256'].items()}
checks['profile_probe.py']=hashlib.sha256((SRC/'profile_probe.py').read_bytes()).hexdigest()==reg['script_SHA256']
checks['profile_registration.json']=hashlib.sha256((SRC/'profile_registration.json').read_bytes()).hexdigest()==saved['registration_SHA256']
assert all(checks.values())
def divergence(x,p):
    return (x*math.log(x/p) if x else 0)+((1-x)*math.log((1-x)/(1-p)) if x<1 else 0)
def interval(mean,cap,n,alpha,side):
    x=mean/cap;target=math.log(1/alpha)/n
    if side=='lower':
        if x<=0:return 0.
        lo,hi=0.,x
        for _ in range(100):
            z=(lo+hi)/2
            if divergence(x,z)>target:lo=z
            else:hi=z
        return lo*cap
    if x>=1:return cap
    lo,hi=x,1.
    for _ in range(100):
        z=(lo+hi)/2
        if divergence(x,z)>target:hi=z
        else:lo=z
    return hi*cap

det=saved['deterministic']; hg=det['h_grid']; alpha=saved['alpha_each']
assert len(hg)==88 and hg[0]==0 and hg[-1]>=det['log_H2']
assert math.isclose(alpha*4*len(hg),.05,rel_tol=1e-15)
n=saved['source']['n']
second_logs=[n*math.log(x)for x in det['one_coordinate_second_moment_bounds']]
assert max(abs(a-b)for a,b in zip(second_logs,det['log_second_moment_bounds_by_query']))==0
top=max(second_logs)
logh=top+math.log(math.fsum(math.exp(x-top)for x in second_logs))
maxdiff=0.; models=[]
for model in saved['models']:
    profile=model['profile']; rows=[]
    for p in profile:
        for side in ['lower','upper']:
            v=interval(p[side]['mean'],p['contribution_cap'],model['N'],alpha,side)
            maxdiff=max(maxdiff,abs(v-p[side]['CI']))
    lo=max(p['lower']['CI']for p in profile)
    interval_ups=[]
    for p,qn in zip(profile,profile[1:]):
        ratio=qn['lambda']/p['lambda']
        clipped=min(p['upper']['CI'],p['lambda'],det['H2']/p['lambda'])
        interval_ups.append(min(ratio*clipped,det['H2']/p['lambda']))
    up=max(1.,max(interval_ups))
    assert lo==model['complete_weak_lower_CI']
    assert up==model['complete_weak_upper_CI']
    loc=max(profile,key=lambda x:x['lower']['CI'])
    models.append(dict(G=model['G'],lower=lo,upper=up,
                       lower_argmax_h=loc['h'],lower_argmax_ESS=loc['lower']['ESS'],
                       mass_mean=model['mass_mean'],mass_SE=model['mass_SE_descriptive']))
out=dict(scope='Read saved summaries only; no MC or query-function rerun; no interval FP certification.',
         file_hash_checks=checks,profile_points_per_model=len(hg),CI_family_count=4*len(hg),
         alpha_total=alpha*4*len(hg),max_CI_recalculation_absolute_difference=maxdiff,
         log_H2_recalculation_difference=logh-det['log_H2'],models=models,
         weak_bands_separate=False)
Path(__file__).with_name('profile_saved_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
