"""All 1008 final CI inversions in Decimal; upstream upper formula reviewed separately."""
from pathlib import Path
from decimal import Decimal as D,getcontext
import json,hashlib,math,sys,importlib.util
import numpy as np
HERE=Path(__file__).resolve().parent;POST=HERE.parent/'round041_post';REP=HERE.parent/'round041_replay';OLD=HERE.parent/'round040_stream'
getcontext().prec=60;sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dec(x):return D.from_float(float(x))
def endpoint(mean,cap,n,upper):
    x=mean/cap;target=D(20160).ln()/D(n)
    if x==0:return cap*(1-(-target).exp())if upper else D(0)
    if x==1:return cap if upper else cap*(-target).exp()
    assert 0<x<1;a,b=(x,D(1))if upper else(D(0),x)
    ent=x*x.ln()+(1-x)*(1-x).ln()
    for _ in range(100):
        p=(a+b)/2;kl=ent-x*p.ln()-(1-x)*(1-p).ln()
        if upper:
            if kl>target:b=p
            else:a=p
        else:
            if kl>target:a=p
            else:b=p
    return cap*(a+b)/2
plan=json.loads((POST/'FINAL_PLAN.json').read_text());result=json.loads((POST/'FINAL_RESULTS.json').read_text())
assert result['plan_SHA256']==sha(POST/'FINAL_PLAN.json')
for field,path in [('script_SHA256',POST/'aggregate.py'),('uniform_upper_SHA256',POST/'uniform_upper.py'),
 ('uniform_proof_SHA256',HERE.parent/'round041_mcd/uniform_PROOF.md'),
 ('method_uniform_CI_proof_SHA256',HERE/'METHOD_UNIFORM_CI.md'),
 ('reconstruction_receipt_SHA256',REP/'RECONSTRUCTION_RECEIPT.json'),('previous_lower_plan_SHA256',POST/'PLAN.json')]:assert plan[field]==sha(path)
for name,h in plan['inputs'].items():assert sha(REP/name)==h
for name,h in json.loads((POST/'PLAN.json').read_text())['input_SHA256'].items():assert sha(OLD/name)==h
spec=importlib.util.spec_from_file_location('reviewed_upper6',POST/'uniform_upper.py')
uppermod=importlib.util.module_from_spec(spec);spec.loader.exec_module(uppermod)
count=0;maxci=maxmean=maxpair=0.;rows=[];floor_total=0
for case in result['cases']:
    R=case['R'];meta=json.loads((OLD/f'R{R}_RESULTS.json').read_text())
    with np.load(REP/f'R{R}_upper_inputs.npz',allow_pickle=False)as p:
        lm=p['logmean_upper'];lv=p['logvariance_old'];sm=p['max_table_logproduct'];s0=p['table0_logproduct'];h=p['h'];logD=p['logD']
    with np.load(OLD/f'R{R}_paired.npz',allow_pickle=False)as p:
        oldlower=p['lower'];oldupper=p['upper'];assert np.array_equal(h,p['h'])and np.array_equal(logD,p['logD'])
    n=meta['N'];density=[dec(z).exp()for z in logD];local_ci=0.;floors=0;sharplow=np.zeros_like(oldlower);sharpup=np.zeros_like(oldupper)
    for hi,hr in enumerate(h):
        conditional,fc=uppermod.uniform_upper(lm,lv,sm,s0,float(hr),R);floors+=fc
        lam=dec(hr).exp()
        for mi,name in enumerate(['fixed0','joint','prior']):
            row=case['models'][name]['profile'][hi];cap=dec(row['cap'])
            lists=[[],[],[],[]]
            for i in range(n):
                x=dec(oldlower[mi,i,hi]);b=x*density[i]/lam
                sharp=3*x/(1+2*b)if x else D(0)
                upper=dec(conditional[mi,i])*lam/density[i]
                lists[0].append(x);lists[1].append(sharp);lists[2].append(dec(oldupper[mi,i,hi]));lists[3].append(upper)
                assert sharp<=upper+D('1e-9')and upper<=cap+D('1e-9')
                sharplow[mi,i,hi]=float(sharp);sharpup[mi,i,hi]=float(upper)
            for vals,key,up in zip(lists,['old_lower','sharp_lower','old_upper','sharp_upper'],[False,False,True,True]):
                avg=sum(vals,D(0))/n;ci=endpoint(avg,cap,n,up);count+=1
                err=abs(float(ci)-row[key]['CI']);maxci=max(maxci,err);local_ci=max(local_ci,err)
                maxmean=max(maxmean,abs(float(avg)-row[key]['mean']));assert err<1e-8
    pairs=[]
    for hi,hr in enumerate(h):
        delta=sharplow[1,:,hi]-sharpup[2,:,hi];avg=math.fsum(map(float,delta))/n
        old=case['pairs'][hi];maxpair=max(maxpair,abs(avg-old['sharp_joint_lower_minus_sharp_prior_upper']))
        assert abs(avg-old['sharp_joint_lower_minus_sharp_prior_upper'])<1e-10
        gap=case['models']['joint']['profile'][hi]['sharp_lower']['CI']-case['models']['prior']['profile'][hi]['sharp_upper']['CI']
        assert gap==old['simultaneous_endpoint_gap']and avg<0 and gap<0
        pairs.append((float(hr),avg,gap))
    assert floors==case['upper_exp_floor_applications'];floor_total+=floors
    rows.append({'R':R,'N':n,'CI_count':252,'max_Decimal_CI_difference':local_ci,
                 'largest_same_h_MC_gap':max(pairs,key=lambda r:r[1]),
                 'largest_same_h_CI_gap':max(pairs,key=lambda r:r[2]),'positive_same_h_count':0,
                 'upward_exponent_floor_applications':floors})
assert count==1008 and floor_total==result['upper_exp_floor_count']
out={'Decimal_precision':60,'bisection_iterations':100,'alpha_exact':'1/20160','CI_count':count,
     'max_Decimal_CI_difference':maxci,'max_mean_difference':maxmean,'max_paired_mean_difference':maxpair,
     'cases':rows,'all_hashes_match':True,'new_independent_samples':0,
     'upper_function_scope':'Calls frozen reviewed factor6 formula; independent per-scale formula audit delegated to other reviewer. CI inversion/exp weighting/sharp-lower recovered independently in Decimal.',
     'qualification':'No certificate for upstream source/BLAS/erfc/pointwise bounds; no new MC'}
p=HERE/'final_ci_checks.json'
if p.exists():raise RuntimeError('Refusing replacement')
p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
