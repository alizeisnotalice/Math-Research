#!/usr/bin/env python3
"""One registered numerical top3 pilot; no rigorous integration claim."""
from pathlib import Path
from itertools import product,combinations,permutations
from collections import Counter
import hashlib,importlib.util,json,os,sys,time
import numpy as np
from scipy.special import betaincinv,gamma
from scipy.stats import qmc
from flint import ctx

HERE=Path(__file__).resolve().parent
MOMENT_SOURCE=HERE.parent/'upper_order002'/'top2_certificate.py'
spec=importlib.util.spec_from_file_location('certified_stable_moments',MOMENT_SOURCE)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def determinant(matrix):
    n=len(matrix);ans=0
    for p in permutations(range(n)):
        inversions=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        term=(-1)**inversions
        for i,j in enumerate(p):term*=matrix[i][j]
        ans+=term
    return ans

def gram(vectors):
    return [[sum(a*b for a,b in zip(v,w)) for w in vectors] for v in vectors]

def geometry():
    directions=[(1,)+t for t in product((-1,1),repeat=4)];x=(1,1,1,1,0)
    classes=Counter();ranks=Counter();gram_dets=Counter()
    for triple in combinations(directions,3):
        g=[];q=[]
        for mask in range(8):
            vectors=[triple[j] for j in range(3) if mask>>j&1]
            g.append(determinant(gram(vectors)) if vectors else 1)
            q.append(determinant(gram(vectors+[x])))
        assert g[7]>0 and min(q)>=0
        ranks[3]+=1;gram_dets[g[7]]+=1;variants=[]
        for p in permutations(range(3)):
            masks=[sum(1<<p[j] for j in range(3) if mask>>j&1) for mask in range(8)]
            gg=[g[k] for k in masks];qq=[q[k] for k in masks]
            variants.append((gg[3],gg[5],gg[6],gg[7],qq[1],qq[2],qq[4],qq[3],qq[5],qq[6],qq[7]))
        classes[min(variants)]+=1
    assert len(classes)==28 and sum(classes.values())==560 and min(gram_dets)==48
    return classes,gram_dets

def polynomial(c,z,key,augmented=False):
    A=16*c;h=z-c[:,None]
    h1,h2,h3=h.T;p12=h1*h2;p13=h1*h3;p23=h2*h3;p123=p12*h3
    if not augmented:
        return A**3+5*A*A*(h1+h2+h3)+A*(key[0]*p12+key[1]*p13+key[2]*p23)+key[3]*p123
    return 4*A**3+A*A*(key[4]*h1+key[5]*h2+key[6]*h3)+A*(key[7]*p12+key[8]*p13+key[9]*p23)+key[10]*p123

def kernel(m,a,z,key):
    Dm=polynomial(m,z,key);Da=polynomial(a,z,key);Na=polynomial(a,z,key,True)
    Q=Na/(16*a*Da)
    assert np.all(Dm>0) and np.all(Q>0)
    return (16*m)**(-1)*Dm**(-.5)*Q**(-1.75)

def main():
    start=time.monotonic();ctx.prec=256;pid=os.getpid()
    hashes={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'moment_source_sha256':hashlib.sha256(MOMENT_SOURCE.read_bytes()).hexdigest()}
    state={'status':'RUNNING_NUMERICAL_PILOT','pid':pid,**hashes,'next_table_index':0,'next_class_index':0,'started_unix':time.time()}
    ngrid=4096;N=16384;logz=np.linspace(np.log(.25),40*np.log(2),ngrid+1)
    table=np.full((ngrid+1,3),np.nan);accum=np.zeros(N)
    if '--resume' in sys.argv:
        old=json.loads((HERE/'checkpoint.json').read_text())
        assert all(old[k]==v for k,v in hashes.items())
        cache=np.load(HERE/'checkpoint_arrays.npz');table=cache['table'];accum=cache['accum']
        state['next_table_index']=old['next_table_index'];state['next_class_index']=old['next_class_index']
    (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps({'pid':pid,'status':'started','source_sha256':hashes['source_sha256']}),flush=True)
    def checkpoint(phase):
        state.update(status=phase,seconds_this_checkpoint=time.monotonic()-start)
        np.savez_compressed(HERE/'checkpoint_arrays.npz',table=table,accum=accum)
        (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
        if time.monotonic()-start>120:
            state['status']='INCOMPLETE_SAVED_RESUMABLE'
            (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
            raise SystemExit(75)
    classes,gram_dets=geometry();keys=sorted(classes)
    for i in range(state['next_table_index'],ngrid+1):
        # Table arguments are fixed positive floating numbers, each converted
        # to an exact Arb point. The interpolation later remains numerical.
        value=.25 if i==0 else float(np.exp(logz[i]))
        table[i]=[float(v) for v in mod.stats(mod.arb(value))]
        state['next_table_index']=i+1
        if i%512==0:checkpoint('BUILDING_FIXED_TABLE')
    assert np.all(np.diff(table[:,0])>0)
    points=qmc.Sobol(d=3,scramble=True,seed=20261026).random_base2(14)
    minimum=betaincinv(14,3,points[:,0])
    upper=minimum[:,None]+(1-minimum[:,None])*points[:,1:]
    uniforms=np.column_stack((minimum,np.sort(upper,axis=1)))
    assert uniforms.min()>table[0,0] and uniforms.max()<table[-1,0], 'No quantile extrapolation allowed'
    z=np.exp(np.interp(uniforms,table[:,0],logz))
    m=np.exp(np.interp(np.log(z[:,0]),logz,np.log(table[:,1])))
    a=np.exp(np.interp(np.log(z[:,0]),logz,np.log(table[:,2])))
    for j in range(state['next_class_index'],len(keys)):
        key=keys[j];value=np.zeros(N)
        for p in permutations(range(3)):value+=kernel(m,a,z[:,p],key)
        accum+=classes[key]*value
        state['next_class_index']=j+1
        checkpoint('ACCUMULATING_GEOMETRY_CLASSES')
    kappa=3*gamma(1.75)/(2**3.5*np.pi**2.5)
    values=kappa*accum/3360
    blocks=np.array_split(values,8);means=[float(v.mean()) for v in blocks]
    result={'status':'NUMERICAL_PILOT_ONLY_NO_STRICT_QUADRATURE_CERTIFICATE','point':[1,1,1,1,0],'pid':pid,'classes':len(classes),'triples':sum(classes.values()),'Gram_determinant_counts':dict(gram_dets),
      'geometry_classes':[{'invariants':list(k),'count':classes[k]} for k in keys],
      'sobol':{'points':N,'dimension':3,'seed':20261026,'scramble':True,'block_count':8},
      'mean_Green_formula':float(values.mean()),'consecutive_block_means':means,'block_min':min(means),'block_max':max(means),'largest_sample_contribution':float(values.max()),
      'CDF_quantiles_range':[float(uniforms.min()),float(uniforms.max())],'stable_retained_values_range':[float(z.min()),float(z.max())],'third_largest_values_range':[float(z[:,0].min()),float(z[:,0].max())],
      'relative_to128e_minus6':float(values.mean()/0.000128),'relative_to1e_minus4':float(values.mean()/0.0001),
      'scope':'Numerical Sobol estimate with fixed-table interpolation. Neither the mean nor block extrema are rigorous bounds on the integral oractualG. No full-domain claim.',
      'seconds':time.monotonic()-start,'python':sys.executable,**hashes}
    (HERE/'pilot_result.json').write_text(json.dumps(result,indent=2)+'\n')
    state.update(status='COMPLETE_EXIT0',seconds_this_checkpoint=time.monotonic()-start,exit_code=0)
    (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='geometry_classes'},indent=2),flush=True)

if __name__=='__main__':main()
