#!/usr/bin/env python3
"""Two frozen inputs only; fresh source-pair/cone diagnostic, no optimization."""
from pathlib import Path
import datetime,hashlib,json,math,types
import numpy as np
HERE=Path(__file__).resolve().parent
module=types.ModuleType('cell_residual_functions');module.__file__=str(HERE/'cell_residual_shell_search.py')
exec(compile(Path(module.__file__).read_text(),module.__file__,'exec'),module.__dict__)
FILES=['validation_20261006T173636280680Z.json','cell_residual_validation_20261006T174331920727Z.json']
M=1024;K=257;SEED=202610076612


def geometry(record):
    atoms=np.array(record['atoms']);n=record['parameters']['n'];N=len(atoms)
    result=dict(n=n,N=N,full_source_linf_diameter=float(np.max(np.ptp(atoms,axis=0))))
    if record['parameters']['family']=='multires_packet_network':
        packets=[]
        for g in range((N+7)//8):
            a=atoms[g*8:min(N,(g+1)*8)]
            cell=np.floor(n*a).astype(np.int64)
            width=np.ptp(a,axis=0)
            packets.append(dict(packet=g,atoms=len(a),source_mass=sum(record['weights'][g*8:min(N,(g+1)*8)]),
                linf_diameter=float(width.max()),linf_diameter_over_cell_side=float(n*width.max()),
                distinct_origin0_cells=len(np.unique(cell,axis=0)),coordinates_crossing_cell_boundaries=int(np.count_nonzero(np.ptp(cell,axis=0))),
                min_coordinate_span_over_cell_side=float(n*width.min()),mean_coordinate_span_over_cell_side=float(n*width.mean())))
        result['packet_geometry']=packets
        result['interpretation']='These are small 1/n-scale packets, potentially straddling several origin0 fine cells; fixed-cell cooperation does not certify robustness under grid shifts, or exclude geometry close to a point packet. No shifted-grid search is performed.'
    return result


def probe(record,index,is_residual):
    atoms=np.array(record['atoms']);w=np.array(record['weights']);n=record['parameters']['n'];N=len(w)
    _,cells=np.unique(np.floor(n*atoms).astype(np.int64),axis=0,return_inverse=True)
    samples=[]
    for side in range(2):
        rng=np.random.default_rng(SEED+index*104729+side*1000003)
        source=rng.choice(N,M,p=w);omega=rng.uniform(-.5,.5,(M,n));face=rng.integers(0,2*n,M)
        omega[np.arange(M),face//2]=np.where(face%2,.5,-.5)
        samples.append((source,omega,face,atoms[source,None,:]-atoms[None,:,:]))
    ts=np.linspace(0,n*math.log(2),K)
    names=['total','mutual','nonmutual','same_face_same_sign','mutual_same_face_same_sign','mutual_other_faces']
    trajectories=np.zeros((M,len(names),K));thin_violations=0;max_thin_excess=0.
    ll=record['selected_log_lambda'];rows=np.arange(M)
    for ki,t in enumerate(ts):
        r=math.exp(t/n);events=[];Qs=[]
        for source,omega,face,diff in samples:
            q=2*np.max(np.abs(diff/r+omega[:,None,:]),axis=2)
            assert np.all(q[rows,source]==1.)
            qs,lq,mass,_,_,ix,lu=module.engine.window_maximal(q,w,n,t)
            delay=n*lq[rows,ix];event=(lu>ll+math.log(2))&(lu<=ll+math.log(4))&(delay>=0)&(delay<=1)
            if is_residual:
                cq,cap=module.captured_cell_maxima(q,w,cells,t,n)
                assert np.array_equal(cq,qs)
                event&=cap[rows,ix]<=.5*mass[rows,ix]
            events.append(event);Qs.append(qs[rows,ix])
        sA,omegaA,faceA,_=samples[0];sB,omegaB,faceB,_=samples[1]
        delta=atoms[sB]-atoms[sA]
        crossA=2*np.max(np.abs(delta/r-omegaA),axis=1)
        crossB=2*np.max(np.abs(-delta/r-omegaB),axis=1)
        both=events[0]&events[1];mutual=(crossA<=Qs[0])&(crossB<=Qs[1]);sameface=faceA==faceB
        masks=[both,both&mutual,both&~mutual,both&sameface,both&mutual&sameface,both&mutual&~sameface]
        for j,mask in enumerate(masks):trajectories[:,j,ki]=mask
        relevant=both&mutual&sameface
        if np.any(relevant):
            facecoord=faceA//2
            bound=r/2*math.expm1(1/n)
            excess=np.abs(delta[rows,facecoord])-bound
            thin_violations+=int(np.count_nonzero(relevant&(excess>1e-10)))
            max_thin_excess=max(max_thin_excess,float(np.max(excess[relevant])))
    values=np.trapezoid(trajectories,ts,axis=2)
    components={}
    totals=values[:,0];mean_total=float(totals.mean())
    # All paired trajectories are iid across row, despite CRN across t.
    # Ordinary SE is diagnostic; finite-node uniform Hoeffding remains wide.
    epsilon=math.sqrt(math.log(2*2*len(names)*K/.05)/(2*M))
    for j,name in enumerate(names):
        y=values[:,j];mean=float(y.mean());fraction=mean/mean_total if mean_total else None
        ratio_se=(float(np.std(y-fraction*totals,ddof=1)/math.sqrt(M)/mean_total) if mean_total else None)
        components[name]=dict(energy_finite_trapezoid_mean=mean,ordinary_MC_standard_error=float(y.std(ddof=1)/math.sqrt(M)),
            fraction_of_total=fraction,delta_method_fraction_standard_error=ratio_se,
            event_sample_node_count=int(trajectories[:,j].sum()),mean_at_nodes=trajectories[:,j].mean(axis=0).tolist())
    assert np.allclose(values[:,1]+values[:,2],values[:,0])
    assert np.allclose(values[:,4]+values[:,5],values[:,1])
    return dict(candidate_id=record['candidate_id'],input_sha256=record['input_sha256'],log_lambda=ll,
        model='fixed-cell eta=.5 cooperative residual' if is_residual else 'full hard relaxation',
        n=n,N=N,pairs=M,fine_nodes=K,seed_base=SEED+index*104729,
        components=components,t=ts.tolist(),finite_node_uniform_Hoeffding_epsilon=epsilon,
        thin_same_face_sign_violation_count=thin_violations,max_observed_thin_bound_excess=max_thin_excess,
        geometry=geometry(record))


def main():
    out=[]
    for i,file in enumerate(FILES):
        doc=json.loads((HERE/file).read_text())
        r=max(doc['records'],key=lambda r:r['quantities']['Gamma_v1']['energy_cross_fine'])
        out.append(probe(r,i,i==1))
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    dest=HERE/f'frozen_pair_diagnostic_{stamp}.json'
    data=dict(records=out,source_files=FILES,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        protocol='Exactly two previously frozen input/lambda pairs; new independent pair source labels and cone directions; no input/lambda optimization. All same-face labels encode both coordinate and sign. Actual geometric delay and complete hardband are evaluated with floating candidate enumeration.',
        geometric_mechanism='For mutual capture with both delays<=1, same face/sign forces |delta_face|<=r/2*(exp(1/n)-1); other-face and nonmutual cases get no such thin-face consequence. Testing this identity is not an energy bound.',
        confidence='Ordinary SE and delta-method fraction SE are diagnostics under iid pair-row sampling; ratio SE is unreliable for very rare components. Hoeffding epsilon gives a global 95% finite-node joint band over 2 inputs x 6 bounded categories x 257 nodes, assuming floating classification is correct. No continuum error certificate; no actual geom gates.',
        limits=['One fresh pair probe per frozen input; no independent implementation audit','Small source packets may straddle origin0 grid boundaries; no shifted-grid robustness claim','Mutual/nonmutual proportions describe these finite inputs only; no asymptotic inference'])
    with dest.open('x') as f:json.dump(data,f,ensure_ascii=False,indent=2,allow_nan=False)
    for r in out:
        print(r['model'],r['candidate_id'],{k:(v['energy_finite_trapezoid_mean'],v['fraction_of_total'],v['ordinary_MC_standard_error']) for k,v in r['components'].items()},flush=True)
    print('OUTPUT',dest,flush=True)


if __name__=='__main__':main()
