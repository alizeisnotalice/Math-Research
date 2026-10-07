#!/usr/bin/env python3
"""Frozen full-source packet partition; no new source/threshold search."""
from pathlib import Path
import datetime,hashlib,json,math,time,types
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'adversarial_shell_pressure_20261007/cell_residual_validation_20261006T174331920727Z.json'
engine=types.ModuleType('winner_engine')
exec(compile((HERE/'winner_engine_snapshot.py').read_text(),str(HERE/'winner_engine_snapshot.py'),'exec'),engine.__dict__)
ETAS=[.5,.25,.125]
ROUNDS=[(1,512,129,202610078012),(2,1024,257,202610079012),(3,2048,513,202610080012)]


def exact_mass_certificate(weights,packets):
    ratios=[float(x).as_integer_ratio() for x in weights]
    exponent=max(d.bit_length()-1 for _,d in ratios)
    integers=[a*(2**exponent//b) for a,b in ratios]
    subsets=[]
    for packet in range(int(packets.max())+1):
        sums={0}
        for j in np.flatnonzero(packets==packet):sums|={v+integers[int(j)] for v in tuple(sums)}
        subsets.append(sums)
    intersection=set.intersection(*subsets)-{0}
    return dict(binary_denominator_exponent=exponent,
       exact_atom_mass_numerators=[str(v) for v in integers],
       exact_source_mass_numerator=str(sum(integers)),
       packet_subset_sum_counts=[len(s) for s in subsets],
       packet_subset_sum_sha256=[hashlib.sha256('\n'.join(map(str,sorted(s))).encode()).hexdigest() for s in subsets],
       positive_all_packet_subset_mass_intersection_numerators=[str(v) for v in sorted(intersection)],
       quarter_empty_certificate=len(subsets)==4 and not intersection,
       scope='Exact binary-rational interpretation of the saved positive atom weights only. For four packets eta=1/4 forces every captured packet mass equal; empty positive subset-mass intersection rules that out for ANY geometric capture subset. This does not certify ideal analytic pre-rounding generator weights.'),integers


def grouped_maxima(q,weights,groups,t,n):
    M,N=q.shape;rows=np.arange(M)
    order=np.argsort(q,axis=1,kind='stable')
    original=np.take_along_axis(q,order,axis=1);r=math.exp(t/n)
    qs=np.concatenate([original,np.broadcast_to([1/r,2/r],(M,2))],axis=1)
    ww=np.concatenate([weights[order],np.zeros((M,2))],axis=1)
    labels=[np.concatenate([g[order],np.zeros((M,2),dtype=int)],axis=1) for g in groups]
    perm=np.argsort(qs,axis=1,kind='stable')
    qs=np.take_along_axis(qs,perm,axis=1);ww=np.take_along_axis(ww,perm,axis=1)
    labels=[np.take_along_axis(g,perm,axis=1) for g in labels]
    last=np.concatenate([qs[:,:-1]!=qs[:,1:],np.ones((M,1),dtype=bool)],axis=1)
    maxima=[];active_counts=[]
    for lab,g in zip(labels,groups):
        state=np.zeros((M,int(g.max())+1));peak=np.zeros(M);active=np.zeros(M)
        values=np.empty_like(qs);counts=np.empty_like(qs)
        for j in range(N+2):
            active+=(state[rows,lab[:,j]]==0)&(ww[:,j]>0)
            state[rows,lab[:,j]]+=ww[:,j]
            peak=np.maximum(peak,state[rows,lab[:,j]])
            values[:,j]=peak;counts[:,j]=active
        maxima.append(np.minimum.accumulate(np.where(last,values,np.inf)[:,::-1],axis=1)[:,::-1])
        active_counts.append(np.minimum.accumulate(np.where(last,counts,np.inf)[:,::-1],axis=1)[:,::-1])
    return qs,maxima,active_counts


def exact_winner_packet_gate(q,Q,weights_int,packets,denominator,rows_to_check):
    answers={}
    for row in rows_to_check:
        masses=[0]*(int(packets.max())+1)
        for j in np.flatnonzero(q[row]<=Q[row]):masses[int(packets[j])]+=weights_int[int(j)]
        answers[int(row)]=denominator*max(masses)<=sum(masses)
    return answers


def validate():
    q=np.ones((1,4));w=np.full(4,.25)
    qs,maxima,active=grouped_maxima(q,w,[np.array([0,0,1,2]),np.array([0,1,2,3])],0.,8)
    assert np.all(maxima[0]==.5) and np.all(maxima[1]==.25)
    assert np.all(active[1]==4)
    cert,_=exact_mass_certificate(np.full(32,1/32),np.arange(32)//8)
    assert not cert['quarter_empty_certificate'] and len(cert['positive_all_packet_subset_mass_intersection_numerators'])==8
    # Strict paid / closed residual equality and energy reconstruction.
    assert not (.25>.25*1.) and .25<=.25*1.
    P_A=np.array([.2,.1]);R_A=np.array([.4,.3]);P_B=np.array([.3,.2]);R_B=np.array([.1,.2])
    assert np.allclose((P_A+R_A)*(P_B+R_B),P_A*P_B+R_A*R_B+P_A*R_B+R_A*P_B)
    return ['completed distance groups aggregate cell and packet masses','uniform four-packet quarter equality remains possible','strict paid/closed residual boundary','all squared cross terms reconstruct total']


def profile_stats(raw,lo,hi,ts,epsilon):
    a,b=raw;mean=(a+b)/2;lower_mean=lo.mean(axis=0);upper_mean=hi.mean(axis=0)
    cross=a*b;fine=float(np.trapezoid(cross,ts));coarse=float(np.trapezoid(cross[::2],ts[::2]))
    lower=np.maximum(0,lower_mean-epsilon);upper=np.minimum(1,upper_mean+epsilon)
    return dict(rep_A=a.tolist(),rep_B=b.tolist(),mean=mean.tolist(),
       screened_lower_mean=lower_mean.tolist(),screened_upper_mean=upper_mean.tolist(),
       finite_node_95_lower_profile=lower.tolist(),finite_node_95_upper_profile=upper.tolist(),
       energy_cross_fine=fine,energy_cross_coarse=coarse,nested_grid_difference=abs(fine-coarse),
       energy_empirical_square=float(np.trapezoid(mean**2,ts)),
       energy_finite_node_95_lower=float(np.trapezoid(lower**2,ts)),
       energy_finite_node_95_upper=float(np.trapezoid(upper**2,ts)),
       integral_mean=float(np.trapezoid(mean,ts)),positive_fine_nodes=int(np.count_nonzero(mean)))


def evaluate(record,round_,M,K,seed):
    atoms=np.array(record['atoms']);w=np.array(record['weights']);n=record['parameters']['n'];N=len(w)
    digest=hashlib.sha256(atoms.astype('<f8').tobytes()+w.astype('<f8').tobytes()).hexdigest()
    assert digest==record['input_sha256']
    packets=np.arange(N)//8;P=int(packets.max())+1
    cells,cell_ids=np.unique(np.floor(n*atoms).astype(np.int64),axis=0,return_inverse=True)
    certificate,integers=exact_mass_certificate(w,packets)
    impossible=[eta<1/P or (P==4 and eta==.25 and certificate['quarter_empty_certificate']) for eta in ETAS]
    ts=np.linspace(0,n*math.log(2),K);rows=np.arange(M);ll=record['selected_log_lambda']
    # replicate, base (full hard / fixed-cell residual), quantity (Gamma1/Theta),
    # part: total, paid/residual at eta=.5,.25,.125; no renormalization.
    raw=np.zeros((2,2,2,7,K));lower=np.zeros_like(raw);upper=np.zeros_like(raw)
    ambiguity=np.zeros((2,2,7),dtype=np.int64)
    exact_gate_checks=np.zeros(3,dtype=np.int64);forced_empty_float_disagreements=np.zeros(3,dtype=np.int64)
    diagnostics={k:np.zeros(K) for k in ['winner_mean_packet_max_share','winner_mean_active_packet_count','source_membership','winner_mean_mass','winner_mean_atom_count']}
    for rep in range(2):
        rng=np.random.default_rng(seed+rep*1000003)
        source=rng.choice(N,M,p=w/w.sum());omega=rng.uniform(-.5,.5,(M,n));face=rng.integers(0,2*n,M)
        omega[rows,face//2]=np.where(face%2,.5,-.5)
        diff=atoms[source,None,:]-atoms[None,:,:]
        for ki,t in enumerate(ts):
            r=math.exp(t/n);q=2*np.max(np.abs(diff/r+omega[:,None,:]),axis=2)
            assert np.all(q[rows,source]==1.)
            qs,lq,mass,valid,response,ix,lu=engine.window_maximal(q,w,n,t)
            cq,maxima,active=grouped_maxima(q,w,[cell_ids,packets],t,n)
            assert np.array_equal(cq,qs)
            max_cell,max_packet=maxima;packet_active=active[1]
            tol=(n+1)*1e-10;near=valid&(response>=lu[:,None]-tol)
            delay=n*lq;d=delay[rows,ix];Q=qs[rows,ix];m=mass[rows,ix];peak=max_packet[rows,ix]
            assert np.all(peak/m>=1/P-1e-12)
            diagnostics['winner_mean_packet_max_share'][ki]+=float(np.mean(peak/m))/2
            diagnostics['winner_mean_active_packet_count'][ki]+=float(packet_active[rows,ix].mean())/2
            diagnostics['source_membership'][ki]+=float(np.mean(Q>=1))/2
            diagnostics['winner_mean_mass'][ki]+=float(m.mean())/2
            diagnostics['winner_mean_atom_count'][ki]+=float((q<=Q[:,None]).sum(axis=1).mean())/2
            band=(lu>ll+math.log(2))&(lu<=ll+math.log(4))
            blo=(lu-tol>ll+math.log(2))&(lu+tol<=ll+math.log(4))
            bhi=(lu+tol>ll+math.log(2))&(lu-tol<=ll+math.log(4))
            capture_lo=qs>=1.;capture_hi=qs>=1.-tol
            quantity_lo=[capture_lo&(delay<=1.-tol),capture_lo]
            quantity_hi=[capture_hi&(delay<=1.+tol),capture_hi]
            quantity_raw=[(d>=0)&(d<=1),d>=0]
            values=[np.ones_like(qs),np.exp(-np.maximum(delay,0))]
            values_raw=[np.ones(M),np.exp(-np.maximum(d,0))]
            cell_margin=.5*mass-max_cell
            cell_lo=cell_margin>tol*np.maximum(mass,1e-300)
            cell_hi=cell_margin>=-tol*np.maximum(mass,1e-300)
            cell_raw=cell_margin[rows,ix]>=0
            base_lo=[np.ones_like(valid),cell_lo];base_hi=[np.ones_like(valid),cell_hi];base_raw=[np.ones(M,dtype=bool),cell_raw]
            gates_raw=[np.ones(M,dtype=bool)];gates_lo=[np.ones_like(valid)];gates_hi=[np.ones_like(valid)]
            for ei,eta in enumerate(ETAS):
                margin=eta*mass-max_packet
                residual=margin[rows,ix]>=0
                if impossible[ei]:
                    forced_empty_float_disagreements[ei]+=np.count_nonzero(residual)
                    residual=np.zeros(M,dtype=bool)
                    residual_lo=np.zeros_like(valid);residual_hi=np.zeros_like(valid)
                    paid_lo=np.ones_like(valid);paid_hi=np.ones_like(valid)
                else:
                    close=np.flatnonzero(np.abs(margin[rows,ix])<=tol*m)
                    exact_gate_checks[ei]+=len(close)
                    answers=exact_winner_packet_gate(q,Q,integers,packets,int(round(1/eta)),close)
                    for row,answer in answers.items():residual[row]=answer
                    residual_lo=margin>tol*np.maximum(mass,1e-300)
                    residual_hi=margin>=-tol*np.maximum(mass,1e-300)
                    paid_lo=margin<-tol*np.maximum(mass,1e-300)
                    paid_hi=margin<tol*np.maximum(mass,1e-300)
                gates_raw.extend([~residual,residual]);gates_lo.extend([paid_lo,residual_lo]);gates_hi.extend([paid_hi,residual_hi])
            for bi in range(2):
                for qi in range(2):
                    for part in range(7):
                        raw[rep,bi,qi,part,ki]=np.mean(values_raw[qi]*quantity_raw[qi]*band*base_raw[bi]*gates_raw[part])
                        yl=values[qi]*(quantity_lo[qi]&blo[:,None]&base_lo[bi]&gates_lo[part])
                        yu=values[qi]*(quantity_hi[qi]&bhi[:,None]&base_hi[bi]&gates_hi[part])
                        a=np.min(np.where(near,yl,np.inf),axis=1);b=np.max(np.where(near,yu,-np.inf),axis=1)
                        lower[rep,bi,qi,part,ki]=a.mean();upper[rep,bi,qi,part,ki]=b.mean()
                        ambiguity[bi,qi,part]+=np.count_nonzero(b-a>1e-14)
    # 15 frozen input-round records, 2 bases, 2 quantities, 7 parts,
    # 2 screening RVs, 2 tails: union factor 15*2*2*7*2*2*K.
    epsilon=math.sqrt(math.log(15*2*2*7*2*2*K/.05)/(4*M))
    bases={}
    reconstruction_max=0.
    for bi,bname in enumerate(['full_hard','fixed_cell_eta_half_residual']):
        totals={qname:profile_stats(raw[:,bi,qi,0],lower[:,bi,qi,0],upper[:,bi,qi,0],ts,epsilon) for qi,qname in enumerate(['Gamma_v1','Theta_hard'])}
        splits=[]
        for ei,eta in enumerate(ETAS):
            paid_index,res_index=1+2*ei,2+2*ei;parts={}
            for qi,qname in enumerate(['Gamma_v1','Theta_hard']):
                A=raw[0,bi,qi];B=raw[1,bi,qi]
                error=float(np.max(np.abs(A[0]-A[paid_index]-A[res_index])))
                error=max(error,float(np.max(np.abs(B[0]-B[paid_index]-B[res_index]))))
                reconstruction_max=max(reconstruction_max,error)
                cross=float(np.trapezoid(A[paid_index]*B[res_index]+A[res_index]*B[paid_index],ts))
                paid=profile_stats(raw[:,bi,qi,paid_index],lower[:,bi,qi,paid_index],upper[:,bi,qi,paid_index],ts,epsilon)
                residual=profile_stats(raw[:,bi,qi,res_index],lower[:,bi,qi,res_index],upper[:,bi,qi,res_index],ts,epsilon)
                recon=paid['energy_cross_fine']+residual['energy_cross_fine']+cross
                assert math.isclose(recon,totals[qname]['energy_cross_fine'],abs_tol=1e-12)
                # Joint node bands also bound 2*integral phi_paid*phi_residual.
                cross_lo=float(np.trapezoid(2*np.array(paid['finite_node_95_lower_profile'])*np.array(residual['finite_node_95_lower_profile']),ts))
                cross_hi=float(np.trapezoid(2*np.array(paid['finite_node_95_upper_profile'])*np.array(residual['finite_node_95_upper_profile']),ts))
                parts[qname]=dict(paid=paid,residual=residual,paid_residual_cross_energy=cross,
                    cross_energy_finite_node_95_lower=cross_lo,cross_energy_finite_node_95_upper=cross_hi,
                    reconstructed_total_energy=recon,profile_reconstruction_max_error=error,
                    ambiguity_sample_nodes=dict(total=int(ambiguity[bi,qi,0]),paid=int(ambiguity[bi,qi,paid_index]),residual=int(ambiguity[bi,qi,res_index])))
            splits.append(dict(eta=eta,analytic_empty_for_stored_model=impossible[ei],quantities=parts))
        bases[bname]=dict(total=totals,packet_splits=splits)
    assert reconstruction_max<1e-12
    return dict(candidate_id=record['candidate_id'],round=round_,n=n,N=N,packet_count=P,input_sha256=digest,
       atoms=record['atoms'],weights=record['weights'],complete_source_mass=float(w.sum()),
       sampling='sample original complete mu/W; no branch/source renormalization; w retained unchanged in all maximal responses',
       frozen_log_lambda=ll,packet_assignment=packets.tolist(),fixed_cell_ids=cell_ids.tolist(),
       sample_seed=seed,samples_per_replicate=M,independent_replicates=2,fine_nodes=K,coarse_nodes=(K+1)//2,t=ts.tolist(),
       finite_node_joint_Hoeffding_epsilon=epsilon,bases=bases,
       exact_dyadic_mass_certificate=certificate,pigeonhole_lower_share=1/P,
       forced_empty_float_comparison_disagreements=forced_empty_float_disagreements.tolist(),
       exact_winner_gate_near_threshold_checks=exact_gate_checks.tolist(),
       diagnostics={k:v.tolist() for k,v in diagnostics.items()},
       profile_reconstruction_max_error=reconstruction_max)


def save(path,data):
    with path.open('x') as f:json.dump(data,f,ensure_ascii=False,indent=2,allow_nan=False)


def main():
    start=time.perf_counter();stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    original=json.loads(SOURCE.read_text());records=original['records'];assert len(records)==5
    plan=dict(source_file=str(SOURCE),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
       rounds=ROUNDS,etas=ETAS,packet_assignment='fixed floor(original atom_index/8), no new grouping, no source/weight/lambda/winner optimization',
       bases='Primary: existing fixed-cell eta=.5 residual. Also full original hardband base to expose the complete source partition.',
       quantities='Gamma_v1=E[hardband,base,0<=beta-t<=1]; Theta_hard=E[hardband,base,beta>=t]*exp(t-beta)',
       decomposition='Profiles total=paid+residual at every node; energies total=paid^2+residual^2+2 paid*residual. Independent replicate energy uses PA*PB+RA*RB+PA*RB+RA*PB.',
       confidence='Unbiased cross-replicate products for fixed floating-node quantities; empirical square is upward biased. Hoeffding joint bands conditional on frozen input/lambda and heuristic arithmetic screening across 15 input-round records, both bases, 2 quantities, 7 parts, lower/upper RVs, both tails. Trapezoid bounds only finite-node targets, not continuum error. Never interpret a Monte Carlo zero as analytic zero.',
       arithmetic='True full continuous [1,2] candidate enumeration using frozen complete mu, endpoints and complete floating distance ties; source q=1 preserved. Near response/band/branch threshold tolerance (n+1)*1e-10 is not exact-real geometry certification. Near packet-share winner decisions use exact dyadic subset masses for stored weights; branch is never used to reselect the winner.',
       pigeonhole='max captured packet share >=1/packet_count. N32 has 4 packets: eta=.125 impossible; eta=.25 forces exact 4-way equality, checked by positive subset-mass intersection in stored dyadic model. N128 has 16 packets: do not extrapolate 4-packet emptiness.',
       jacobian='normalized cube cone boundary [-.5,.5]^n gives dx=r^n dt dnu; source-cone expectations count complete weighted labels, not unweighted Lebesgue samples',
       L02_skill='/Users/zhengzhihao/.codex/skills/math-l02-adversarial-input-search/SKILL.md',
       L02_references_read=['provenance.md','method.md','cube-interface.md','verification-log-20261006.md','audit-sol_hn.md','examples/positive.md','examples/negative.md'],
       skill_use='Reconstruct frozen feasible input and independently recompute actual response; no external optimization theorem or toy weak-ratio theorem is invoked.',
       validation=validate(),limits=['genuine hard continuous-window partition, not actual FIRST/GOOD/CP/GP/early-jump gates','No original arbitrary finite J or diffuse L1 conclusion','No asymptotic fitting, no new source search','Exact mass certificates concern saved dyadic weights; ordinary numeric zeros remain numeric zeros'],
       engine_sha256=hashlib.sha256((HERE/'winner_engine_snapshot.py').read_bytes()).hexdigest())
    save(HERE/f'plan_{stamp}.json',plan)
    allresults=[]
    for round_,M,K,seedbase in ROUNDS:
        current=[]
        for index,record in enumerate(records):
            r=evaluate(record,round_,M,K,seedbase+index*104729);current.append(r)
            base=r['bases']['fixed_cell_eta_half_residual'];split=base['packet_splits'][0]['quantities']['Gamma_v1']
            print('ROUND',round_,'n',r['n'],'N',r['N'],'Gamma total/paid/res/cross',base['total']['Gamma_v1']['energy_cross_fine'],split['paid']['energy_cross_fine'],split['residual']['energy_cross_fine'],split['paid_residual_cross_energy'],flush=True)
        out=HERE/f'round_{round_}_{stamp}.json';save(out,dict(records=current,runtime_seconds=time.perf_counter()-start));allresults.extend(current)
    final_summary=[]
    for r in allresults:
        base=r['bases']['fixed_cell_eta_half_residual'];s=[]
        for split in base['packet_splits']:
            s.append(dict(eta=split['eta'],analytic_empty_for_stored_model=split['analytic_empty_for_stored_model'],
             quantities={q:{k:v[k] for k in ['paid_residual_cross_energy','reconstructed_total_energy']}|{'paid_energy':v['paid']['energy_cross_fine'],'residual_energy':v['residual']['energy_cross_fine'],'residual_node_95_upper':v['residual']['energy_finite_node_95_upper'],'residual_positive_nodes':v['residual']['positive_fine_nodes']} for q,v in split['quantities'].items()}))
        final_summary.append(dict(candidate_id=r['candidate_id'],round=r['round'],n=r['n'],N=r['N'],packet_count=r['packet_count'],
            primary_base_totals={q:v['energy_cross_fine'] for q,v in base['total'].items()},splits=s,
            quarter_mass_certificate_empty=r['exact_dyadic_mass_certificate']['quarter_empty_certificate']))
    dest=HERE/f'summary_{stamp}.json';save(dest,dict(records=final_summary,runtime_seconds=time.perf_counter()-start,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        conclusion='Partition audit of five frozen atomic inputs only; no source/threshold search and no general residual estimate. Exact mass/pigeonhole empty branches are separately labeled from sampled zero branches.'))
    print('OUTPUT',dest,'RUNTIME',time.perf_counter()-start,flush=True)


if __name__=='__main__':main()
