#!/usr/bin/env python3
"""Exact original-space continuous cube arrival and cone-relative tilt guards.
No tail MC; all evaluated receivers are registered finite rational probes.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter,defaultdict
from itertools import product
from math import prod,comb
import gzip,hashlib,json,time,sys
import numpy as np

BASE=Path(__file__).resolve().parent
PREFIX='critical_source_tail_arrival_20261007'
sys.set_int_max_str_digits(0)  # Bounded n64 exact RN fractions exceed 4300 digits.


def fs(x):return str(x)
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def source(n,A=None,ks=(1,2,4,8),weights=None):
    A=16*n if A is None else A
    weights=[F(1,2),F(1,6),F(1,6),F(1,6)] if weights is None else weights
    assert len(ks)==len(weights) and sum(weights)==1 and all(w>0 for w in weights)
    return dict(n=n,A=A,ks=list(ks),weights=weights,qs=[2*A*k+1 for k in ks])


def axis_count(x,R,k,A):
    low=max(-A*k,ceil(k*(x-R/2)));high=min(A*k,floor(k*(x+R/2)))
    return max(0,high-low+1)


def direct_tensor_response(data,x,R):
    masses=[F(prod(axis_count(xi,R,k,data['A']) for xi in x),q**data['n']) for k,q in zip(data['ks'],data['qs'])]
    return sum(w*m for w,m in zip(data['weights'],masses))/R**data['n']


def arrival_oracle(data,x):
    n=data['n'];events=defaultdict(Counter)
    counts=[[axis_count(xi,F(1),k,data['A']) for xi in x] for k in data['ks']]
    for l,k in enumerate(data['ks']):
        for i,xi in enumerate(x):
            low=max(-data['A']*k,ceil(k*(xi-1)))
            high=min(data['A']*k,floor(k*(xi+1)))
            for m in range(low,high+1):
                R=2*abs(xi-F(m,k))
                if 1<R<=2:events[R][(l,i)]+=1
    candidates=sorted({F(1),F(2),*events.keys()})
    nonzero=[prod(c for c in row if c) for row in counts]
    zeros=[sum(c==0 for c in row) for row in counts]
    qpowers=[q**n for q in data['qs']]
    best=None;bestR=None;trace=[];component_best=[None]*len(data['ks']);component_R=[None]*len(data['ks'])
    for R in candidates:
        for (l,i),increment in events.get(R,{}).items():
            old=counts[l][i];new=old+increment
            if old:nonzero[l]//=old
            else:zeros[l]-=1
            nonzero[l]*=new;counts[l][i]=new
        mass=[F(0) if zeros[l] else F(nonzero[l],qpowers[l]) for l in range(len(data['ks']))]
        component=[m/R**n for m in mass]
        full=sum(w*u for w,u in zip(data['weights'],component))
        if best is None or full>best:best,bestR=full,R
        for l,u in enumerate(component):
            if component_best[l] is None or u>component_best[l]:component_best[l],component_R[l]=u,R
        trace.append(dict(R=fs(R),response=fs(full),group_event_multiplicity=sum(events.get(R,{}).values()),
            distinct_updated_pairs=len(events.get(R,{}))))
    check_indices=sorted({0,len(candidates)-1,len(candidates)//2,candidates.index(bestR),min(1,len(candidates)-1),max(0,len(candidates)-2)})
    for i in check_indices:assert direct_tensor_response(data,x,candidates[i])==F(trace[i]['response'])
    return dict(winner_R=fs(bestR),M=fs(best),candidate_count=len(candidates),
        simultaneous_event_groups=sum(t['group_event_multiplicity']>1 for t in trace),
        total_arrival_events=sum(t['group_event_multiplicity'] for t in trace),
        largest_group_multiplicity=max(t['group_event_multiplicity'] for t in trace),
        component_winner_R=[fs(v) for v in component_R],component_max_response=[fs(v) for v in component_best],
        direct_tensor_selected_candidate_checks=len(check_indices),trace=trace)


def all_atom_catalog(data):
    masses=defaultdict(F)
    for k,q,w in zip(data['ks'],data['qs'],data['weights']):
        axis=[F(m,k) for m in range(-data['A']*k,data['A']*k+1)]
        for point in product(axis,repeat=data['n']):masses[point]+=w/F(q**data['n'])
    assert sum(masses.values())==1
    return masses


def scalar_all_atom_oracle(data,x):
    catalog=all_atom_catalog(data);groups=defaultdict(F)
    for y,mass in catalog.items():groups[2*max(abs(xi-yi) for xi,yi in zip(x,y))]+=mass
    candidates=sorted({F(1),F(2),*(R for R in groups if 1<R<=2)})
    mass=sum(m for R,m in groups.items() if R<=1);best=None;bestR=None
    for R in candidates:
        if R>1:mass+=groups.get(R,F(0))
        value=mass/R**data['n']
        if best is None or value>best:best,bestR=value,R
    return bestR,best,len(catalog)


def high_sets(yi,width,A):
    """Actual finite k=1 reference high count=2 at R*=3/2."""
    lo=yi-width/2;hi=yi+width/2;cuts={lo,hi}
    for m in range(max(-A,ceil(lo-F(3,4))),min(A,floor(hi+F(3,4)))+1):
        for z in [F(m)-F(3,4),F(m)+F(3,4)]:
            if lo<z<hi:cuts.add(z)
    cuts=sorted(cuts);segments=[];length=F(0)
    for left,right in zip(cuts,cuts[1:]):
        high=axis_count((left+right)/2,F(3,2),1,A)==2
        segments.append((left,right,high))
        if high:length+=right-left
    assert 0<=length<=width
    return segments,length/width


def status_law(probabilities):
    law=[F(1)]
    for p in probabilities:
        out=[F(0)]*(len(law)+1)
        for i,mass in enumerate(law):out[i]+=mass*(1-p);out[i+1]+=mass*p
        law=out
    assert sum(law)==1
    return law


def tilt_guard(data,y,width,face,sign,T,t):
    free=[i for i in range(data['n']) if i!=face]
    sets=[];ps=[]
    for i in free:
        segments,p=high_sets(y[i],width,data['A']);sets.append(segments);ps.append(p)
    Z=[1-p+p*T for p in ps];normalizer=prod(Z)
    prior=status_law(ps);tilted=status_law([p*T/z for p,z in zip(ps,Z)])
    eta=F(1,4);proposal=[eta*p+(1-eta)*q for p,q in zip(prior,tilted)]
    RN=[1/(eta+(1-eta)*T**h/normalizer) for h in range(len(prior))]
    assert all(0<w<=4 for w in RN)
    total=sum(q*w for q,w in zip(proposal,RN));assert total==1
    mean=sum(F(h)*q*w for h,(q,w) in enumerate(zip(proposal,RN)));assert mean==sum(ps)
    face_high=0 if face is None else int(axis_count(y[face]+sign*width/2,F(3,2),1,data['A'])==2)
    prior_tail=sum(p for h,p in enumerate(prior) if h+face_high>=t)
    weighted_tail=sum(q*w for h,(q,w) in enumerate(zip(proposal,RN)) if h+face_high>=t)
    assert prior_tail==weighted_tail
    second=sum(q*w*w for q,w in zip(proposal,RN));assert second<=4
    return dict(width=fs(width),face=face,sign=sign,face_high=face_high,
        source_y=[fs(v) for v in y],p_i=[fs(p) for p in ps],Z_i=[fs(z) for z in Z],
        high_set_intervals=[[[fs(a),fs(b),h] for a,b,h in segments] for segments in sets],
        prior_status_law=[fs(p) for p in prior],proposal_status_law=[fs(p) for p in proposal],RN_by_count=[fs(p) for p in RN],
        RN_mass=fs(total),RN_first_high_moment=fs(mean),prior_critical_tail=fs(prior_tail),weighted_critical_tail=fs(weighted_tail),RN_second_moment=fs(second),passed=True)


def primes(count):
    values=[];v=17
    while len(values)<count:
        if all(v%d for d in range(2,int(v**.5)+1)):values.append(v)
        v+=1
    return values


def source_point(data,rng):
    l=int(rng.choice(len(data['ks']),p=[float(w) for w in data['weights']]))
    k=data['ks'][l];A=data['A']
    m=rng.integers(-A*k,A*k+1,size=data['n'])
    return [F(int(v),k) for v in m],l


def summarize_round(plan,payload,profile,resumed=False):
    oracle_records=payload['oracle_records'];serial_data=payload['input']
    return dict(n=plan['n'],seed=plan['seed'],critical_t=plan['t_n'],critical_tau=payload['critical_tau'],binomial_tail=payload['binomial_tail'],
        input_sha256=hashlib.sha256(json.dumps(serial_data,sort_keys=True).encode()).hexdigest(),profile_path=str(profile),profile_sha256=sha(profile),
        oracle_probes=len(oracle_records),critical_exceeding_probes=sum(z['critical_exceedance'] for z in oracle_records),
        candidate_count_range=[min(z['candidate_count'] for z in oracle_records),max(z['candidate_count'] for z in oracle_records)],
        simultaneous_groups=sum(z['simultaneous_event_groups'] for z in oracle_records),
        largest_group_multiplicity=max(z['largest_group_multiplicity'] for z in oracle_records),
        direct_tensor_checks=sum(z['direct_tensor_selected_candidate_checks'] for z in oracle_records),
        component_vs_full_winner_different_probes=sum(any(z['winner_R']!=cr for cr in z['component_winner_R']) for z in oracle_records),
        tilt_configs=len(payload['tilt_records']),RN_checks_passed=all(z['passed'] for z in payload['tilt_records']),
        loaded_saved_profile_without_rerun=resumed)


def main():
    regpath=BASE/(PREFIX+'_registration.json');reg=json.loads(regpath.read_text())
    resultpath=BASE/(PREFIX+'_results.json');assert not resultpath.exists()
    start=time.perf_counter();scalar_records=[]
    for n in [1,2,3]:
        data=source(n,A=1,ks=(1,2,4),weights=[F(1,2),F(1,4),F(1,4)])
        for v in [F(0),F(1,2),F(1,4),F(-1,4),F(3,4),F(2,7)]:
            x=[v if i%2==0 else -v for i in range(n)]
            R,M,N=scalar_all_atom_oracle(data,x);fast=arrival_oracle(data,x)
            assert F(fast['winner_R'])==R and F(fast['M'])==M
            scalar_records.append(dict(n=n,x=[fs(z) for z in x],distinct_atoms=N,winner_R=fs(R),M=fs(M),passed=True))
    rounds=[]
    for plan in reg['critical_threshold']['round_parameters']:
        n,t,T=plan['n'],plan['t_n'],F(plan['tilt_odds_T']);data=source(n)
        tail=F(sum(comb(n,k) for k in range(t,n+1)),2**n)
        assert fs(tail)==plan['binomial_tail'] and tail<=F(1,n*n)
        assert t==n//2+1 or F(sum(comb(n,k) for k in range(t-1,n+1)),2**n)>F(1,n*n)
        tau=F(1,2)*F(2**t,data['qs'][0]**n)*F(2,3)**n
        profile=BASE/(PREFIX+f'_n{n}_exact_profiles.json.gz')
        if profile.exists():
            with gzip.open(profile,'rt',encoding='utf8') as f:saved=json.load(f)
            assert saved['critical_tau']==fs(tau) and saved['binomial_tail']==fs(tail)
            assert saved['input']['n']==n and saved['input']['ks']==data['ks'] and saved['input']['qs']==data['qs']
            rounds.append(summarize_round(plan,saved,profile,resumed=True))
            print(json.dumps(dict(n=n,status='loaded_saved_profile_without_rerun',profile_sha256=sha(profile))),flush=True)
            continue
        rng=np.random.default_rng(plan['seed']);oracle_records=[];ys=[];pv=primes(n)
        for j in range(8):
            y,label=source_point(data,rng);ys.append(y)
            width=[F(5,4),F(3,2),F(7,4),F(3,2)][j%4]
            face=j%n;sign=1 if j%2==0 else -1
            omega=[F(int(rng.integers(-7,8)),16) if j<4 else F(int(rng.integers(-p+1,p)),2*p) for p in pv]
            x=[yi+width*om for yi,om in zip(y,omega)];x[face]=y[face]+sign*width/2
            if j==7:
                for i in range(n):
                    if i==face:continue
                    intervals,p=high_sets(y[i],width,data['A']);high=[z for z in intervals if z[2]]
                    if high:
                        left,right,_=high[0];x[i]=left+F(pv[i]//2,pv[i])*(right-left)
            out=arrival_oracle(data,x)
            out.update(source_y=[fs(v) for v in y],receiver_x=[fs(v) for v in x],source_component=label,
                ray_width=fs(width),face=face,sign=sign,probe_type='dyadic_ties' if j<4 else 'prime_high_set' if j==7 else 'prime_free',
                critical_exceedance=F(out['M'])>tau)
            oracle_records.append(out)
        tilt=[]
        for j,(width,face,sign) in enumerate([(F(1),None,0),(F(5,4),0,1),(F(7,4),n-1,-1)]):
            tilt.append(tilt_guard(data,ys[j],width,face,sign,T,t))
        assert not profile.exists()
        serial_data=dict(n=n,A=data['A'],ks=data['ks'],qs=data['qs'],weights=[fs(w) for w in data['weights']],source_union_atoms=data['qs'][-1]**n)
        payload=dict(input=serial_data,critical_tau=fs(tau),binomial_tail=fs(tail),oracle_records=oracle_records,tilt_records=tilt)
        with gzip.open(profile,'wt',encoding='utf8') as f:json.dump(payload,f,separators=(',',':'))
        rounds.append(summarize_round(plan,payload,profile))
        print(json.dumps(rounds[-1]),flush=True)
    payload=dict(status='passed',registration_sha256=sha(regpath),script_sha256=sha(Path(__file__)),rounds=rounds,
        scalar_all_atom_checks=scalar_records,elapsed_seconds=time.perf_counter()-start,
        exact_arithmetic='Fraction on realized rational input, no tolerance grouping or floating response',
        no_MC_tail_estimate=True,no_source_moment_estimate=True,no_asymptotic_evidence=True,
        scope='Complete original finite atomic mixtures, continuous [1,2] exact winner and actual source-cone-relative RN cross-checks')
    resultpath.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(status='passed',scalar_checks=len(scalar_records),oracle_probes=24,tilt_configs=9,elapsed_seconds=payload['elapsed_seconds'])),flush=True)


if __name__=='__main__':main()
