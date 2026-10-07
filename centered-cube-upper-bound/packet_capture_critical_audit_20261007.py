#!/usr/bin/env python3
"""Exact capture postprocessing only; no arrival oracle import or execution."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import gzip,hashlib,json,time,sys
sys.set_int_max_str_digits(0)
BASE=Path(__file__).resolve().parent
PREFIX='packet_capture_critical_audit_20261007'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)


def count_axis(x,R,k,A):
    left=max(-A*k,ceil(k*(x-R/2)));right=min(A*k,floor(k*(x+R/2)))
    return max(0,right-left+1)


def main():
    start=time.perf_counter()
    regpath=BASE/(PREFIX+'_registration.json');reg=json.loads(regpath.read_text())
    out=BASE/(PREFIX+'_results.json');assert not out.exists()
    original_result=BASE/'critical_source_tail_arrival_20261007_results.json'
    saved=json.loads(original_result.read_text());records=[];rounds=[];input_hashes={original_result.name:sha(original_result)}
    for row in saved['rounds']:
        profile=Path(row['profile_path']);assert sha(profile)==row['profile_sha256']
        input_hashes[profile.name]=sha(profile)
        with gzip.open(profile,'rt') as f:data=json.load(f)
        inp=data['input'];n=inp['n'];theta=F(reg['theta_by_n'][str(n)])
        assert theta*theta==F(1,n)
        assert hashlib.sha256(json.dumps(inp,sort_keys=True).encode()).hexdigest()==row['input_sha256']
        weights=[F(w) for w in inp['weights']];assert weights==[F(1,2),F(1,6),F(1,6),F(1,6)] and sum(weights)==1
        bounds=[F(2*k+1,q)**n for k,q in zip(inp['ks'],inp['qs'])]
        for bound in bounds:assert bound<F(1,100)
        indices=[]
        for point,record in enumerate(data['oracle_records']):
            x=[F(v) for v in record['receiver_x']];R=F(record['winner_R']);M=F(record['M']);tau=F(data['critical_tau'])
            assert 1<=R<=2
            E=M>tau;assert E==record['critical_exceedance']
            counts=[[count_axis(xi,R,k,inp['A']) for xi in x] for k in inp['ks']]
            fullness=[F(prod(c),q**n) for c,q in zip(counts,inp['qs'])]
            mi=[w*f for w,f in zip(weights,fullness)];m=sum(mi)
            assert m==M*R**n and m>0
            shares=[v/m for v in mi];assert sum(shares)==1
            assert all(f<=bound for f,bound in zip(fullness,bounds))
            common=dict(n=n,source_seed=row['seed'],point_index=point,probe_type=record['probe_type'],
                receiver_x=record['receiver_x'],saved_winner_R=record['winner_R'],saved_full_M=record['M'],saved_tau=data['critical_tau'],
                saved_E=E,theta=str(theta),full_m=str(m),packet_Mi=[str(w) for w in weights],packet_mi=[str(v) for v in mi],
                packet_capture_counts=counts,capture_share_mi_over_m=[str(v) for v in shares],fullness_mi_over_Mi=[str(v) for v in fullness])
            for estr in reg['etas']:
                eta=F(estr);dom=[v>=theta for v in shares];full=[v>=eta for v in fullness]
                both=[d and f for d,f in zip(dom,full)]
                parts=dict(paid=F(0),dominance_only_failure=F(0),fullness_only_failure=F(0),both_failure=F(0))
                for share,d,f in zip(shares,dom,full):
                    key='paid' if d and f else 'dominance_only_failure' if not d and f else 'fullness_only_failure' if d and not f else 'both_failure'
                    parts[key]+=share
                assert sum(parts.values())==1
                qraw=parts['paid'];qdom=qraw if E else F(0)
                assert qraw==0 and not any(full)  # Uniform bound proves it, not sample averages.
                indices.append(len(records))
                records.append(dict(**common,eta=estr,dominance_pass=dom,fullness_pass=full,qualified_packet=both,
                    packet_count_dominance_pass=sum(dom),packet_count_fullness_pass=sum(full),packet_count_qualified=sum(both),
                    raw_qualified_capture_fraction=str(qraw),q_dom_on_E=str(qdom),
                    captured_fraction_partition={k:str(v) for k,v in parts.items()},
                    dominance_failure_capture_fraction=str(parts['dominance_only_failure']+parts['both_failure']),
                    fullness_failure_capture_fraction=str(parts['fullness_only_failure']+parts['both_failure'])))
        summary=[]
        for eta in reg['etas']:
            at=[records[i] for i in indices if records[i]['eta']==eta];ev=[v for v in at if v['saved_E']]
            summary.append(dict(eta=eta,points=len(at),saved_E_points=len(ev),fullness_fail_incidences=4*len(at),
                saved_E_fullness_fail_incidences=4*len(ev),dominance_pass_incidences=sum(v['packet_count_dominance_pass'] for v in at),
                saved_E_dominance_pass_incidences=sum(v['packet_count_dominance_pass'] for v in ev),
                q_dom_all_exact_zero=True,
                saved_E_dominance_failure_share_range=[str(min(F(v['dominance_failure_capture_fraction']) for v in ev)),str(max(F(v['dominance_failure_capture_fraction']) for v in ev))],
                fullness_failure_share_exact='1',point_statistics_are_not_Lebesgue_coverage=True))
        rounds.append(dict(n=n,theta=str(theta),packet_fullness_uniform_upper_bounds=[str(b) for b in bounds],
            uniform_all_receivers_fullness_failure_for_all_registered_eta=True,record_indices=indices,summary=summary))
        print(json.dumps(dict(n=n,saved_E_points=sum(v['critical_exceedance'] for v in data['oracle_records']),records=len(indices),fullness_all_failed=True)),flush=True)
    payload=dict(status='complete',registration_sha256=sha(regpath),script_sha256=sha(Path(__file__)),
        input_hashes_sha256=input_hashes,rounds=rounds,records=records,elapsed_seconds=time.perf_counter()-start,
        oracle_rerun=False,new_source_or_receiver_samples=False,full_denominator_preserved=True,
        no_probability_or_I_dom_estimate=True,no_general_coverage_or_square_theorem_verification=True,
        conclusion='These four global layer packets fail fullness uniformly for R in [1,2] and eta>=1/100. Other lawful pre-fixed spatial allocations are not excluded.')
    out.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(status='complete',rows=len(records),packet_eta_incidences=4*len(records),elapsed_seconds=payload['elapsed_seconds'])),flush=True)


if __name__=='__main__':main()
