#!/usr/bin/env python3
"""Continue metadata join only. Never rerun completed455 endpoint checks.

The first audit expected158 kind='evaluated';the frozen158 producer writes
'evaluated_single_source'.This fixes that vocabulary at the metadata boundary.
"""
from pathlib import Path
from collections import Counter
import argparse,hashlib,json,os,time,resource
import audit_binding as A
H=Path(__file__).resolve().parent
def baseline(plan,registry):
    q=A.read(A.B/'audit_reuse168/run001/result.json');assert q['status']=='all_actual158_corners_firstpass_best_and_complete_domain_independently_qualified' and q['complete_old111201_domain_preserved'] and q['pending115']==1120
    bind=A.read(A.B/'audit_reuse168/candidate_binding001.json');prod=A.read(plan['base158_candidate_path']);assert prod['full115_passed']==110081 and prod['target115_pending']==1120
    for p in (plan['base158_candidate_path'],plan['base158_updates_path']):
        h=bind['candidate_locks'][p];assert A.sha(p)==(h['sha256'] if isinstance(h,dict) else h)
    base={r['leaf_index']:r for r in A.read(plan['base158_updates_path'])};prior155=[];prior158={};seen=set();receivers=set();used=set();count=Counter()
    for sp in plan['phase2_corner_shards']:
        d=A.read(sp['input_path']);assert d['shard_id']==sp['shard_id']
        for k,r in d['receivers'].items():
            i=int(k);assert i not in receivers and i%15==sp['shard_id'];receivers.add(i);assert r['base_best158']==base[r['leaf_index']]['best'] and r['node']==base[r['leaf_index']]['node']
        for p in d['pairs']:
            cid=p['canonical_id'];assert p['pair_index'] not in seen and cid in registry;seen.add(p['pair_index']);used.add(cid)
            assert p['node']==d['receivers'][str(p['receiver_index'])]['node'] and d['sources'][cid]==registry[cid]['source171_metadata']
            r=p['prior_test_same_old155_q']
            if r:
                assert r['new_q_strict_decrease_required'];count[r['stage']]+=1
                if r['stage']=='155':prior155.append((p,r))
                else:assert r['stage']=='158';prior158.setdefault(r['journal_path'],[]).append((p,r))
    assert seen==set(range(35872)) and receivers==set(range(1120)) and used==set(registry) and count=={'155':13,'158':58}
    needed={r['source155_leaf_index'] for _,r in prior155};leafs={r['leaf_index']:r for r in A.lines(plan['source155_margin_part']) if r['leaf_index'] in needed}
    for p,r in prior155:
        row=leafs[r['source155_leaf_index']];assert row['leaf_index']==p['leaf_index'] and row['node']==p['node'] and row['canonical_id']==p['canonical_id'] and row['margin110_new']==r['margin110'] and row['margin115']==r['margin115']
    for path,items in prior158.items():
        assert A.sha(path)==items[0][1]['journal_sha256'];events={r['seq']:r for r in A.lines(path)}
        for p,r in items:
            e=events[r['event_seq']];assert e['kind']=='evaluated_single_source' and e['hash']==r['event_hash'] and e['canonical_id']==p['canonical_id']
            c=e['certificate'];assert c['node']==p['node'] and c['margin110']==r['margin110'] and c['margin115']==r['margin115']
    return {'status':'PASS_168_155_158_metadata_join','receivers':1120,'pairs':35872,'old155_guards':13,'old158_guards':58,'new_or_old_corner_recalculations':0}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan-sha',required=True);a=ap.parse_args();assert A.sha(H/'continuation_plan.json')==a.plan_sha;c=A.read(H/'continuation_plan.json');start=time.monotonic()
    for p,h in c['locks'].items():assert A.sha(p)==h
    prior=Path(c['successful_source_prefix_directory']);r=A.read(prior/'result.json');assert r['status']=='FAULT_prefix_preserved' and len(r['workers'])==15 and all(w['exitcode']==0 for w in r['workers'])
    original=A.read(H/'plan.json');registry=A.read(original['registry_path']);seen=set();rows=[]
    for i in range(15):
        result=A.read(prior/f'worker{i:02d}/result.json');assert result==r['outcomes'][str(i)] and result['status']=='complete'
        p=prior/f'worker{i:02d}/source_checks.jsonl';assert A.sha(p)==result['source_checks_sha256']
        for row in A.lines(p):assert row['status']=='PASS_identity_and_existing_qualification_binding' and row['canonical_id'] not in seen;seen.add(row['canonical_id']);rows.append(row)
    assert seen==set(registry) and len(rows)==455
    finite=A.read(original['finite_plan_path']);b=baseline(finite,registry)
    for p,h in c['locks'].items():assert A.sha(p)==h
    for p,h in original['source_locks'].items():assert A.sha(p)==h
    u=resource.getrusage(resource.RUSAGE_SELF)
    result={'status':'PASS_source455_and_metadata_join_complete','pid':os.getpid(),'seconds':time.monotonic()-start,'cpu_user':u.ru_utime,'cpu_system':u.ru_stime,'peak_rss_bytes':u.ru_maxrss,
      'continuation_plan_sha256':a.plan_sha,'source_audit_master_pid':r['master_pid'],'source_workers':r['workers'],'source_stage_counts':dict(Counter(x['stage'] for x in rows)),
      'original_H64_numeric_stage_counts':dict(Counter(x['original_jet_numeric_source_id'] for x in rows)),'sources':455,'positive_endpoint_intervals_checked':455*73,'legacy_raw_endpoint_arrays_exact_equal':411,
      'no_source_prefix_rerun':True,'source_audit_fault_retained':'metadata kind vocabulary expected evaluated,actual evaluated_single_source;no donor/jet failure',
      'baseline':b,'new_q_bandproducts_corners_H_log_gradient_calls':0,'all1076_source_locks_after_unchanged':True,'shared_H_numeric_truth_boundary':True}
    A.put(H/'binding_result.json',result);print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
