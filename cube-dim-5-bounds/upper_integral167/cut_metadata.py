#!/usr/bin/env python3
"""New certificate-cut metadata selection; no factor/mean/chord recalculation.

Reads saved candidates, compares saved means/self fields, and checks local
child-sum selection using saved effective integers. Does not recalculate UH.
"""
from pathlib import Path
from collections import Counter
import gzip,hashlib,json,os
H=Path(__file__).resolve().parent
B=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_integral156')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def write(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def main():
    assert not (H/'cut_registration.json').exists()
    reg=read(S/'preflight002/registration.json');run=S/'run001'
    files=[B/'upper_path112/numeric_hierarchy.json.gz',B/'upper_path112/nodes.json.gz',B/'upper_path112/node_certificates.json.gz',
      S/'preflight002/registration.json',Path(reg['reuse150']['path']),run/'effective_changes.json',run/'results.json']
    results=[]
    for i in range(15):
        p=run/f'worker{i:02d}/result.json';files.append(p);wr=read(p);assert wr['status']=='complete';results.append(wr)
        for x in wr['chunks']:
            assert sha(x['path'])==x['sha256'];files.append(Path(x['path']))
    locks={str(p):sha(p) for p in files}
    write(H/'cut_registration.json',dict(task='167_candidate_cut_metadata_only',pid=os.getpid(),code_sha256=sha(__file__),source_locks=locks,
        exact_mean_factor_chord_power_gain_UH_recalculation=False,new_H_R_calls=0,cut_tie_policy='choose_children',candidate_acceptance='not inferred'))
    tree=read(B/'upper_path112/numeric_hierarchy.json.gz');ns=tree['nodes'];rows={r['key']:r for r in read(B/'upper_path112/nodes.json.gz')};certs={r['key']:r for r in read(B/'upper_path112/node_certificates.json.gz')}
    strict={k for k,n in ns.items() if n['selfU']==certs[k]['new_selfU']<certs[k]['old_selfU']}
    special='5:31,31,31,31';assert set(ns)-strict=={special}
    updates={};mi_count=Counter();missing=[]
    for wr in results:
        count=0
        for ch in wr['chunks']:
            rr=read(ch['path']);assert len(rr)==ch['count'];count+=len(rr)
            for r in rr:
                k=r['key'];assert k not in updates and k in strict
                assert r['old_mean_i']==certs[k]['new_M_i'] and r['new_mean_i']<r['old_mean_i']
                assert r['old_selfU']==ns[k]['selfU'] and r['new_selfU']<=r['old_selfU']
                assert r['packet_sha256']==rows[k]['packet_sha256']
                updates[k]=r['new_selfU'];mi_count['156_new_band_strict_means']+=1
        assert count==wr['completed']
    reused=read(reg['reuse150']['path'])
    for r in reused:
        k=r['key'];assert k not in updates and k in strict and not ns[k]['children']
        assert r['old_mean_i']==certs[k]['new_M_i'] and r['new_mean_i']<r['old_mean_i']
        assert r['old_selfU']==ns[k]['selfU'] and r['new_selfU']<=r['old_selfU']
        updates[k]=r['new_selfU'];mi_count['150_reused_strict_means']+=1
    assert set(updates)==strict and len(updates)==97971
    changed=read(run/'effective_changes.json');eff={r['key']:r['new_effU'] for r in changed};old_eff={r['key']:r['old_effU'] for r in changed}
    assert len(eff)==len(ns)==97972
    branch={};local=Counter()
    for k,n in ns.items():
        assert old_eff[k]==n['effU']
        selfU=updates.get(k,n['selfU']);assert eff[k]<=n['effU']
        if n['children']:
            for c in n['children']:assert ns[c]['parent']==k and ns[c]['level']==n['level']+1
            childsum=sum(eff[c] for c in n['children'])
            if eff[k]==childsum:
                assert childsum<=selfU;branch[k]='children'
            else:
                assert eff[k]==selfU<childsum;branch[k]='self'
        else:
            assert eff[k]==selfU;branch[k]='self'
        local[branch[k]]+=1
    selected=set();expanded=set();visited=set()
    def visit(k):
        assert k not in visited;visited.add(k)
        if branch[k]=='children':
            expanded.add(k)
            for c in ns[k]['children']:visit(c)
        else:selected.add(k)
    for k in tree['roots']:visit(k)
    assert special in expanded and special not in selected and selected<=strict
    assert ns[special]['parent'] is None
    all_desc=set()
    def descendants(k):
        if not ns[k]['children']:all_desc.add(k)
        else:
            for c in ns[k]['children']:descendants(c)
    for k in tree['roots']:descendants(k)
    assert all_desc==set(tree['frontier'])
    info=dict(status='candidate_stored_cut_metadata_pass_not_mathematical_acceptance',pid=os.getpid(),tree_nodes=len(ns),strict112_new_chord_self_nodes=len(strict),new_self_updates=len(updates),
      current_means_strictly_below_m112=dict(mi_count),all_updated_means_source_is_new_band_sum=True,legacy_mean_min_branch_inactive=True,
      local_saved_effective_branch_counts=dict(local),selected_cut_self_nodes=len(selected),selected_cut_expanded_nodes=len(expanded),selected_cut_levels=dict(Counter(ns[k]['level'] for k in selected)),
      selected_cut_uses_only_strict112_self=True,selected_cut_has_no_legacy080_082_088_self=True,selected_cut_drops_all_descendants_of_selected_self=True,
      special=dict(key=special,branch=branch[special],old_selfU=ns[special]['selfU'],old_effU=ns[special]['effU'],new_effU=eff[special],
        children=len(ns[special]['children']),new_childsum=sum(eff[c] for c in ns[special]['children']),new_eff_below_old_self=eff[special]<ns[special]['selfU']),
      selected_cut_key_sha256=hashlib.sha256(json.dumps(sorted(selected),separators=(',',':')).encode()).hexdigest(),
      all_source_locks_unchanged=all(sha(p)==h for p,h in locks.items()),global_UH_recomputed=False,new_H_R_factor_mean_chord_power_gain_calls=0,old_successful_baseline_replays=0)
    write(H/'cut_metadata.json',info);print(json.dumps(info,sort_keys=True))
if __name__=='__main__':main()
