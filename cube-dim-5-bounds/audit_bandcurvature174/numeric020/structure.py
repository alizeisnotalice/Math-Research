"""Consume new independently checked results; no old corner/q re-evaluation."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from common import read,lines,sha,atomic
def finish(plan,O):
    O=Path(O);P=Path(plan['producer_output']);old=read(plan['finite_plan_path']);C=F(read(plan['common_log_path'])['upper']);q=read(O/'independent_q_registry.json');actualq=read(P/'source_q_registry.json');assert set(q)==set(actualq) and len(q)==455
    for cid,a in q.items():assert a['producer_row_hash']==actualq[cid]['row_hash'] and F(a['q'])==F(actualq[cid]['certificate']['q_new'])
    for phase in ('mass','corners'):
        r=read(P/f'actual_attempt0001_{phase}_exit_receipt.json');assert r['complete'] and not r['fault'] and len(r['workers'])==15 and all(p['exitcode']==0 and not p['alive'] for p in r['workers'])
    best={};receivers={};seenpairs=set();counts={'evaluated':0,'skip_after_firstPASS':0,'skip_oldpair_q_unchanged':0}
    for sp in plan['corner_shards']:
        data=read(sp['input_path']);b=read(O/'corners'/f'shard{sp["audit_shard"]:02d}'/'receiver_best.json');s=read(O/'corners'/f'shard{sp["audit_shard"]:02d}'/'summary.json');assert s['status']=='complete'
        for k,v in b.items():assert int(k) not in best;best[int(k)]=v;receivers[int(k)]=data['receivers'][k]
        for p in data['pairs']:assert p['pair_index'] not in seenpairs;seenpairs.add(p['pair_index'])
        for k in counts:counts[k]+=s['counts'][k]
    assert seenpairs==set(range(35872)) and set(best)==set(receivers)==set(range(1120)) and sum(counts.values())==35872
    assert counts['evaluated']==plan['actual_corner_cases'] and counts['skip_after_firstPASS']==1702 and counts['skip_oldpair_q_unchanged']==0
    base={r['leaf_index']:r for r in read(old['base158_updates_path'])};basecandidate=read(old['base158_candidate_path']);by_leaf={receivers[i]['leaf_index']:(receivers[i],b) for i,b in best.items()};assert len(by_leaf)==1120
    seen=set();roots=defaultdict(list);minimum=None;base_min=None;pending=[];updates=[];passed=improved=0
    for row in lines(old['source155_margin_part']):
        i=row['leaf_index'];root,path=row['node'];assert i not in seen and isinstance(path,str) and set(path)<={'0','1'};seen.add(i);roots[root].append(path);v=F(row['margin110_new']);source={'kind':'original155','canonical_id':row['canonical_id']}
        if i in base:
            b=base[i];assert b['node']==row['node'] and F(b['best']['best_margin110'])>=v;v=F(b['best']['best_margin110']);source={'kind':b['best']['best_source_kind'],'canonical_id':b['best']['best_canonical_id']}
        assert v>0;base_min=v if base_min is None else min(base_min,v)
        if i in by_leaf:
            r,b=by_leaf[i];assert r['node']==row['node'] and r['base_best158']==base[i]['best'] and F(r['base_best158']['best_margin110'])==v
            assert F(b['best_margin110'])>=v and F(b['best_margin115'])==F(b['best_margin110'])-C;improved+=F(b['best_margin110'])>v;v=F(b['best_margin110']);source={'kind':b['best_source_kind'],'canonical_id':b['best_canonical_id']};passed+=b['first_pass115'] is not None;updates.append({'leaf_index':i,'node':row['node'],'best':b})
            if b['first_pass115'] is None:pending.append({'leaf_index':i,'node':row['node'],'best_canonical_id':b['best_canonical_id'],'best_margin110':b['best_margin110'],'best_margin115':b['best_margin115']})
        else:assert v-C>0
        key=(v,i);minimum=(key,row['node'],source) if minimum is None or key<minimum[0] else minimum
    assert seen==set(range(111201)) and set(roots)==set(range(52360)) and passed+len(pending)==1120 and F(basecandidate['minimum_margin110_updated'])==base_min and minimum[0][0]>=base_min
    # Exact Kraft plus a different bottom-up binary-sibling closure check.
    for paths in roots.values():
        cells=set(paths);assert len(cells)==len(paths);depth=max(map(len,cells));assert sum(1<<(depth-len(p)) for p in cells)==1<<depth
        while cells!={''}:
            p=max(cells,key=len);assert p;other=p[:-1]+('1' if p[-1]=='0' else '0');assert other in cells and p[:-1] not in cells;cells.remove(p);cells.remove(other);cells.add(p[:-1])
    expected={'full111201_closed_domain_retained':True,'unchanged158_leaves':110081,'updated_receivers':1120,'best_rule':'max(base158,all actually tested new single-source wholebox margins including failed)','base158_minimum_margin110':str(base_min),'minimum_margin110_updated':str(minimum[0][0]),'minimum_ge_base158':True,'witness_leaf_index':minimum[0][1],'witness_node':minimum[1],'witness_source':minimum[2],'minimum_margin115':str(minimum[0][0]-C),'candidate_linear_floor':str(F(11,100000)*(1+minimum[0][0])),'new115_passes':passed,'target115_pending':len(pending),'full115_passed':110081+passed,'improved_receivers':improved,'source_q_strict_decreases':sum(a['strict'] for a in q.values()),'evaluated_pairs':counts['evaluated'],'firstPASS_skips':counts['skip_after_firstPASS'],'oldpair_q_unchanged_skips':counts['skip_oldpair_q_unchanged'],'main_P8_down192_separate':True,'formal_g_C5_update':False}
    assert expected==read(P/'complete_domain_candidate.json') and sorted(updates,key=lambda r:r['leaf_index'])==read(P/'receiver_best_updates.json') and sorted(pending,key=lambda r:r['leaf_index'])==read(P/'target115_pending.json')
    atomic(O/'independent_complete_domain.json',expected)
    r={'status':'PASS_fixed173v3_new_q_corners_full_domain','sources':455,'band_products':455*49,'corner_cases':counts['evaluated'],'corner_points':16*counts['evaluated'],'all_pair_events':35872,'firstPASS_skips':counts['skip_after_firstPASS'],'old71_guards_inherited_source_gate_and_actual_skip_rule_checked':True,'complete_leaves':111201,'closed_roots_Kraft_and_binary_merges':52360,'candidate':expected,'shared_old_H_endpoint_and_physical_jet_truth':True,'new_H_log_gradient_or_old_corner_q_replay':0,'formal_C5_acceptance':False};return r
