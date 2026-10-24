"""Exact audit of saved models, inheritance, event ledger and frontier partition.

Standard library only. No flint/checker import, model, integral or point call.
Designed for this normally returned finite run, including an incomplete cover.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from posixpath import normpath
from time import monotonic

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'saved_verification.json'
BASE_SHA = '21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
OLD_SHA = '6028e48efe444e27c5259858cf1eeb8b41e0081b6e39f06102439ddb198c4cac'
ORIGIN18_SHA = '1b545178105b1b344557e1b0540f0ce5af126555df62cd3cd021ca95556ab6b4'
SIGNED_SHA = 'd8e9e38f0c7d61c9b49ea7d5875ab823c85c02c58ba6ef242c8315471f7c0cdd'
SUP_SHA = '1cb0a6baed805ce17162179988ba6004acc49df49ef811c8c19a53bc1f7f3298'
FIXED_SHA = '4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05'
GRID = tuple(Q(s) for s in ('0','1/256','1/128','1/64','1/32','1/16','1/8','1/4',
    '3/8','1/2','5/8','3/4','7/8','15/16','31/32','63/64','127/128','255/256','1'))


def main():
    if OUT.exists():
        raise SystemExit('Existing saved audit; no overwrite or repeated verification.')
    start = monotonic()
    def guard():
        if monotonic() - start > 90:
            raise TimeoutError('90 second saved-data audit limit')
    def read(p):
        guard()
        return json.loads(p.read_text())
    def digest(p):
        guard()
        return sha256(p.read_bytes()).hexdigest()
    j, cover, run = (read(ROOT / name) for name in ('result.json','leaf_cover.json','extend.run.json'))
    assert run['status'] == 'terminal' and run['exit_code'] == 0 and run['elapsed_seconds'] < 300
    assert run['source_sha256_before'] == run['source_sha256_after'] == j['source_sha256'] == digest(ROOT/'extend_signed.py')
    assert run['log_sha256'] == digest(ROOT/'extend.log')
    assert j['cover_sha256'] == digest(ROOT/'leaf_cover.json')
    for inp in j['inputs']:
        assert digest(ROOT.parent/inp['file']) == inp['sha256']
    base_path, old_path, sup_path = (ROOT.parent/s for s in
        ('round01/scaled_certificate.json','round20/leaf_cover.json','round19/new_model.json'))
    assert digest(base_path) == BASE_SHA and digest(old_path) == OLD_SHA and digest(sup_path) == SUP_SHA
    base, old, supplemental = read(base_path), read(old_path), read(sup_path)
    c1 = Q(base['c'])
    assert c1 > 0 and cover['baseline_sha256'] == BASE_SHA and cover['parents'] == [[52,1]]
    assert tuple(map(Q,cover['grid'])) == GRID and cover['roots'] == old['roots'] and len(cover['roots']) == 324
    assert cover['source_cover_sha256'] == OLD_SHA and cover['source_Psi615_sha256'] == SUP_SHA
    assert cover['source_signed661_sha256'] == SIGNED_SHA
    assert cover['max_depth'] == 4 and cover['precision_bits'] == 192
    assert cover['historical_round18_Taylor_calls'] == 0 and 'Taylor_calls' not in cover
    nodes, old_nodes = cover['nodes'], old['nodes']
    assert cover['initial_inherited_nodes'] == len(old_nodes)
    assert all(n['id'] == i for i,n in enumerate(nodes))
    incoming = Counter()
    rects = {}
    immutable = ('id','parent','local','physical','depth','parent_id','root_id')
    old_pending = set(old['pending'])
    for i, n in enumerate(nodes):
        guard()
        assert n['parent'] == [52,1] and 0 <= n['depth'] <= 4
        x0,x1,y0,y1 = map(Q,n['local'])
        assert 0 <= x0 < x1 <= 1 and 0 <= y0 < y1 <= 1
        phys = ((52+x0)/32,(52+x1)/32,(1+y0)/32,(1+y1)/32)
        assert tuple(map(Q,n['physical'])) == phys and phys[0]+phys[2] > Q(33,20)
        rects[i] = (x0,x1,y0,y1)
        assert n['status'] in ('accepted','split','pending','unresolved_depth4')
        if i < len(old_nodes):
            before = old_nodes[i]
            assert all(n[k] == before[k] for k in immutable)
            if i not in old_pending and i != 661:
                assert n['status'] == before['status']
                assert n.get('children') == before.get('children')
                if 'model_file' in before:
                    assert n['model_sha256'] == before['model_sha256']
                    assert n['model_file'] == normpath('../round20/'+before['model_file'])
            if 'partial_model_file' in before:
                assert n['historical_partial_optional'] is True
                assert n['partial_model_sha256'] == before['partial_model_sha256']
                # Historical partial files are metadata only, never read or required.
                assert n['partial_model_file'] == normpath('../round20/'+before['partial_model_file'])
        if n['status'] == 'split':
            assert n['depth'] < 4 and len(n['children']) == 4
            xm,ym = (x0+x1)/2,(y0+y1)/2
            expected = ((x0,xm,y0,ym),(xm,x1,y0,ym),(x0,xm,ym,y1),(xm,x1,ym,y1))
            for cid, box in zip(n['children'],expected):
                z = nodes[cid]
                assert tuple(map(Q,z['local'])) == box and z['parent_id'] == i
                assert z['depth'] == n['depth']+1 and z['root_id'] == n['root_id']
                incoming[cid] += 1
        else:
            assert 'children' not in n
        if n['status'] == 'unresolved_depth4':
            assert n['depth'] == 4
    for k,rid in enumerate(cover['roots']):
        i,jj = divmod(k,18)
        assert rects[rid] == (GRID[i],GRID[i+1],GRID[jj],GRID[jj+1])
        assert nodes[rid]['depth'] == 0 and nodes[rid]['parent_id'] is None and nodes[rid]['root_id'] == rid
    roots = set(cover['roots'])
    assert all(incoming[i] == (0 if i in roots else 1) for i in range(len(nodes)))
    assert set(cover['pending']) == {n['id'] for n in nodes if n['status']=='pending'}
    assert len(cover['pending']) == len(set(cover['pending']))
    assert set(cover['unresolved']) == {n['id'] for n in nodes if n['status']=='unresolved_depth4'}
    assert cover['in_flight'] is None and j['in_flight'] == []

    model_cache, rows_cache, modes = {}, {}, Counter()
    rebuilt_rows = 0
    def bern(model):
        assert len(model['bernstein']) == 16
        e = Q(model['remainder_upper'])
        assert e >= 0
        if 'remainder_lower' in model:
            assert 0 <= Q(model['remainder_lower']) <= e
        vals = []
        for k,b in enumerate(model['bernstein']):
            assert b['index'] == list(divmod(k,4)) and Q(b['lower']) <= Q(b['upper'])
            vals.append(Q(b['lower'])-e)
        return vals
    def match_rows(rows, a, b):
        assert len(rows) == 16
        normalized = []
        for k,row in enumerate(rows):
            guard()
            assert row['bernstein_index'] == list(divmod(k,4))
            rho = c1*a[k]+b[k]
            assert (Q(row['c_column']),Q(row['s_column']),Q(row['baseline_slack'])) == (a[k],b[k],rho)
            normalized.append((a[k],b[k],rho))
        return normalized
    def signed_deduction(data, node, derivative_key, deduction_key):
        assert data['one_sided_only'] is True
        u0,u1,v0,v1 = map(Q,node['physical'])
        assert len(data[derivative_key]) == 2
        uppers = []
        for dr,(a,b) in zip(data[derivative_key],((u0,v0),(v0,u0))):
            lo,hi,au = map(Q,(dr['lower'],dr['upper'],dr['absolute_upper']))
            assert (Q(dr['a']),Q(dr['b']),dr['kind']) == (a,b,'Hab')
            assert lo <= hi and au >= max(abs(lo),abs(hi))
            uppers.append(au)
        muv = sum(uppers,Q(0));du,dv = (u1-u0)/2,(v1-v0)/2
        assert Q(data['Muv_upper']) == muv and Q(data['du']) == du and Q(data['dv']) == dv
        deduction = muv*du*dv
        assert deduction >= 0 and Q(data[deduction_key]) == deduction
        return deduction
    def signed_coefficients(data,node):
        assert data['unscaled_by_c'] and data['Bernstein_degree'] == [3,3]
        assert 'remainder_upper' not in data
        deduction = signed_deduction(data,node,'mixed_derivatives','one_sided_lower_deduction')
        assert set(data['center_values']) == {'F0','Fu','Fv'}
        for v in data['center_values'].values(): assert Q(v['lower']) <= Q(v['upper'])
        assert len(data['bernstein']) == 16
        a = []
        for k,v in enumerate(data['bernstein']):
            assert v['index'] == list(divmod(k,4)) and Q(v['lower']) <= Q(v['upper'])
            a.append(Q(v['lower'])-deduction)
        return a
    for n in nodes:
        guard()
        if 'model_file' not in n:
            assert n['status'] == 'pending'
            continue
        p = (ROOT/n['model_file']).resolve()
        assert digest(p) == n['model_sha256']
        m = read(p)
        assert m['node_id'] == n['id'] and m['parent'] == n['parent'] and m['depth'] == n['depth']
        assert m['local'] == n['local'] and m['physical'] == n['physical'] and m['baseline_sha256'] == BASE_SHA
        b = bern(m['AS'])
        if m['mode'] == 'const':
            assert n['id'] < len(old_nodes) and n['id'] not in old_pending
            assert Q(m['const']['lower']) <= Q(m['const']['upper'])
            a = [Q(m['const']['lower'])]*16
        elif m['mode'] == 'outside_poly':
            assert m['fixed_hlcert_sha256'] == FIXED_SHA
            assert m['Psi']['unscaled_by_c'] and m['Psi']['Bernstein_degree'] == [3,3]
            a = bern(m['Psi'])
        else:
            assert m['mode'] == 'outside_signed' and m['fixed_hlcert_sha256'] == FIXED_SHA
            a = signed_coefficients(m['Psi'],n)
        reconstructed = match_rows(m['rows'],a,b)
        assert m['accepted'] == all(rho > 0 for _,_,rho in reconstructed)
        if m['mode'] == 'const':
            assert m['accepted'] == (n['id'] != 615 and old_nodes[n['id']]['status']=='accepted')
        elif m['mode'] == 'outside_poly':
            assert m['accepted'] == (old_nodes[n['id']]['status']=='accepted')
        else:
            assert m['accepted'] == (n['status']=='accepted')
        model_cache[n['id']] = m
        rows_cache[n['id']] = reconstructed
        modes[m['mode']] += 1
        rebuilt_rows += 16
    n615 = nodes[615]
    assert n615['status'] == 'accepted' and n615['mode'] == 'outside_poly_reused'
    assert n615['previous_status'] == 'unresolved_depth4'
    assert digest(ROOT/n615['supplemental_Psi_file']) == n615['supplemental_Psi_sha256'] == SUP_SHA
    assert supplemental['source_AS_model_sha256'] == n615['model_sha256']
    assert supplemental['source_cover_sha256'] == ORIGIN18_SHA and supplemental['source_hlcert_sha256'] == FIXED_SHA
    assert supplemental['baseline_sha256'] == BASE_SHA
    assert supplemental['local'] == n615['local'] and supplemental['physical'] == n615['physical']
    assert supplemental['unscaled_by_c'] and supplemental['Bernstein_degree'] == [3,3]
    ep = Q(supplemental['Psi_remainder_upper'])
    assert 0 <= Q(supplemental['Psi_remainder_lower']) <= ep
    a = []
    for k,b in enumerate(supplemental['Psi_bernstein']):
        assert b['index'] == list(divmod(k,4)) and Q(b['lower']) <= Q(b['upper'])
        a.append(Q(b['lower'])-ep)
    rows_cache[615] = match_rows(supplemental['rows'],a,bern(model_cache[615]['AS']))
    assert all(rho > 0 for _,_,rho in rows_cache[615])
    rebuilt_rows += 16
    n661 = nodes[661]
    assert n661['status'] == 'accepted' and n661['mode'] == 'outside_signed_reused'
    assert n661['previous_status'] == 'unresolved_depth4'
    signed_path = ROOT/n661['supplemental_signed_file']
    assert digest(signed_path) == n661['supplemental_signed_sha256'] == SIGNED_SHA
    signed = read(signed_path)
    assert signed['source_model_sha256'] == n661['model_sha256']
    assert signed['source_cover_sha256'] == OLD_SHA and signed['fixed_hlcert_sha256'] == FIXED_SHA
    assert signed['baseline_sha256'] == BASE_SHA
    assert signed['local'] == n661['local'] and signed['physical'] == n661['physical']
    assert signed['reused_T1_Bernstein'] == model_cache[661]['Psi']['bernstein']
    deduction = signed_deduction(signed,n661,'derivatives','Psi_one_sided_lower_deduction')
    a = [Q(v['lower'])-deduction for v in model_cache[661]['Psi']['bernstein']]
    rows_cache[661] = match_rows(signed['rows'],a,bern(model_cache[661]['AS']))
    assert all(rho > 0 for _,_,rho in rows_cache[661])
    rebuilt_rows += 16
    leaf_ids = [leaf['node_id'] for leaf in cover['leaves']]
    assert len(leaf_ids) == len(set(leaf_ids)) and set(leaf_ids) == {n['id'] for n in nodes if n['status']=='accepted'}
    for leaf in cover['leaves']:
        guard()
        n = nodes[leaf['node_id']]
        for field in ('parent','local','physical','depth','mode','model_file','model_sha256'):
            assert leaf[field] == n[field]
        expected = rows_cache[n['id']]
        assert all(rho > 0 for _,_,rho in expected)
        actual = [(Q(v['c_column']),Q(v['s_column']),Q(v['baseline_slack'])) for v in leaf['rows']]
        assert actual == expected and [v['bernstein_index'] for v in leaf['rows']] == [list(divmod(k,4)) for k in range(16)]
        if n['id'] == 615:
            assert leaf['supplemental_Psi_file'] == n['supplemental_Psi_file']
            assert leaf['supplemental_Psi_sha256'] == SUP_SHA
        if n['id'] == 661:
            assert leaf['supplemental_signed_file'] == n['supplemental_signed_file']
            assert leaf['supplemental_signed_sha256'] == SIGNED_SHA
    for leaf in old['leaves']:
        inherited = next(z for z in cover['leaves'] if z['node_id']==leaf['node_id'])
        assert inherited['rows'] == leaf['rows'] and inherited['mode'] == leaf['mode']

    assert digest(ROOT/'calls.jsonl') == j['calls_sha256']
    active, started, terminal, by_kind = [], {}, {}, Counter()
    hint_by_node, integral_by_node, processing_order = Counter(), Counter(), []
    events = 0
    with (ROOT/'calls.jsonl').open() as stream:
        for line in stream:
            guard();event = json.loads(line);events += 1
            cid,kind,nid = event['id'],event['kind'],event['node_id']
            assert kind in ('AS','Psi','Hint','acb_integral')
            if event['status'] == 'started':
                assert cid == len(started) and cid not in started
                if kind in ('AS','Psi'):
                    assert not active
                    if kind == 'AS': processing_order.append(nid)
                elif kind == 'Hint':
                    assert active and started[active[-1]]['kind']=='Psi' and started[active[-1]]['node_id']==nid
                    hint_by_node[nid] += 1
                else:
                    assert active and started[active[-1]]['kind']=='Hint'
                    parent = started[active[-1]]
                    assert parent['node_id']==nid and parent['Hint_id']==event['Hint_id']
                    integral_by_node[nid] += 1
                started[cid] = event;active.append(cid);by_kind[kind] += 1
            else:
                assert event['status']=='completed' and active.pop()==cid and cid not in terminal
                assert all(event[k]==v for k,v in started[cid].items() if k!='status')
                assert event['elapsed_seconds'] >= 0
                terminal[cid] = event
    assert not active and set(started)==set(terminal) and events==2*len(started)
    assert len(started)==j['calls_started_total']
    for kind in ('AS','Psi','Hint','acb_integral'):
        assert j['counts'][kind] == dict(started=by_kind[kind],completed=by_kind[kind],failed=0,interrupted=0)
    per_node = defaultdict(Counter)
    for event in started.values(): per_node[event['node_id']][event['kind']] += 1
    assert all(v['AS']==v['Psi']==1 and v['Hint']<=8 and v['acb_integral']<=v['Hint'] for v in per_node.values())
    pending = list(old['pending'])
    for nid in processing_order:
        assert pending.pop() == nid
        if nodes[nid]['status']=='split': pending.extend(reversed(nodes[nid]['children']))
    assert pending == cover['pending'] and len(processing_order)==len(set(processing_order))
    integral_per_hint = Counter(e['Hint_id'] for e in started.values() if e['kind']=='acb_integral')
    for nid in processing_order:
        u0,u1,v0,v1 = map(Q,nodes[nid]['physical']);uc,vc = (u0+u1)/2,(v0+v1)/2
        expected = [(uc,vc,'H'),(vc,uc,'H'),(uc,vc,'Ha'),(vc,uc,'Hb'),
                    (uc,vc,'Hb'),(vc,uc,'Ha'),(u0,v0,'Hab'),(v0,u0,'Hab')]
        hints = [e for e in started.values() if e['kind']=='Hint' and e['node_id']==nid]
        assert len(hints) == 8
        for e,(a,b,kind) in zip(hints,expected):
            assert (Q(e['a']),Q(e['b']),e['derivative_kind']) == (a,b,kind)
            assert integral_per_hint[e['Hint_id']] == (1 if b < Q(33,20) else 0)
    new_models = {i for i,m in model_cache.items() if m['mode']=='outside_signed'}
    assert new_models == set(processing_order)
    assert j['new_boxes']['started']==j['new_boxes']['completed']==len(new_models)
    assert j['new_boxes']['failed']==j['new_boxes']['interrupted']==0
    assert j['new_boxes']['accepted']==sum(nodes[i]['status']=='accepted' for i in new_models)
    assert j['new_boxes']['splits']==sum(nodes[i]['status']=='split' for i in new_models)
    assert j['AS_models_saved']==j['Psi_models_saved']==len(new_models)
    assert len(new_models)<=1500 and by_kind['Hint']<=12000 and by_kind['acb_integral']<=12000
    assert all(j[k]==0 for k in ('full_absolute_Psi_calls','new_const_calls','new_T_calls','new_point_calls','LP_calls','full_checker_calls'))
    assert j['accepted_leaves']==len(leaf_ids) and j['pending_count']==len(cover['pending'])
    assert j['unresolved_count']==len(cover['unresolved']) and j['complete']==cover['complete']

    frontier = [n['id'] for n in nodes if n['status']!='split']
    xs = sorted({Q(0),Q(1)} | {x for i in frontier for x in rects[i][:2]})
    for left,right in zip(xs,xs[1:]):
        guard();cursor = Q(0)
        intervals = sorted((y0,y1) for i in frontier for x0,x1,y0,y1 in [rects[i]] if x0<=left and right<=x1)
        for bottom,top in intervals:
            assert bottom == cursor
            cursor = top
        assert cursor == 1
    area = lambda ids: sum(((rects[i][1]-rects[i][0])*(rects[i][3]-rects[i][2]) for i in ids),Q(0))
    assert area(frontier) == 1
    accepted_area = area(leaf_ids)
    if cover['complete']:
        assert not cover['pending'] and not cover['unresolved'] and accepted_area==1
    if j['status']=='depth_cap_unresolved':
        assert len(cover['unresolved'])==1 and processing_order[-1]==cover['unresolved'][0]
    worst_new = min(((rho,i,k) for i in new_models for k,(_,_,rho) in enumerate(rows_cache[i])),default=None)
    cross_R = sum(Q(nodes[i]['physical'][0])<Q(33,20)<Q(nodes[i]['physical'][1]) or
        Q(nodes[i]['physical'][2])<Q(33,20)<Q(nodes[i]['physical'][3]) for i in new_models)
    output = dict(status='passed_saved_arithmetic_tree_frontier_audit',
        scope='Exact saved enclosures/rows, immutable origins, ledger and frontier partition; no independent model or integral evaluation.',
        complete=cover['complete'],old_const_models_checked=modes['const'],old_absolute_poly_models_checked=modes['outside_poly'],new_signed_models_checked=modes['outside_signed'],
        supplemental_models_checked=2,rows_reconstructed=rebuilt_rows,accepted_leaves=len(leaf_ids),
        accepted_rows_checked=16*len(leaf_ids),accepted_area=str(accepted_area),
        pending_area=str(area(cover['pending'])),unresolved_area=str(area(cover['unresolved'])),
        frontier_nodes=len(frontier),frontier_area='1',slabs_checked=len(xs)-1,
        status_counts=dict(Counter(n['status'] for n in nodes)),call_counts=dict(by_kind),trace_events=events,
        new_boxes_crossing_coordinate_R=cross_R,
        minimum_accepted_slack=str(min(rho for i in leaf_ids for _,_,rho in rows_cache[i])),
        worst_new_row=None if worst_new is None else dict(slack=str(worst_new[0]),node=worst_new[1],index=list(divmod(worst_new[2],4))),
        elapsed_seconds=monotonic()-start,result_sha256=digest(ROOT/'result.json'),
        cover_sha256=digest(ROOT/'leaf_cover.json'),run_sha256=digest(ROOT/'extend.run.json'),
        calls_sha256=j['calls_sha256'],verifier_sha256=digest(Path(__file__)))
    OUT.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output))


if __name__ == '__main__':
    main()
