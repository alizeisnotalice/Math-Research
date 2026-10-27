"""Repo-only parameter reconstruction; numerical verification is explicit opt-in.

Importing this module performs no mathematical evaluation. Parameter views use
JSON fields, string concatenation, lookup, and SHA only. Raw-state traceability
is inherited from the frozen actual reader31006, not replayed here.
"""
from collections import defaultdict
import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_BASE = HERE.parents[1]


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        while part := handle.read(8 << 20):
            h.update(part)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_bytes())


def rows(path):
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            yield json.loads(line)


def unique(rows_, field):
    result = {}
    for row in rows_:
        key = row[field]
        assert key not in result
        result[key] = row
    return result


def repository_parameters(base=REPO_BASE):
    """Pure field reconstruction. No Fraction, exp, q, division, or corners."""
    base = Path(base)
    portable = base/'upper_global151'
    manifest = load(portable/'manifest.json')

    def dataset(name):
        for part in manifest['datasets'][name]:
            yield from rows(portable/part['path'])

    jets = unique(dataset('jets'), 'jet_id')
    coefs = unique(dataset('coefficients'), 'canonical_id')
    aliases151 = unique(dataset('provenance'), 'source_alias_index')
    leaves = unique(dataset('leaves'), 'leaf_index')
    by_alias = {}
    for leaf in leaves.values():
        key = leaf['source_alias_index']
        assert key not in by_alias
        by_alias[key] = leaf
    supports = {}
    for support in dataset('supports'):
        key = tuple(support['node'])
        assert key not in supports
        supports[key] = support
    alias_choice = {}
    aliases149 = {}
    # Original155 selection: first alias, replaced only by the first curvature
    # alias when the existing choice has no curvature proof. Dataset151 retained
    # original149 row order; full equality is checked against the frozen inputs.
    for key, alias in aliases151.items():
        restored = dict(alias)
        restored['saved_g110_margin'] = by_alias[key]['claimed_minimum_log_margin']
        restored['target110'] = coefs[alias['canonical_id']]['core']['target110']
        aliases149[key] = restored
        cid = alias['canonical_id']
        if cid not in alias_choice or (alias_choice[cid]['curvature_proof'] is None and alias['curvature_proof'] is not None):
            alias_choice[cid] = alias
    qrows = unique(rows(base/'upper_retarget155/packet/new_q-0001.jsonl.gz'), 'canonical_id')
    marginrows = unique(rows(base/'upper_retarget155/packet/leaf_margins-0001.jsonl.gz'), 'leaf_index')
    common = load(base/'upper_retarget155/packet/common_log_certificate.json')
    source_views = {}
    for cid, coefficient in coefs.items():
        core = coefficient['core']
        jet = jets[coefficient['jet_id']]
        alias = alias_choice[cid]
        original_alias = aliases149[alias['source_alias_index']]
        mass = support = None
        if core['curvature_rule_class'] == '129_anchor_fixed_Jensen':
            support = supports[tuple(alias['curvature_proof']['node'])]
            mass = {
                'lower':'21/200000','lower_epsilon':'1/64','upper_epsilon':'1/128',
                'upper':jet['values'][0][1],'point':core['point'],'support_node':support['node'],
                'support_source_leaf_sha256':support['source_leaf_sha256'],
                'support_source_state_sha256':support['source_state_sha256'],
                'support_certificate_sha256':digest(support['certificate']),
                'original_support_certificate_sha256':support['original_certificate_sha256'],
                'source_alias_index':alias['source_alias_index'],
                'source_alias_sha256':digest(original_alias),'same_point_same_model':True,
                'point_in_support_closedbox':True,'no_new_raw_lower_lookup':True,
                'global_H128_floor_not_used_as_H64':True}
        source_views[cid] = {
            'canonical_id':cid,'core':core,'jet':jet,
            'old_H_interval':jet['values'][0],
            'typed_cache_key':[core['epsilon'],','.join(core['point'])],
            'raw_record_sha256':coefficient['raw_record_sha256'],
            'q_new':qrows[cid]['q_new'],'saved_q155':qrows[cid],
            'mass_qualification':mass,'support':support,
            'proof_alias151':alias,'proof_alias149':original_alias,
            'common_log_upper':common['upper']}
    leaf_views = {}
    for idx, leaf in leaves.items():
        leaf_views[idx] = {'leaf151':leaf,'saved_leaf155':marginrows[idx],
                           'alias151':aliases151[leaf['source_alias_index']]}
    support_views = {key:{'support':support,'jet':jets[support['jet_id']]}
                     for key,support in supports.items()}
    assert len(source_views) == 107530 and len(leaf_views) == 111201 and len(support_views) == 3812
    return source_views, leaf_views, support_views, common


def original_source_view(task, common_upper):
    """Select exactly the input fields consumed by mathematical verification."""
    s = task['source155']
    coefficient = task['coefficient']
    assert coefficient['core'] == s['old_core']
    return {'canonical_id':s['canonical_id'],'core':s['old_core'],'jet':task['jet'],
            'old_H_interval':s['old_H_interval'],'typed_cache_key':s['typed_cache_key'],
            'raw_record_sha256':s['raw_record_sha256'],
            'q_new':task['saved_q155']['q_new'],'saved_q155':task['saved_q155'],
            'mass_qualification':s['mass_qualification'],'support':task['support'],
            'proof_alias151':task['proof_alias151'],'proof_alias149':task['proof_alias149'],
            'common_log_upper':common_upper}


def original_leaf_view(task):
    old, leaf = task['source_leaf155'], task['leaf151']
    for field in ('leaf_index','node','canonical_id','source_alias_index','epsilon','box_lower','box_upper','closed_box'):
        assert old[field] == leaf[field]
    assert old['saved_g110_margin'] == leaf['claimed_minimum_log_margin']
    return {key:task[key] for key in ('leaf151','saved_leaf155','alias151')}


def aggregate(views):
    h = hashlib.sha256()
    for key in sorted(views):
        h.update(encoded(views[key])+b'\n')
    return h.hexdigest()


def mathematical_worker(spec):
    """Future opt-in only: reuse actual180 arithmetic functions verbatim."""
    from fractions import Fraction as F
    p = HERE/'evidence/reader180.py'
    modspec = importlib.util.spec_from_file_location('frozen_actual180_pure_math',p)
    core = importlib.util.module_from_spec(modspec)
    modspec.loader.exec_module(core)
    sources, leaves, supports, common = spec
    upper = core.common_log(common)
    table = core.roots()
    for task in supports:
        s, jet = task['support'], task['jet']
        c = core.support_core(s,jet)
        point,delta,gamma,_ = core.numeric_core(c,jet,F(21,200000),'Delta_lower')
        margin,_ = core.full_corners(point,delta,gamma,F(256),core.box(s['node'],table))
        assert margin == F(s['certificate']['minimum_log_margin']) > 0
    parameters = {}
    used = set()
    for v in sources:
        c, jet = v['core'],v['jet']
        assert core.digest(c) == v['canonical_id'] and jet['model'] == core.MODEL
        assert c['epsilon'] == jet['epsilon'] and v['old_H_interval'] == jet['values'][0]
        assert v['raw_record_sha256'] == jet['raw_record_sha256']
        point,delta,gamma,U = core.numeric_core(c,jet,F(11,100000),'Delta110_lower')
        e,old = F(c['epsilon']),F(c['penalty_coefficient'])
        A = 1/(16*e*e)
        mass = v['mass_qualification']
        if mass is None:
            assert c['curvature_rule_class'] == '075_full_frame' and old == A
            L = F(0)
        else:
            assert e == F(1,128) and c['curvature_rule_class'] == '129_anchor_fixed_Jensen'
            L = F(mass['lower'])
            assert L == F(21,200000) and F(mass['upper']) == U and mass['point'] == c['point']
            assert old == A*(1-F(3,4)*L/U)
            lo,hi = core.box(mass['support_node'],table)
            assert all(a<=x<=b for a,x,b in zip(lo,point,hi))
            used.add(tuple(mass['support_node']))
        assert 0 <= L < U
        q = min(old,F(2,3)*A,A*(F(7,10)-F(8,15)*L/U))
        assert q == F(v['q_new']) and q > 0
        parameters[v['canonical_id']] = point,delta,gamma,q,c['epsilon']
    minimum = None
    passed = 0
    for task in leaves:
        leaf,saved = task['leaf151'],task['saved_leaf155']
        point,delta,gamma,q,eps = parameters[leaf['canonical_id']]
        assert leaf['closed_box'] and eps == leaf['epsilon']
        lo,hi = core.box(leaf['node'],table)
        assert lo == list(map(F,leaf['box_lower'])) and hi == list(map(F,leaf['box_upper']))
        margin,bits = core.full_corners(point,delta,gamma,q,(lo,hi))
        assert margin == F(saved['margin110_new']) and bits == saved['worst_corner_bits']
        assert margin >= F(leaf['claimed_minimum_log_margin']) > 0
        assert F(saved['margin115']) == margin-upper and saved['passed115'] == (margin-upper>0)
        candidate = (margin,leaf['leaf_index'],leaf['node'],eps)
        if minimum is None or candidate[:2] < minimum[:2]:
            minimum = candidate
        passed += margin-upper>0
    return {'sources':len(sources),'supports':len(supports),'leaves':len(leaves),
            'minimum':[str(minimum[0]),*minimum[1:]],'passed115':passed,'used_supports':[list(n) for n in sorted(used)]}


def verify_mathematics(base):
    """Explicit future run; never invoked by archival/equivalence operations."""
    from fractions import Fraction as F
    import multiprocessing as mp
    import time
    start = time.monotonic()
    assert sha(HERE/'evidence/reader180.py') == 'd80c831a253fa3f97bd961da3b83bceb60b7bf4effe4caf5dd1d65be0f8d4900'
    lock = load(HERE/'repository_locks.json')
    for relative, expected in lock.items():
        assert sha(Path(base)/relative) == expected['sha256']
    receipt = load(HERE/'equivalence_result.json')
    assert receipt['status'] == 'all_parameter_views_exactly_equal_no_new_mathematics'
    sources,leaves,supports,common = repository_parameters(base)
    for label,views in [('sources',sources),('leaves',leaves),('supports',supports)]:
        assert aggregate(views) == receipt['parameter_view_sha256'][label]
    jobs = [([],[],[],common) for _ in range(15)]
    owners = {}
    for cid,view in sources.items():
        sid = int(cid[:16],16)%15
        owners[cid] = sid
        jobs[sid][0].append(view)
    for view in leaves.values():
        jobs[owners[view['leaf151']['canonical_id']]][1].append(view)
    for i,view in enumerate(supports.values()):
        jobs[i%15][2].append(view)
    with mp.get_context('spawn').Pool(15) as pool:
        results = pool.map(mathematical_worker,jobs)
    expected = load(HERE/'evidence/run001/result.json')
    minimum = min((F(r['minimum'][0]),*r['minimum'][1:]) for r in results)
    assert [str(minimum[0]),*minimum[1:]] == expected['minimum']
    assert sum(r['sources'] for r in results) == 107530
    assert sum(r['supports'] for r in results) == 3812
    assert sum(r['leaves'] for r in results) == 111201
    assert sum(r['passed115'] for r in results) == 109961
    used = {tuple(n) for r in results for n in r['used_supports']}
    assert len(used) == 3812
    modspec = importlib.util.spec_from_file_location('frozen_actual180_forest',HERE/'evidence/reader180.py')
    core = importlib.util.module_from_spec(modspec)
    modspec.loader.exec_module(core)
    forest = core.forest({'portable':str(Path(base)/'upper_global151')})
    print(json.dumps({'status':'repo_only_mathematical_verification_complete',
        'minimum':[str(minimum[0]),*minimum[1:]],'forest':forest,'workers':15,
        'seconds':time.monotonic()-start,'raw_state_trace':'inherited actual31006; not replayed',
        'independent_original_H_enclosures':False,'P8_down192':'inherited178; not replayed'}))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-base',type=Path,default=REPO_BASE)
    ap.add_argument('--verify-mathematics',action='store_true')
    args = ap.parse_args()
    if args.verify_mathematics:
        verify_mathematics(args.repo_base)
    else:
        print('No mathematics executed. Future explicit verification: python3 portable180.py --verify-mathematics')
