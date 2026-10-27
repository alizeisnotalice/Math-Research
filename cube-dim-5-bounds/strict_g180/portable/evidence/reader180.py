"""Second G-certificate reader. Stdlib only; no historical math imports.

The preparation mode only copies/repartitions saved JSON fields and hashes.
Only --actual executes the finite mathematical verification.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import queue
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
DESKTOP = Path('/Users/zhengzhihao/Desktop/T')
MAIN = DESKTOP / 'Math-Research-cube-dim-5/cube-dim-5-bounds'
MODEL = 'n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
WORKERS = 15
CHUNK = 256


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


def write_new(path, value):
    with Path(path).open('xb') as handle:
        handle.write(encoded(value) + b'\n')


def atomic(path, value):
    tmp = Path(str(path) + '.tmp')
    with tmp.open('wb') as handle:
        handle.write(encoded(value) + b'\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def rows(path):
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            yield json.loads(line)


def prepare():
    """No interval arithmetic, q, logs, exp, corners, or numeric comparisons."""
    out = HERE / 'inputs'
    out.mkdir(exist_ok=True)
    assert not (HERE/'plan.json').exists()
    portable = BASE / 'upper_global151'
    m = load(portable / 'manifest.json')
    meta = load(portable / 'metadata.json')
    paths = set([portable / 'manifest.json', portable / 'metadata.json'])

    def dataset(kind):
        for part in m['datasets'][kind]:
            p = portable / part['path']
            paths.add(p)
            yield from rows(p)

    jets = {r['jet_id']: r for r in dataset('jets')}
    coefs = {r['canonical_id']: r for r in dataset('coefficients')}
    supports = {tuple(r['node']): r for r in dataset('supports')}
    provenance = {r['source_alias_index']: r for r in dataset('provenance')}
    leaves = list(dataset('leaves'))
    qpath = BASE / 'upper_retarget155/packet/new_q-0001.jsonl.gz'
    lpath = BASE / 'upper_retarget155/packet/leaf_margins-0001.jsonl.gz'
    saved_q = {r['canonical_id']: r for r in rows(qpath)}
    saved_l = {r['leaf_index']: r for r in rows(lpath)}
    paths.update([qpath, lpath])
    aliaspath = DESKTOP / 'cube-finite-dim-scratch/upper_retarget149/source_provenance.jsonl'
    aliases149 = {}
    with aliaspath.open() as handle:
        for line in handle:
            r = json.loads(line)
            assert r['source_alias_index'] not in aliases149
            aliases149[r['source_alias_index']] = r
    paths.add(aliaspath)
    finitepath = DESKTOP / 'cube-finite-dim-scratch/upper_retarget155/finite_plan.json'
    finite = load(finitepath)
    paths.add(finitepath)
    sources = {}
    source_leaves = {}
    fullq = {}
    for part in finite['shards']:
        p = Path(part['input_path'])
        assert sha(p) == part['input_sha256']
        paths.add(p)
        data = load(p)
        qjournal = p.parent / 'run001' / f'shard{part["shard_id"]:02}' / 'q_updates.jsonl'
        paths.add(qjournal)
        previous = '0'*64
        with qjournal.open() as handle:
            for seq, line in enumerate(handle, 1):
                event = json.loads(line)
                unhashed = {k:v for k,v in event.items() if k != 'hash'}
                assert event['previous'] == previous and event['seq'] == seq and digest(unhashed) == event['hash']
                previous = event['hash']
                cid = event['body']['canonical_id']
                assert cid not in fullq
                fullq[cid] = event
        for r in data['sources']:
            cid = r['canonical_id']
            assert cid not in sources
            sources[cid] = r
        for r in data['leaves']:
            idx = r['leaf_index']
            assert idx not in source_leaves
            source_leaves[idx] = r
    shards = [{'sources': [], 'leaves': [], 'supports': []} for _ in range(WORKERS)]
    assignment = {cid: int(cid[:16], 16) % WORKERS for cid in coefs}
    assert len(jets) == 109446 and len(coefs) == len(sources) == len(saved_q) == len(fullq) == 107530
    assert len(leaves) == len(source_leaves) == len(saved_l) == len(provenance) == len(aliases149) == 111201
    assert len(supports) == 3812
    for cid, c in coefs.items():
        s = sources[cid]
        assert c['core'] == s['old_core']
        mass = s['mass_qualification']
        support = None if mass is None else supports[tuple(mass['support_node'])]
        shards[assignment[cid]]['sources'].append({
            'coefficient': c, 'jet': jets[c['jet_id']], 'source155': s,
            'saved_q155': saved_q[cid], 'support': support,
            'full_q_event155': fullq[cid],
            'proof_alias151': provenance[s['old_proof_alias_index']],
            'proof_alias149': aliases149[s['old_proof_alias_index']]})
    for i, (_, r) in enumerate(sorted(supports.items())):
        shards[i % WORKERS]['supports'].append({'support': r, 'jet': jets[r['jet_id']]})
    for r in leaves:
        idx = r['leaf_index']
        shards[assignment[r['canonical_id']]]['leaves'].append({
            'leaf151': r, 'source_leaf155': source_leaves[idx],
            'saved_leaf155': saved_l[idx], 'alias151': provenance[r['source_alias_index']]})
    specs = []
    for sid, data in enumerate(shards):
        p = out / f'chunk-input-{sid:02}.json.gz'
        if p.exists():
            with gzip.open(p,'rt') as handle:
                assert json.load(handle) == data
        else:
            with p.open('xb') as raw:
                with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as z:
                    z.write(encoded(data) + b'\n')
        paths.add(p)
        specs.append({'worker': sid, 'path': str(p),
                      'sources': len(data['sources']), 'leaves': len(data['leaves']),
                      'supports': len(data['supports'])})
    for stage in meta['numeric_stages'].values():
        paths.add(Path(stage['production_state']))
    paths.update([
        BASE / 'upper_retarget155/packet/common_log_certificate.json',
        MAIN / 'round018/primary_curvature155_acceptance.json',
        MAIN / 'round018/primary_curvature155_proof.md',
        MAIN / 'round019/primary_upper166_acceptance.json',
        BASE / 'audit_retarget157/integer157.py',
        BASE / 'audit_retarget157/run001/public_log_check.json',
        BASE / 'upper_retarget155/packet/theory_sources.json'])
    paths.add(portable/'data/documents-0001.jsonl.gz')
    for directory, relative in [
        ('Math-Research-cube5-astra-lower','audit_frame075/report.md'),
        ('Math-Research-cube5-auditcurvature129-sol','audit_curvature129/report.md'),
        ('Math-Research-cube5-curvature128-sol','upper_curvature153/report.md'),
        ('Math-Research-cube5-auditcurvature129-sol','audit_curvature154/report.md')]:
        paths.add(DESKTOP/directory/'cube-dim-5-bounds'/relative)
    inherit = DESKTOP / 'Math-Research-cube5-secondaudit178-sol/cube-dim-5-bounds/second_audit178'
    paths.update(inherit / p for p in ['inspect166178.py', 'inspection166_registration.json', 'inspection166_result.json'])
    locks = {str(p): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(paths)}
    plan = {
        'task': 'strict_g180', 'workers': WORKERS, 'chunk_size': CHUNK,
        'code_sha256': sha(__file__), 'specs': specs, 'sources': locks,
        'portable': str(portable), 'numeric_stages': meta['numeric_stages'],
        'mathematics_started': False, 'new_H_DD_gradient_calls': 0,
        'no_mathematical_import_from_old_code': True, 'wallcap': None,
        'automatic_retry': False,
        'contract': {'sources': 107530, 'supports': 3812, 'leaves': 111201,
                     'jets_raw_bindings': 109446, 'ordered_roots': 52360, 'splits': 58841,
                     'full_leaf_corner_values': 1779216, 'support_corner_values': 60992,
                     'delta_exp_checks': 111342, 'exp_fixed_bits': 1024, 'exp_degree': 512,
                     'delta_absolute_bound': '16', 'public_log_objects': 1},
        'algorithm': 'raw state→typed jet→H-positive interval division and exp(Delta) certificate→153 q; four coordinate endpoint factors→16 whole-anchor sums; binary sibling merging per ordered root',
        'scope': 'H128 on F only, then G>=H128 and symmetry/-4 homogeneity; no fixed-epsilon homogeneity',
        'shared_dependency': 'Original rigorous H/physical derivative enclosures are inherited, including H64_inherited052; this reader checks all needed raw bytes and certificate arithmetic, not a second enclosure engine.',
        'exp_proof': 'd=0 exact; d>0 ceil dyadic Taylor terms, tail<=term512*d/513/(1-d/514); d<0 floor positive Taylor for exp(-d), reciprocal is exp(d) upper. |d|>16 or unresolved comparison fails and preserves prefix; never drops sources.',
        'cost': '15 workers; 111342 fixed 512-step integer exp checks, 107530 Fraction q/gradient checks, 115013 boxes with 16 factored Fraction sums. Sequential raw JSON/hash pass about 1.62GB plus portable/shards; no deadline, no new grid or kernel. Peak memory dominated by one original 765MB JSON at a time plus typed jets; workers each load one ~1/15 shard.',
        'inherited_scalar': '178 inspect166178.py actual PID25410 exit0, exact saved P8 recurrence/down192 and U166 power comparisons, bound by three source locks; read its code/result, do not execute it.'}
    write_new(HERE / 'plan.json', plan)
    print(json.dumps({'prepared_plan_sha256': sha(HERE/'plan.json'),
                      'code_sha256': sha(__file__), 'math_started': False,
                      'specs': specs, 'locked_files': len(locks)}), flush=True)


def exact(value):
    assert not isinstance(value, (float, bool))
    return F(value)


def ceildiv(a, b):
    assert a >= 0 and b > 0
    return (a + b - 1) // b


def qualify_delta(delta, ratio):
    """exp(delta)<=ratio; independently directed integer Taylor, no log."""
    assert ratio > 0 and abs(delta) <= 16, ('finite_delta_contract', str(delta))
    if delta == 0:
        assert ratio >= 1
        return
    scale = 1 << 1024
    positive = delta > 0
    t = abs(delta)
    a, b = t.numerator, t.denominator
    term = total = scale
    for k in range(1, 513):
        term = ceildiv(term*a, b*k) if positive else (term*a)//(b*k)
        total += term
    if positive:
        tail = ceildiv(term*a*514, 513*(b*514-a))
        assert (total+tail)*ratio.denominator <= scale*ratio.numerator, ('Delta_not_qualified', str(delta), str(ratio))
    else:
        assert total*ratio.numerator >= scale*ratio.denominator, ('negative_Delta_not_qualified', str(delta), str(ratio))


def numeric_core(core, jet, target, delta_field):
    assert core['epsilon'] == jet['epsilon'] and core['point'] == jet['point']
    assert len(core['point']) == 5 and exact(core['point'][0]) == 1
    values = [[exact(x) for x in pair] for pair in jet['values']]
    assert len(values) == 5 and all(a <= b for a, b in values)
    lower, upper = values[0]
    assert 0 < lower <= upper
    gamma = [[exact(x) for x in pair] for pair in core['loggradient_directed']]
    assert len(gamma) == 4
    for j, ((a, b), (gl, gu)) in enumerate(zip(values[1:], gamma)):
        assert gl <= gu
        quotients = [a/lower, a/upper, b/lower, b/upper]
        assert gl <= min(quotients) <= max(quotients) <= gu, {
            'gate': 'physical_gamma_division', 'coordinate': j,
            'saved': list(map(str, [gl, gu])),
            'independent_quotient': list(map(str, [min(quotients), max(quotients)])),
            'action': 'preserve_and_diagnose_correlation_or_arithmetic_no_clip'}
    delta = exact(core[delta_field])
    qualify_delta(delta, lower/target)
    return [exact(x) for x in core['point'][1:]], delta, gamma, upper


def roots():
    # Ordered lower grid indices; permutation symmetry extends their cover to F.
    return [(a, b, c, d) for a in range(32) for b in range(a+1)
            for c in range(b+1) for d in range(c+1)]


def box(node, root_table):
    rid, path = node
    assert type(rid) is int and 0 <= rid < 52360 and isinstance(path, str)
    assert all(bit in '01' for bit in path)
    lo = [F(a, 32) for a in root_table[rid]]
    hi = [a+F(1, 32) for a in lo]
    for depth, bit in enumerate(path):
        axis = depth % 4
        middle = (lo[axis]+hi[axis])/2
        if bit == '0':
            hi[axis] = middle
        else:
            lo[axis] = middle
    return lo, hi


def full_corners(point, delta, gamma, q, bounds):
    assert q > 0
    lo, hi = bounds
    factors = []
    for j in range(4):
        pair = []
        for endpoint in (lo[j], hi[j]):
            displacement = endpoint-point[j]
            pair.append(min(gamma[j][0]*displacement,
                            gamma[j][1]*displacement)-q*displacement**2)
        factors.append(pair)
    values = [(delta, ())]
    for pair in factors:
        values = [(value+pair[b], bits+(b,)) for value, bits in values for b in (0, 1)]
    assert len(values) == 16
    winner = min(values)
    assert winner[0] == delta+sum(min(pair) for pair in factors)
    return winner[0], list(winner[1])


def support_core(support, jet):
    c = support['certificate']
    assert c['epsilon'] == '1/64' and c['target'] == '21/200000'
    assert c['checked_points'] == 16 and c['passed'] is True
    assert c['point_key'] == ','.join(jet['point'])
    assert c['penalty_coefficient'] == '256'
    return {'epsilon': c['epsilon'], 'point': jet['point'],
            'Delta_lower': c['Delta_lower'], 'loggradient_directed': c['loggradient_directed']}


def common_log(c):
    # Reverse Horner of a polynomial in z^2, rather than integer-denominator sum.
    assert c['ratio'] == '23/22' and c['terms'] == 64
    z = F(1, 45)
    poly = F(0)
    for k in range(63, -1, -1):
        poly = F(1, 2*k+1)+z*z*poly
    lower = 2*z*poly
    tail = 2*z**129/(129*(1-z*z))
    assert exact(c['lower']) == lower and exact(c['tail_upper']) == tail
    upper = exact(c['upper'])
    assert upper == lower+tail and upper > 0
    return upper


def validate_raw(plan):
    portable = Path(plan['portable'])
    manifest = load(portable/'manifest.json')
    jets = {}
    counts = Counter()
    for part in manifest['datasets']['jets']:
        for jet in rows(portable/part['path']):
            view = dict(jet)
            h = view.pop('minimal_projection_sha256')
            assert digest(view) == h and jet['jet_id'] not in jets
            assert jet['model'] == MODEL and jet['epsilon'] in ('1/64', '1/128')
            pk = ','.join(jet['point'])
            assert digest({'epsilon': jet['epsilon'], 'point_key': pk,
                           'raw_record_sha256': jet['raw_record_sha256']}) == jet['jet_id']
            jets[jet['jet_id']] = jet
    assert len(jets) == 109446
    bystage = defaultdict(dict)
    for jet in jets.values():
        bystage[jet['numeric_source_id']][','.join(jet['point'])] = jet
    supports = {}
    for part in manifest['datasets']['supports']:
        for s in rows(portable/part['path']):
            n = tuple(s['node'])
            assert n not in supports
            supports[n] = s
    assert len(supports) == 3812
    support_found = set()
    for name, stage in plan['numeric_stages'].items():
        assert plan['sources'][stage['production_state']]['sha256'] == stage['production_state_sha256']
        state = load(stage['production_state'])
        field = state[stage['cache_field']]
        cache = field if stage['cache_field'] in ('cache', 'new_cache') else field[stage['epsilon']]
        for key, jet in bystage[name].items():
            record = cache[key]
            assert digest(record) == jet['raw_record_sha256']
            assert record['point'] == jet['point'] and record['bits'] == jet['bits']
            assert record['values'] == jet['values'] and stage['epsilon'] == jet['epsilon']
            assert ('model' in record) == jet['raw_model_field_present']
            assert ('epsilon' in record) == jet['raw_epsilon_field_present']
            if 'model' in record:
                assert record['model'] == MODEL
            else:
                assert jet['model_tag_source'] == '6cb840b9c1f9690792578b15d553216c0ea2061e044c2702e582f810f97b63b4'
            if 'epsilon' in record:
                assert exact(record['epsilon']) == exact(jet['epsilon'])
            else:
                assert stage['epsilon'] == '1/64' and jet['epsilon_tag_source'] == 'immutable source outertypedcache 1/64'
            counts[name] += 1
        if name == 'H128_old098':
            assert sha(stage['production_state']) == stage['production_state_sha256']
            for row in state['passed']:
                n = tuple(row['node'])
                if n not in supports:
                    continue
                s = supports[n]
                assert n not in support_found and digest(row) == s['source_leaf_sha256']
                assert row['certificate'] == s['original_certificate']
                assert digest(row['certificate']) == s['original_certificate_sha256']
                assert s['source_state_sha256'] == stage['production_state_sha256']
                assert row['leaf_epsilon'] == '1/64' and row['target'] == '21/200000'
                c, original = s['certificate'], row['certificate']
                for key in ('Delta_lower', 'loggradient_directed', 'minimum_log_margin', 'checked_points', 'passed'):
                    assert c[key] == original[key]
                pk = original.get('point_key', original.get('pair_cache_key', [None,None])[1]) or original.get('anchor')
                assert s['typed_cache_key'] == ['1/64', pk]
                raw = state['cache_by_epsilon']['1/64'][pk]
                assert digest(raw) == s['raw_record_sha256']
                assert raw['values'] == jets[s['jet_id']]['values']
                support_found.add(n)
        del cache, field, state
    assert sum(counts.values()) == 109446 and len(support_found) == 3812
    return {'jets_raw_bound': sum(counts.values()), 'stage_counts': dict(counts),
            'original_support_leaves_bound': len(support_found),
            'independent_H_enclosure_engine': False}


def forest(plan):
    portable = Path(plan['portable'])
    manifest = load(portable/'manifest.json')
    sets = defaultdict(set)
    allnodes = set()
    indices = set()
    for part in manifest['datasets']['leaves']:
        for leaf in rows(portable/part['path']):
            node = tuple(leaf['node'])
            assert node not in allnodes and leaf['closed_box'] is True
            assert leaf['leaf_index'] not in indices
            allnodes.add(node)
            sets[node[0]].add(node[1])
            indices.add(leaf['leaf_index'])
    assert indices == set(range(111201)) and len(sets) == 52360
    splits = 0
    for rid in range(52360):
        paths = sets[rid]
        ancestors = {p[:j] for p in paths for j in range(len(p))}
        assert not paths & ancestors
        # Independent bottom-up binary sibling collapse, no Kraft sum.
        for parent in sorted(ancestors, key=lambda p: (-len(p), p)):
            assert parent+'0' in paths and parent+'1' in paths
            paths.remove(parent+'0')
            paths.remove(parent+'1')
            paths.add(parent)
            splits += 1
        assert paths == {''}
    assert splits == 58841
    return {'leaves': 111201, 'ordered_roots': 52360, 'splits': splits,
            'method': 'all binary sibling pairs collapsed to all ordered roots',
            'extension': 'tail permutation symmetry; clamp floor(32 y_j) at31 includes upper face'}


def worker(spec, output, logupper_text, result_queue):
    out = Path(output)
    out.mkdir()
    write_new(out/'registration.json', {'pid': os.getpid(), 'worker': spec['worker'],
               'input_sha256': sha(spec['path']), 'counts': spec, 'wallcap': None})
    start = time.monotonic()
    completed = {'supports': 0, 'sources': 0, 'leaves': 0}
    chain = '0'*64
    journal = (out/'chunks.jsonl').open('x')
    records = []
    try:
        with gzip.open(spec['path'], 'rt') as handle:
            data = json.load(handle)
        table = roots()
        logupper = exact(logupper_text)
        qtable = {}
        minimum = None
        below256 = negative_delta = outside_box = passed115 = 0
        used_supports = set()

        def checkpoint(phase, payload, force=False):
            nonlocal chain, records
            records.append(payload)
            if len(records) < CHUNK and not force:
                return
            body = {'previous': chain, 'phase': phase, 'completed': dict(completed), 'records': records}
            chain = digest(body)
            journal.write(json.dumps({'body': body, 'sha256': chain}, sort_keys=True, separators=(',', ':'))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
            atomic(out/'checkpoint.json', {'completed': dict(completed), 'chain': chain,
                                          'pid': os.getpid(), 'status': 'running', 'automatic_retry': False})
            records = []

        for task in data['supports']:
            s, jet = task['support'], task['jet']
            assert s['jet_id'] == jet['jet_id'] and s['raw_record_sha256'] == jet['raw_record_sha256']
            c = support_core(s, jet)
            point, delta, gamma, _ = numeric_core(c, jet, F(21,200000), 'Delta_lower')
            margin, bits = full_corners(point, delta, gamma, F(256), box(s['node'], table))
            assert margin == exact(s['certificate']['minimum_log_margin']) > 0
            completed['supports'] += 1
            checkpoint('supports', {'node': s['node'], 'margin': str(margin), 'bits': bits})
        checkpoint('supports', {'phase_finished': True}, True)
        for task in data['sources']:
            cobj, jet, s = task['coefficient'], task['jet'], task['source155']
            core, cid = cobj['core'], cobj['canonical_id']
            assert digest(core) == cid and cid not in qtable and cid == s['canonical_id']
            assert core == s['old_core'] and cobj['jet_id'] == jet['jet_id']
            assert cobj['raw_record_sha256'] == s['raw_record_sha256'] == jet['raw_record_sha256']
            assert s['typed_cache_key'] == [core['epsilon'], ','.join(core['point'])]
            assert s['old_H_interval'] == jet['values'][0] and core['target110'] == '11/100000'
            alias = task['proof_alias151']
            original_alias = task['proof_alias149']
            assert digest(original_alias) == s['old_proof_alias_sha256']
            assert {k:v for k,v in original_alias.items() if k not in ('saved_g110_margin','target110')} == alias
            assert alias['canonical_id'] == cid and alias['raw_record_sha256'] == jet['raw_record_sha256']
            assert alias['source_alias_index'] == s['old_proof_alias_index']
            point, delta, gamma, upper = numeric_core(core, jet, F(11,100000), 'Delta110_lower')
            negative_delta += delta < 0
            e = exact(core['epsilon'])
            A = 1/(16*e*e)
            old = exact(core['penalty_coefficient'])
            mass = s['mass_qualification']
            if mass is None:
                assert core['curvature_rule_class'] == '075_full_frame' and old == A
                lower = F(0)
                assert task['support'] is None
            else:
                assert e == F(1,128) and core['curvature_rule_class'] == '129_anchor_fixed_Jensen'
                support = task['support']
                node = tuple(mass['support_node'])
                assert tuple(support['node']) == node and alias['curvature_proof']['node'] == list(node)
                assert mass['lower_epsilon'] == '1/64' and mass['upper_epsilon'] == '1/128'
                assert mass['point'] == core['point'] and exact(mass['upper']) == upper
                lower = exact(mass['lower'])
                assert lower == F(21,200000) and old == A*(1-F(3,4)*lower/upper)
                assert mass['support_source_leaf_sha256'] == support['source_leaf_sha256']
                assert mass['support_source_state_sha256'] == support['source_state_sha256']
                assert mass['support_certificate_sha256'] == digest(support['certificate'])
                assert mass['original_support_certificate_sha256'] == support['original_certificate_sha256']
                assert mass['source_alias_index'] == s['old_proof_alias_index']
                assert mass['source_alias_sha256'] == digest(original_alias)
                assert mass['point_in_support_closedbox'] and mass['no_new_raw_lower_lookup'] and mass['global_H128_floor_not_used_as_H64']
                lo, hi = box(node, table)
                assert all(a <= x <= b for a,x,b in zip(lo, point, hi))
                used_supports.add(node)
            assert 0 <= lower < upper
            q = min(old, F(2,3)*A, A*(F(7,10)-F(8,15)*lower/upper))
            assert 0 < q <= old
            saved = task['saved_q155']
            event = task['full_q_event155']
            assert digest({k:v for k,v in event.items() if k != 'hash'}) == event['hash'] == saved['source_event_hash']
            fullq = event['body']
            assert digest(fullq) == saved['source_body_sha256']
            for key,value in [('canonical_id',cid),('epsilon',core['epsilon']),('point',core['point']),
                              ('Delta110_lower',core['Delta110_lower']),('loggradient_directed',core['loggradient_directed']),
                              ('mass_qualification',mass),('raw_record_sha256',jet['raw_record_sha256'])]:
                assert fullq[key] == value
            assert exact(fullq['Delta115_lower']) == delta-logupper
            assert exact(fullq['q_new']) == q and exact(fullq['q_old']) == old
            assert saved['canonical_id'] == cid and exact(saved['q_new']) == q
            assert saved['q_new_rule'] == '153_time_shape_uniform_or_original_g105_mass'
            qtable[cid] = (point, delta, gamma, q, core['epsilon'])
            below256 += q < 256
            completed['sources'] += 1
            checkpoint('sources', {'canonical_id': cid, 'q': str(q), 'Delta110': str(delta),
                                   'Delta115': str(delta-logupper)})
        checkpoint('sources', {'phase_finished': True}, True)
        for task in data['leaves']:
            leaf, sourceleaf, saved, alias = (task[k] for k in ['leaf151','source_leaf155','saved_leaf155','alias151'])
            idx, cid = leaf['leaf_index'], leaf['canonical_id']
            assert saved['leaf_index'] == idx and saved['canonical_id'] == cid and saved['node'] == leaf['node']
            for key in ['leaf_index','canonical_id','node','source_alias_index','epsilon','box_lower','box_upper','closed_box']:
                assert sourceleaf[key] == leaf[key]
            assert sourceleaf['saved_g110_margin'] == leaf['claimed_minimum_log_margin']
            assert alias['source_alias_index'] == leaf['source_alias_index'] and alias['canonical_id'] == cid and alias['node'] == leaf['node']
            point, delta, gamma, q, eps = qtable[cid]
            assert eps == leaf['epsilon'] and leaf['closed_box'] is True
            lo, hi = box(leaf['node'], table)
            assert list(map(exact,leaf['box_lower'])) == lo and list(map(exact,leaf['box_upper'])) == hi
            outside_box += not all(a <= x <= b for a,x,b in zip(lo,point,hi))
            margin, bits = full_corners(point,delta,gamma,q,(lo,hi))
            assert margin == exact(saved['margin110_new']) and bits == saved['worst_corner_bits']
            assert margin >= exact(leaf['claimed_minimum_log_margin']) > 0
            shifted = margin-logupper
            assert exact(saved['margin115']) == shifted and saved['passed115'] == (shifted > 0)
            candidate = (margin,idx,leaf['node'],eps)
            if minimum is None or candidate[:2] < minimum[:2]:
                minimum = candidate
            passed115 += shifted > 0
            completed['leaves'] += 1
            checkpoint('leaves', {'leaf_index': idx, 'margin110': str(margin), 'margin115': str(shifted), 'bits': bits})
        checkpoint('leaves', {'phase_finished': True}, True)
        assert all(completed[k] == spec[k] for k in completed)
        result = {'status': 'complete', 'pid': os.getpid(), 'worker': spec['worker'],
                  **completed, 'below256': below256, 'negative_delta': negative_delta,
                  'anchors_outside_leaf': outside_box, 'used_supports': [list(n) for n in sorted(used_supports)],
                  'passed115': passed115, 'minimum': [str(minimum[0]),*minimum[1:]],
                  'seconds': time.monotonic()-start, 'journal_chain': chain}
        atomic(out/'checkpoint.json', {'completed': completed,'chain':chain,'status':'complete','pid':os.getpid()})
        write_new(out/'result.json', result)
        result_queue.put(result)
    except BaseException as error:
        journal.flush()
        failure = {'status':'failed_preserved','pid':os.getpid(),'worker':spec['worker'],
                   'completed':completed,'journal_chain':chain,'pending_records':records,
                   'exception':repr(error),'traceback':traceback.format_exc(),'automatic_retry':False}
        write_new(out/'failure.json',failure)
        result_queue.put(failure)
        raise
    finally:
        journal.close()


def actual(plan_sha, out):
    planpath = HERE/'plan.json'
    assert sha(planpath) == plan_sha
    plan = load(planpath)
    assert sha(__file__) == plan['code_sha256']
    assert out.is_absolute() and out.parent == HERE and not out.exists()
    start = time.monotonic()
    for path, lock in plan['sources'].items():
        assert Path(path).stat().st_size == lock['bytes'] and sha(path) == lock['sha256'], path
    out.mkdir()
    write_new(out/'registration.json', {'pid':os.getpid(),'plan_sha256':plan_sha,
              'code_sha256':sha(__file__),'workers':WORKERS,'wallcap':None,'automatic_retry':False,
              'contract':plan['contract'],'started_unix':time.time()})
    try:
        raw = validate_raw(plan)
        write_new(out/'raw_binding.json', raw)
        cover = forest(plan)
        write_new(out/'forest.json', cover)
        upper = common_log(load(BASE/'upper_retarget155/packet/common_log_certificate.json'))
        write_new(out/'common_log.json', {'upper':str(upper),'method':'reverse Horner exact polynomial+positive geometric tail'})
        ctx = mp.get_context('spawn')
        rq = ctx.Queue()
        processes = []
        for spec in plan['specs']:
            p = ctx.Process(target=worker,args=(spec,str(out/f'worker{spec["worker"]:02}'),str(upper),rq))
            p.start()
            processes.append(p)
        results = {}
        failure_seen = False
        while len(results) < WORKERS:
            try:
                r = rq.get(timeout=.5)
                assert r['worker'] not in results
                results[r['worker']] = r
                failure_seen |= r['status'] != 'complete'
            except queue.Empty:
                # An abnormal dead worker must not trigger an automatic retry.
                if not any(p.is_alive() for p in processes):
                    failure_seen = True
                    break
        if failure_seen:
            for p in processes:
                if p.is_alive():
                    p.join()
            while True:
                try:
                    r = rq.get_nowait()
                    results.setdefault(r['worker'],r)
                except queue.Empty:
                    break
            write_new(out/'failure.json', {'status':'failed_preserved','results':results,
                      'exitcodes':[p.exitcode for p in processes], 'automatic_retry':False,
                      'seconds':time.monotonic()-start})
            raise RuntimeError('finite reader failed; outputs retained; no automatic retry')
        for p in processes:
            p.join()
        assert all(p.exitcode == 0 for p in processes)
        rq.close()
        rq.join_thread()
        assert sum(r['sources'] for r in results.values()) == 107530
        assert sum(r['supports'] for r in results.values()) == 3812
        assert sum(r['leaves'] for r in results.values()) == 111201
        uses = {tuple(n) for r in results.values() for n in r['used_supports']}
        assert len(uses) == 3812
        minimum = min((exact(r['minimum'][0]),*r['minimum'][1:]) for r in results.values())
        accepted = load(MAIN/'round018/primary_curvature155_acceptance.json')
        assert minimum == (exact(accepted['minimum_margin110']),accepted['minimum_witness_leaf_index'],accepted['minimum_witness_node'],accepted['minimum_witness_epsilon'])
        assert sum(r['passed115'] for r in results.values()) == 109961
        for path, lock in plan['sources'].items():
            assert sha(path) == lock['sha256'], path
        write_new(out/'result.json', {
            'status':'complete_finite_G_certificate_reader','pid':os.getpid(),
            'worker_pids':[p.pid for p in processes],'exitcodes':[p.exitcode for p in processes],
            'sources':107530,'supports':3812,'raw_jets':109446,'leaves':111201,
            'full_leaf_corners':1779216,'support_corners':60992,
            'delta_exp_checks':111342,'minimum':[str(minimum[0]),*minimum[1:]],
            'q_below256':sum(r['below256'] for r in results.values()),
            'negative_delta_sources':sum(r['negative_delta'] for r in results.values()),
            'anchors_outside_leaf':sum(r['anchors_outside_leaf'] for r in results.values()),
            'passed115':109961,'pending115':1240,'forest':cover,
            'seconds':time.monotonic()-start,'all_input_hashes_before_after_equal':True,
            'new_H_DD_gradient_calls':0,'new_grid_or_split_calls':0,'old_math_imports':0,
            'independent_original_H_enclosures':False,
            'P8_down192': 'inherits locked178 actual exact field check; no replay',
            'worker_result_sha256':{str(i):sha(out/f'worker{i:02}/result.json') for i in range(WORKERS)}})
        print(json.dumps({'status':'complete','pid':os.getpid(),'seconds':time.monotonic()-start,
                          'result_sha256':sha(out/'result.json')}),flush=True)
    except BaseException as error:
        if not (out/'failure.json').exists():
            write_new(out/'failure.json', {'status':'failed_preserved','pid':os.getpid(),
                      'exception':repr(error),'traceback':traceback.format_exc(),'automatic_retry':False})
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--actual', action='store_true')
    parser.add_argument('--reviewed-plan-sha')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    assert args.prepare != args.actual
    if args.prepare:
        prepare()
    else:
        assert args.reviewed_plan_sha and args.out
        actual(args.reviewed_plan_sha,args.out)
