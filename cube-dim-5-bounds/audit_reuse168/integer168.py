"""Independent existing-source/new-box corners. Negative110 is valid."""
from math import gcd
from rational168 import Q


def lcm(a,b):
    return a//gcd(a,b)*b


def one_source(source,child,existing_log_upper):
    c = source['core'];q = Q(c['q_new']);delta = Q(c['Delta_lower'])
    assert c['target'] == '23/200000' and c['epsilon'] in ('1/64','1/128')
    assert c['q_rule'] == '153_time_shape_uniform_or_original_g105_mass' and q > Q(0)
    assert delta+existing_log_upper == Q(source['saved_Delta110_lower'])
    assert Q(c['point'][0]) == 1
    point = [Q(x) for x in c['point'][1:]]
    low,high = [Q(x) for x in child['box_lower']],[Q(x) for x in child['box_upper']]
    gamma = [(Q(a),Q(b)) for a,b in c['loggradient_directed']]
    assert len(point) == len(low) == len(high) == len(gamma) == 4
    assert all(a <= b for a,b in gamma) and all(a < b for a,b in zip(low,high))
    W = 1
    for x in point+low+high:W = lcm(W,x.d)
    anchor = [x.n*(W//x.d) for x in point]
    sides = [(a.n*(W//a.d),b.n*(W//b.d)) for a,b in zip(low,high)]
    T = lcm(delta.d,q.d*W*W)
    for a,b in gamma:
        T = lcm(T,a.d*W);T = lcm(T,b.d*W)
    di = delta.n*(T//delta.d);qi = q.n*(T//(q.d*W*W))
    gi = [(a.n*(T//(a.d*W)),b.n*(T//(b.d*W))) for a,b in gamma]
    best = None
    for mask in range(16):
        bits = [(mask>>(3-j))&1 for j in range(4)]
        v = [sides[j][bits[j]]-anchor[j] for j in range(4)]
        value = di+sum(min(a*x,b*x) for (a,b),x in zip(gi,v))-qi*sum(x*x for x in v)
        if best is None or value < best[0]:best = value,bits
    margin115 = Q(best[0],T)
    # This add-back is valid only because the saved Delta identity was checked.
    margin110 = margin115+existing_log_upper
    return margin115,margin110,best[1]


def initialize(child):
    old = child['prior155_failed115']
    assert Q(old['margin110_new']) > Q(0) and Q(old['margin115']) <= Q(0)
    return {'best_margin110':old['margin110_new'],'best_margin115':old['margin115'],
            'best_source_kind':'original155','best_canonical_id':old['canonical_id'],
            'best_certificate_ref':{'source155_leaf_index':child['leaf_index']},
            'tested_pairs':0,'first_pass115':None,'status':'old_full_domain_retained'}


def accept_test(best,certificate):
    assert best['first_pass115'] is None
    best['tested_pairs'] += 1
    if Q(certificate['margin110']) > Q(best['best_margin110']):
        best.update({'best_margin110':certificate['margin110'],'best_margin115':certificate['margin115'],
                     'best_source_kind':'new158_single_anchor','best_canonical_id':certificate['canonical_id'],
                     'best_certificate_ref':certificate})
    if certificate['passed115']:
        best['first_pass115'] = certificate;best['status'] = 'first115PASS_stop_receiver'
    else:
        best['status'] = 'best_of_old155_and_all_tested_newsource_wholebox_margins'
    assert Q(best['best_margin110']) > Q(0)
