"""Exact audit of saved models, inheritance, event ledger and frontier partition.

Standard library only. No flint/checker import, model, integral or point call.
Designed for this normally returned finite run, including an incomplete cover.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from math import comb
import json
from pathlib import Path
from posixpath import normpath
from time import monotonic

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'saved_verification.json'
BASE_SHA = '21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
OLD_SHA = 'e78b2bf4e6f5e00974a180cb85dfa9044ed42939930c4add48078e6c4064c8e7'
ORIGIN20_SHA = '6028e48efe444e27c5259858cf1eeb8b41e0081b6e39f06102439ddb198c4cac'
PIECES_SHA = '1921a3b1e009005ecc5e3bea6309cedc8932e5f8b2862401232974d27eef1581'
ALGEBRA_SHA = 'fd16ccf0309a4ae4627049ff730b2bb233be3948013d25cd514824febcc107b2'
SOURCE_SHA = '4bc087786a59772d8cc559aab411ada8446836a2609e67724ee0bf3d39f35c33'
ORIGIN18_SHA = '1b545178105b1b344557e1b0540f0ce5af126555df62cd3cd021ca95556ab6b4'
SIGNED_SHA = 'd8e9e38f0c7d61c9b49ea7d5875ab823c85c02c58ba6ef242c8315471f7c0cdd'
SUP_SHA = '1cb0a6baed805ce17162179988ba6004acc49df49ef811c8c19a53bc1f7f3298'
FIXED_SHA = '4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05'
GRID = tuple(Q(s) for s in ('0','1/256','1/128','1/64','1/32','1/16','1/8','1/4',
    '3/8','1/2','5/8','3/4','7/8','15/16','31/32','63/64','127/128','255/256','1'))


def bernstein_to_power(B):
    """Expand C(3,i)t^i(1-t)^(3-i) in each coordinate."""
    P = [[Q(0) for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            for k in range(i,4):
                for l in range(j,4):
                    P[k][l] += B[i][j]*comb(3,i)*comb(3-i,k-i)*(-1)**(k-i)*comb(3,j)*comb(3-j,l-j)*(-1)**(l-j)
    return P


def affine_substitution(P, a, b, c, d):
    """Substitute old coordinates t=a+(b-a)x, s=c+(d-c)y."""
    out = [[Q(0) for _ in range(4)] for _ in range(4)]
    for k in range(4):
        for l in range(4):
            for i in range(k+1):
                for j in range(l+1):
                    out[i][j] += P[k][l]*comb(k,i)*a**(k-i)*(b-a)**i*comb(l,j)*c**(l-j)*(d-c)**j
    return out


def power_to_bernstein(P):
    """Use t^k=sum_{i>=k} C(i,k)/C(3,k) B_i^3(t)."""
    return [[sum((P[k][l]*Q(comb(i,k),comb(3,k))*Q(comb(j,l),comb(3,l))
        for k in range(i+1) for l in range(j+1)),Q(0)) for j in range(4)] for i in range(4)]


def restrict(B, domain):
    return power_to_bernstein(affine_substitution(bernstein_to_power(B),*domain))


def matrix(rows):
    assert len(rows)==4 and all(len(row)==4 for row in rows)
    return [list(map(Q,row)) for row in rows]


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
    assert run['source_sha256_before'] == run['source_sha256_after'] == j['source_sha256'] == digest(ROOT/'extend_piecewise.py') == SOURCE_SHA
    assert run['log_sha256'] == digest(ROOT/'extend.log')
    assert j['cover_sha256'] == digest(ROOT/'leaf_cover.json')
    for inp in j['inputs']:
        assert digest(ROOT.parent/inp['file']) == inp['sha256']
    base_path, old_path, sup_path = (ROOT.parent/s for s in
        ('round01/scaled_certificate.json','round22/leaf_cover.json','round19/new_model.json'))
    assert digest(base_path) == BASE_SHA and digest(old_path) == OLD_SHA and digest(sup_path) == SUP_SHA
    base, old, supplemental = read(base_path), read(old_path), read(sup_path)
    c1 = Q(base['c'])
    assert (Q(base['alpha']),Q(base['R']),Q(base['h']),base['N'])==(Q(47,40),Q(33,20),Q(1,32),80)
    assert c1 > 0 and cover['baseline_sha256'] == BASE_SHA and cover['parents'] == [[52,1]]
    assert tuple(map(Q,cover['grid'])) == GRID and cover['roots'] == old['roots'] and len(cover['roots']) == 324
    assert cover['source_cover_sha256'] == OLD_SHA and cover['source_Psi615_sha256'] == SUP_SHA
    assert cover['source_pieces663_sha256']==PIECES_SHA
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
            if i not in old_pending and i != 663:
                assert n['status'] == before['status']
                assert n.get('children') == before.get('children')
                if 'model_file' in before:
                    assert n['model_sha256'] == before['model_sha256']
                    assert n['model_file'] == normpath('../round22/'+before['model_file'])
            if 'partial_model_file' in before:
                assert n['historical_partial_optional'] is True
                assert n['partial_model_sha256'] == before['partial_model_sha256']
                # Historical partial files are metadata only, never read or required.
                assert n['partial_model_file'] == normpath('../round22/'+before['partial_model_file'])
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

    model_cache, rows_cache, piece_cache, modes = {}, {}, {}, Counter()
    piece_slabs = 0
    piece_gains = []
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
    def verify_pieces(m, n, pieces):
        """Different restriction algorithm: saved lower Bernstein -> power -> substitute -> Bernstein."""
        nonlocal piece_slabs
        assert len(pieces) == 4
        psi=m['Psi'];D=signed_deduction(psi,n,'mixed_derivatives','one_sided_lower_deduction')
        assert psi['unscaled_by_c'] and psi['Bernstein_degree']==[3,3] and 'remainder_upper' not in psi
        assert set(psi['center_values'])=={'F0','Fu','Fv'}
        for v in psi['center_values'].values():assert Q(v['lower'])<=Q(v['upper'])
        assert len(psi['bernstein'])==len(m['AS']['bernstein'])==16
        t1=[[Q(0) for _ in range(4)] for _ in range(4)]
        AS=[[Q(0) for _ in range(4)] for _ in range(4)]
        for k,(pb,sb) in enumerate(zip(psi['bernstein'],m['AS']['bernstein'])):
            i,jj=divmod(k,4)
            assert pb['index']==sb['index']==[i,jj] and Q(pb['lower'])<=Q(pb['upper']) and Q(sb['lower'])<=Q(sb['upper'])
            t1[i][jj]=Q(pb['lower']);AS[i][jj]=Q(sb['lower'])
        es=Q(m['AS']['remainder_upper']);assert es>=0
        muv=Q(psi['Muv_upper']);assert muv>=0
        u0,u1,v0,v1=map(Q,n['physical']);x0,x1,y0,y1=map(Q,n['local'])
        uc,vc=(u0+u1)/2,(v0+v1)/2
        halves=((Q(0),Q(1,2)),(Q(1,2),Q(1)))
        rectangles=[];out=[]
        for pid,(iu,iv) in enumerate(((0,0),(1,0),(0,1),(1,1))):
            guard();piece=pieces[pid]
            assert piece['piece_id']==pid and piece['proof_piece'] is True and piece['source_tree_depth']==n['depth']
            if 'source_tree_node_id' in piece:assert piece['source_tree_node_id']==n['id']
            assert piece['same_sign']==(iu==iv)
            assert (piece['u_sign'],piece['v_sign'])==('-' if iu==0 else '+','-' if iv==0 else '+')
            xa,xb=halves[iu];ya,yb=halves[iv];domain=(xa,xb,ya,yb)
            a,b=u0+(u1-u0)*xa,u0+(u1-u0)*xb;c,d=v0+(v1-v0)*ya,v0+(v1-v0)*yb
            assert tuple(map(Q,piece['physical']))==(a,b,c,d)
            assert tuple(map(Q,piece['local']))==(x0+(x1-x0)*xa,x0+(x1-x0)*xb,y0+(y1-y0)*ya,y0+(y1-y0)*yb)
            rectangles.append((a,b,c,d))
            tt,ss=restrict(t1,domain),restrict(AS,domain)
            assert matrix(piece['restricted_T1_lower'])==tt and matrix(piece['restricted_AS_lower'])==ss
            cp=[[Q(0) for _ in range(4)] for _ in range(4)]
            if iu!=iv:
                cp[0][0]=muv*(a-uc)*(c-vc);cp[1][0]=muv*(b-a)*(c-vc)
                cp[0][1]=muv*(a-uc)*(d-c);cp[1][1]=muv*(b-a)*(d-c)
            corr=power_to_bernstein(cp)
            assert matrix(piece['bilinear_Bernstein_degree3'])==corr
            corners=[[sum((cp[k][l]*Q(i)**k*Q(jj)**l for k in range(4) for l in range(4)),Q(0)) for jj in range(2)] for i in range(2)]
            assert [list(map(Q,line)) for line in piece['bilinear_Bernstein_degree1']]==corners
            assert len(piece['rows'])==16
            pr=[]
            for k,row in enumerate(piece['rows']):
                guard();i,jj=divmod(k,4);assert row['bernstein_index']==[i,jj]
                cc=tt[i][jj]+corr[i][jj];sc=ss[i][jj]-es;rho=c1*cc+sc
                oldc=tt[i][jj]-D;oldrho=c1*oldc+sc;gain=rho-oldrho
                assert gain==c1*(D+corr[i][jj]) and gain>=0
                assert tuple(Q(row[key]) for key in ('c_column','s_column','baseline_slack','gain'))==(cc,sc,rho,gain)
                if 'restricted_uniform_slack' in row:assert Q(row['restricted_uniform_slack'])==oldrho
                else:
                    assert tuple(Q(row[key]) for key in ('restricted_old_c_column','restricted_old_s_column','restricted_old_slack'))==(oldc,sc,oldrho)
                pr.append((cc,sc,rho));piece_gains.append(gain)
            assert piece['accepted']==all(rho>0 for _,_,rho in pr)
            out.extend(pr)
        # A second partition layer: 4 proof rectangles cover this one ORIGINAL tree leaf.
        xs=sorted({u0,u1}|{x for a,b,c,d in rectangles for x in (a,b)})
        for left,right in zip(xs,xs[1:]):
            cursor=v0
            for bottom,top in sorted((c,d) for a,b,c,d in rectangles if a<=left and right<=b):
                assert bottom==cursor;cursor=top
            assert cursor==v1
        assert sum(((b-a)*(d-c) for a,b,c,d in rectangles),Q(0))==(u1-u0)*(v1-v0)
        piece_slabs+=len(xs)-1
        return out
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
        elif m['mode']=='outside_signed':
            assert m['fixed_hlcert_sha256']==FIXED_SHA
            a=signed_coefficients(m['Psi'],n)
        else:
            assert m['mode']=='outside_quadrants' and m['fixed_hlcert_sha256']==FIXED_SHA
            assert m['proof_restriction_source_sha256']==ALGEBRA_SHA and m['proof_region_count']==4 and 'rows' not in m
            reconstructed=verify_pieces(m,n,m['proof_pieces'])
            piece_cache[n['id']]=m['proof_pieces']
        if m['mode']!='outside_quadrants':reconstructed=match_rows(m['rows'],a,b)
        assert m['accepted'] == all(rho > 0 for _,_,rho in reconstructed)
        if m['mode'] == 'const':
            assert m['accepted'] == (n['id'] != 615 and old_nodes[n['id']]['status']=='accepted')
        elif m['mode'] == 'outside_poly':
            assert m['accepted'] == (n['id']!=661 and old_nodes[n['id']]['status']=='accepted')
        elif m['mode']=='outside_signed':
            assert m['accepted']==(old_nodes[n['id']]['status']=='accepted')
        else:
            assert m['accepted']==(n['status']=='accepted')
        model_cache[n['id']] = m
        rows_cache[n['id']] = reconstructed
        modes[m['mode']] += 1
        rebuilt_rows += len(reconstructed)
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
    assert signed['source_cover_sha256'] == ORIGIN20_SHA and signed['fixed_hlcert_sha256'] == FIXED_SHA
    assert signed['baseline_sha256'] == BASE_SHA
    assert signed['local'] == n661['local'] and signed['physical'] == n661['physical']
    assert signed['reused_T1_Bernstein'] == model_cache[661]['Psi']['bernstein']
    deduction = signed_deduction(signed,n661,'derivatives','Psi_one_sided_lower_deduction')
    a = [Q(v['lower'])-deduction for v in model_cache[661]['Psi']['bernstein']]
    rows_cache[661] = match_rows(signed['rows'],a,bern(model_cache[661]['AS']))
    assert all(rho > 0 for _,_,rho in rows_cache[661])
    rebuilt_rows += 16
    n663=nodes[663]
    assert n663['status']=='accepted' and n663['mode']=='outside_quadrants_reused' and n663['previous_status']=='unresolved_depth4'
    p663=ROOT/n663['supplemental_piece_file']
    assert digest(p663)==n663['supplemental_piece_sha256']==PIECES_SHA
    proof663=read(p663)
    assert proof663['source_model_sha256']==n663['model_sha256']
    assert proof663['source_cover_sha256']==OLD_SHA and proof663['baseline_sha256']==BASE_SHA
    assert proof663['local']==n663['local'] and proof663['physical']==n663['physical']
    assert proof663['source_tree_depth']==4 and proof663['source_tree_changed'] is False
    assert proof663['proof_partition']=='four_center_sign_quadrants'
    for key,source in (('source_T1_lower',model_cache[663]['Psi']),('source_AS_lower',model_cache[663]['AS'])):
        assert matrix(proof663[key])==[[Q(source['bernstein'][4*i+jj]['lower']) for jj in range(4)] for i in range(4)]
    assert Q(proof663['AS_remainder_upper'])==Q(model_cache[663]['AS']['remainder_upper'])
    assert Q(proof663['Muv_upper'])==Q(model_cache[663]['Psi']['Muv_upper'])
    assert Q(proof663['old_uniform_D'])==Q(model_cache[663]['Psi']['one_sided_lower_deduction'])
    rows_cache[663]=verify_pieces(model_cache[663],n663,proof663['pieces'])
    piece_cache[663]=proof663['pieces'];rebuilt_rows+=64
    assert all(piece['accepted'] for piece in proof663['pieces']) and all(rho>0 for _,_,rho in rows_cache[663])
    leaf_ids = [leaf['node_id'] for leaf in cover['leaves']]
    assert len(leaf_ids) == len(set(leaf_ids)) and set(leaf_ids) == {n['id'] for n in nodes if n['status']=='accepted'}
    for leaf in cover['leaves']:
        guard()
        n = nodes[leaf['node_id']]
        for field in ('parent','local','physical','depth','mode','model_file','model_sha256'):
            assert leaf[field] == n[field]
        expected = rows_cache[n['id']]
        assert all(rho > 0 for _,_,rho in expected)
        if n['id'] in piece_cache:
            assert leaf['proof_region_count']==4 and n['proof_region_count']==4 and 'rows' not in leaf
            assert len(leaf['proof_pieces'])==4
            for saved,wanted in zip(leaf['proof_pieces'],piece_cache[n['id']]):
                for key,value in wanted.items():assert saved[key]==value
                assert saved['source_tree_node_id']==n['id']
            actual=[(Q(v['c_column']),Q(v['s_column']),Q(v['baseline_slack'])) for piece in leaf['proof_pieces'] for v in piece['rows']]
            assert actual==expected
            if n['id']==663:
                assert leaf['supplemental_piece_file']==n['supplemental_piece_file'] and leaf['supplemental_piece_sha256']==PIECES_SHA
            else:
                assert leaf['proof_file']==n['model_file'] and leaf['proof_sha256']==n['model_sha256'] and leaf['origin']=='round24'
        else:
            assert leaf['proof_region_count']==1 and 'proof_pieces' not in leaf
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
    new_models = {i for i,m in model_cache.items() if m['mode']=='outside_quadrants'}
    assert new_models == set(processing_order)
    assert j['new_boxes']['started']==j['new_boxes']['completed']==len(new_models)
    assert j['new_boxes']['failed']==j['new_boxes']['interrupted']==0
    assert j['new_boxes']['accepted']==sum(nodes[i]['status']=='accepted' for i in new_models)
    assert j['new_boxes']['splits']==sum(nodes[i]['status']=='split' for i in new_models)
    assert j['AS_models_saved']==j['Psi_models_saved']==len(new_models)
    assert len(new_models)<=1500 and by_kind['Hint']<=12000 and by_kind['acb_integral']<=12000
    assert all(j[k]==0 for k in ('full_absolute_Psi_calls','new_const_calls','new_T_calls','new_point_calls','LP_calls','full_checker_calls'))
    for kind,number,cap in (('pieces',4*len(new_models),6000),('rows',64*len(new_models),96000)):
        assert j['proof_counts'][kind]==dict(started=number,completed=number,failed=0,interrupted=0) and number<=cap
    assert j['proof_in_flight'] is None
    regions=sum(4 if i in piece_cache else 1 for i in leaf_ids)
    assert cover['accepted_tree_leaves']==j['accepted_tree_leaves']==len(leaf_ids)
    assert cover['accepted_proof_regions']==j['accepted_proof_regions']==regions
    assert j['limits']==dict(seconds_hard=300,Arb_bits=192,new_started_boxes=1500,max_original_root_depth=4,AS=1500,Psi=1500,Hint=12000,acb_integral=12000,proof_pieces=6000,new_rows=96000)
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
        scope='Exact saved enclosures/rows, three supplemental origins, ledger, tree-frontier partition and internal proof partitions; no model or integral reevaluation.',
        complete=cover['complete'],old_const_models_checked=modes['const'],old_absolute_poly_models_checked=modes['outside_poly'],old_signed_models_checked=modes['outside_signed'],new_piecewise_models_checked=modes['outside_quadrants'],
        supplemental_models_checked=3,rows_reconstructed=rebuilt_rows,accepted_leaves=len(leaf_ids),
        accepted_proof_regions=regions,accepted_rows_checked=16*regions,proof_piece_slabs_checked=piece_slabs,
        proof_pieces_reconstructed=4*len(piece_cache),
        restriction_algorithm='Bernstein to power, exact affine substitution, power to Bernstein; no de Casteljau calls.',accepted_area=str(accepted_area),
        pending_area=str(area(cover['pending'])),unresolved_area=str(area(cover['unresolved'])),
        frontier_nodes=len(frontier),frontier_area='1',slabs_checked=len(xs)-1,
        status_counts=dict(Counter(n['status'] for n in nodes)),call_counts=dict(by_kind),trace_events=events,
        new_boxes_crossing_coordinate_R=cross_R,
        minimum_accepted_slack=str(min(rho for i in leaf_ids for _,_,rho in rows_cache[i])),
        worst_new_row=None if worst_new is None else dict(slack=str(worst_new[0]),node=worst_new[1],piece=worst_new[2]//16,index=list(divmod(worst_new[2]%16,4))),
        elapsed_seconds=monotonic()-start,result_sha256=digest(ROOT/'result.json'),
        cover_sha256=digest(ROOT/'leaf_cover.json'),run_sha256=digest(ROOT/'extend.run.json'),
        calls_sha256=j['calls_sha256'],verifier_sha256=digest(Path(__file__)))
    OUT.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output))


if __name__ == '__main__':
    main()
