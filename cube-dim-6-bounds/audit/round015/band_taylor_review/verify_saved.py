#!/usr/bin/env python3
"""Independent endpoint audit: eight frozen calibrations and two geometry toys.

Never imports or executes the producer, a heat kernel, or a gradient routine.
Old 384-bit balls are converted to rational intervals; Phi is independently
evaluated as (u/12 + sqrt(1+(u/12)^2))^12 at 768 bits.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product, combinations
import hashlib, json, os, platform, time
import sympy as sp
import flint
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
TOP = Path('/Users/zhengzhihao/Desktop/T')
A = TOP/'Math-Research-cube6-theory/cube-dim-6-bounds'
B = TOP/'Math-Research-cube6-upper-astra/cube-dim-6-bounds'
SOURCE = A/'upper/round015/band_taylor_cover'
INPUTS = {}

def read(path):
    data = path.read_bytes()
    INPUTS[str(path)] = hashlib.sha256(data).hexdigest()
    return json.loads(data)

def endpoint(x, up=False):
    v = x.upper() if up else x.lower()
    assert v.is_exact()
    r = v.fmpq()
    return Q(int(r.numerator), int(r.denominator))

def interval(s):
    ctx.prec = 384
    x = arb(s)
    assert x.is_finite()
    return endpoint(x), endpoint(x, True)

def overlap(x, y):
    return max(x[0], y[0]) <= min(x[1], y[1])

def subtract(x, y):
    return x[0]-y[1], x[1]-y[0]

def aq(x):
    return arb(x.numerator)/x.denominator

def dot_lower(grad, h):
    return sum((min(lo*x, hi*x) for (lo, hi), x in zip(grad, h)), Q(0))

def norm2(v):
    return sum((x*x for x in v), Q(0))

def dot(v, w):
    return sum((a*b for a, b in zip(v, w)), Q(0))

def plane_distance2(x, v, w):
    vv, ww, vw = dot(v,v), dot(w,w), dot(v,w)
    xv, xw = dot(x,v), dot(x,w)
    determinant = vv*ww-vw*vw
    assert determinant > 0
    return norm2(x)-(ww*xv*xv-2*vw*xv*xw+vv*xw*xw)/determinant

def geometry_toys():
    root = [tuple(map(Q, [1]+[1]*j+[0]*(5-j))) for j in range(6)]
    bary = tuple(sum((v[i] for v in root), Q(0))/6 for i in range(6))
    bodies = [(Q(1),)+tuple(map(Q,s)) for s in product((-1,1), repeat=5)]
    pairs = list(combinations(bodies, 2))
    assert len(pairs) == 496
    # This direct all-pair calculation does not assume A012's closest-plane theorem.
    assert all(min(plane_distance2(x,v,w) for v,w in pairs)==0 for x in root)
    assert min(plane_distance2(bary,v,w) for v,w in pairs)==Q(1,9)
    delta = Q(1,100000)
    canonical = [(tuple(map(Q,[1]*k+[0]*(6-k))),
                  tuple(map(Q,[0]*k+[1]*(6-k)))) for k in range(1,6)]
    assert not any(all(plane_distance2(x,v,w)<=delta**2 for x in root)
                   for v,w in canonical)
    tube = Q(1,5000000000)*(1/delta-40000)
    assert tube==Q(3,250000)>Q(6,125)/4096
    # Toy 2: exact union via barycentric identities, with volumes as a crosscheck.
    lengths = {ij:norm2(tuple(a-b for a,b in zip(root[ij[0]],root[ij[1]])))
               for ij in combinations(range(6),2)}
    edge = min(ij for ij,d in lengths.items() if d==max(lengths.values()))
    assert edge==(0,5)
    mid = tuple((a+b)/2 for a,b in zip(root[0],root[5]))
    left, right = root.copy(), root.copy()
    left[5], right[0] = mid, mid
    def volume(vertices):
        matrix = sp.Matrix([[sp.Rational(vertices[j][i]-vertices[0][i])
                             for j in range(1,6)] for i in range(1,6)])
        return abs(matrix.det())/120
    assert volume(root)==sp.Rational(1,120)
    assert volume(left)==volume(right)==sp.Rational(1,240)
    lam = sp.symbols('l0:6')
    for h in range(6):
        parent = sum(lam[i]*root[i][h] for i in range(6))
        common = sum(lam[i]*root[i][h] for i in range(1,5))
        assert sp.expand(parent-(lam[0]-lam[5])*root[0][h]-2*lam[5]*mid[h]-common)==0
        assert sp.expand(parent-(lam[5]-lam[0])*root[5][h]-2*lam[0]*mid[h]-common)==0
    return {'count':2,'root_vertices_in_S2':True,'all_pair_count':496,
            'root_barycenter_distance_squared':'1/9','same_plane_root_rejected':True,
            'T48_local_tube_lower':str(tube),'longest_edge':[0,5],
            'root_volume':'1/120','child_volumes':['1/240','1/240'],
            'exact_barycentric_split_identities':True}

def main():
    start = time.perf_counter()
    meta = read(SOURCE/'input_manifest.json')
    saved = read(SOURCE/'results.json')
    post = read(SOURCE/'postcheck.json')
    for path, expected in meta['input_sha256'].items():
        p = Path(path)
        actual = hashlib.sha256(p.read_bytes()).hexdigest()
        assert actual==expected, path
        INPUTS[path] = actual
    for name in ('calibrate.py','report.md'):
        p = SOURCE/name
        INPUTS[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    assert INPUTS[str(SOURCE/'calibrate.py')]==meta['source_sha256']==saved['source_sha256']
    assert INPUTS[str(SOURCE/'results.json')]==post['results_sha256']
    old = read(A/'upper/round007/g_time_moments/results.json')
    enddata = read(B/'upper/round010/short_minimum_diagnostic/gradient_results.json')
    reference = read(B/'upper/round009/short_potential/postprocess_results.json')
    refutations = read(B/'upper/round010/short_minimum_diagnostic/checked_results.json')
    b15 = meta['B015_context']
    read(Path(b15['path']))
    assert INPUTS[b15['path']]==b15['sha256']
    assert b15['strictly_below_threshold64'] is True
    assert old['center_numerators']==enddata['center_numerators']
    assert old['denominator']==enddata['denominator']
    c = tuple(Q(n,old['denominator']) for n in old['center_numerators'])
    assert c==tuple(map(Q,saved['center']))
    T, target, bare = Q(4,125), Q(1,125000), Q(1,128000)
    assert str(T) in refutations['strict_refutation_times']
    H = {Q(row['tau']):interval(row['H_arb384']) for row in old['time_nodes']}
    ga = {Q(row['tau']):row for row in old['gradients']}
    end = next(row for row in enddata['gradient_passes'] if Q(row['tau'])==T)
    HT = interval(end['H_arb384'])
    assert overlap(HT,H[T])
    gT = list(map(interval,end['gradient_x2_to_x6']))
    oldrow = next(r for r in reference['rows'] if Q(r['T'])==T)
    oldradius = Q(oldrow['local_band_radius_m0'])
    assert oldradius==Q(991,500000)
    assert Q(oldrow['m0_for_cost_125_over_2'])==target
    expected = {(a,r) for a in map(Q,('1/1000','1/500','1/250','1/125'))
                for r in map(Q,('1/250','1/125'))}
    assert len(saved['rows'])==8
    assert {(Q(r['a']),Q(r['radius'])) for r in saved['rows']}==expected
    cases = []
    for row in saved['rows']:
        a, radius = Q(row['a']),Q(row['radius'])
        assert Q(row['T'])==T and row['infinite_time_tail_added'] is False
        assert overlap(interval(ga[a]['H_arb384']),H[a])
        grad = [subtract(interval(x),y) for x,y in zip(ga[a]['gradient_x2_to_x6'],gT)]
        B_interval = subtract(H[a],HT)
        LB = Q(row['B_center_lower_exact'])
        assert LB<=B_interval[0]
        bounds = [(max(Q(0),x-radius), min(Q(1),x+radius)) for x in c[1:]]
        assert bounds==[tuple(map(Q,v)) for v in row['cell_free_bounds']]
        vertices = [(Q(1),)+v for v in product(*bounds)]
        assert len(vertices)==row['cell_vertices']==32
        assert all(all(v[i]>=v[i+1] for i in range(5)) and v[-1]>=0 for v in vertices)
        weights = []
        for bits in product((0,1),repeat=5):
            weight = Q(1)
            for x,(lo,hi),bit in zip(c[1:],bounds,bits):
                w = (x-lo)/(hi-lo)
                weight *= w if bit else 1-w
            weights.append(weight)
        assert all(w>=0 for w in weights) and sum(weights)==1
        assert all(sum(w*v[i] for w,v in zip(weights,vertices))==c[i] for i in range(6))
        hs = [tuple(x-y for x,y in zip(v,c)) for v in vertices]
        assert max(max(map(abs,h)) for h in hs)==radius
        D = max(map(norm2,hs))
        assert D==Q(row['max_displacement_norm_squared'])==5*radius**2
        linear = min(dot_lower(grad,h[1:]) for h in hs)
        ell = Q(row['min_vertex_linear_lower_exact'])
        assert ell<=linear
        assert tuple(map(Q,row['min_linear_vertex'])) in vertices
        grid = sorted(t for t in H if a<=t<=T)
        assert [(Q(z['left']),Q(z['right'])) for z in row['time_bands']]==list(zip(grid,grid[1:]))
        assert grid[0]==a and grid[-1]==T
        ctx.prec = 768
        U = arb(0)
        for band in row['time_bands']:
            t,u = Q(band['left']),Q(band['right'])
            delta = subtract(H[t],H[u])
            upper,lower = Q(band['mass_upper']),Q(band['mass_lower'])
            assert upper>=delta[1] and 0<=lower<=max(Q(0),delta[0])
            assert upper>=lower
            z = aq(radius/t/12)
            phi = (z+(1+z*z).sqrt())**12
            U += aq(upper/(t*t))*phi
        Uupper = Q(row['J2_cell_upper_exact'])
        assert Uupper>=endpoint(U,True)
        loss = D*Uupper/32
        assert loss==Q(row['Taylor_curvature_loss_upper_exact'])
        L = LB+ell-loss
        assert L==Q(row['K_cell_lower_exact'])
        assert (L>target)==row['passes_matched_m0']
        assert (L>bare)==row['passes_bare64_threshold']
        assert Q(row['radius_over_old_Phi_radius_exact'])==radius/oldradius
        cases.append({'a':str(a),'r':str(radius),'strict_matched_pass':L>target,
                      'strict_bare64_pass':L>bare,'all_directed_endpoint_checks':True,
                      'saved_lower_exact':str(L),'display_only':format(float(L),'.14g')})
    assert sum(r['strict_matched_pass'] for r in cases)==4
    assert sum(r['strict_bare64_pass'] for r in cases)==4
    toys = geometry_toys()
    assert saved['global_cover_executed'] is False
    assert saved['unresolved_root'] is True
    assert saved['whole_cube_mT_certified'] is False
    out = {'status':'PASS_LOCAL_THEORY_AND_EIGHT_SAVED_CERTIFICATES_ONLY',
           'pid':os.getpid(),'python':platform.python_version(),'python_flint':flint.__version__,
           'sympy':sp.__version__,'input_ball_bits':384,'independent_Phi_bits':768,
           'Phi_formula':'(u/12+sqrt(1+(u/12)^2))^12',
           'new_kernel_calls':0,'new_gradient_calls':0,'MC_calls':0,'MATLAB_calls':0,
           'calibration_cases':8,'geometry_toys':toys,'rows':cases,
           'passing_geometrically_distinct_cells':1,'passing_radius':'1/250',
           'matched_radius_ratio':'2000/991','global_cover_implemented_or_run':False,
           'global_bound_certified':False,'body_T032_and_T048_excluded':True,
           'input_sha256':INPUTS,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'elapsed_seconds':time.perf_counter()-start}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','pid','elapsed_seconds','calibration_cases',
                                      'passing_geometrically_distinct_cells','new_kernel_calls')}))

if __name__=='__main__':
    main()
