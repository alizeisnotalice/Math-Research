"""Read-only convex-boundary chain certificate and exact rational heat anchors."""
import argparse,ast,gzip,hashlib,json,math,os,sys,time
from collections import Counter
from fractions import Fraction as F
from itertools import permutations,product
from pathlib import Path

HERE=Path(__file__).resolve().parent;RUN='audit_vertex032_fixed_001'
MANIFEST='9084775b6975369b37c6e20b19778e21be54f75b922f1874d1411fe94037e640'
TABLE='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
V=[(1,)+s for s in product((-1,1),repeat=4)]
PERMS=[(p,(-1)**sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))) for p in permutations(range(5))]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(obj):return json.dumps(obj,sort_keys=True,separators=(',',':')).encode()
def project(p):
    p=p.resolve();return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def psi(z):return sum(abs(dot(v,z)) for v in V)
def det5(rows):return sum(sign*math.prod(rows[i][p[i]] for i in range(5)) for p,sign in PERMS)
def rank(rows):
    basis={}
    for row in rows:
        r=list(map(F,row))
        for p in sorted(basis):
            v=r[p]
            if v:r=[x-v*y for x,y in zip(r,basis[p])]
        p=next((j for j,x in enumerate(r) if x),None)
        if p is not None:
            v=r[p];basis[p]=[x/v for x in r]
            if len(basis)==5:return 5
    return len(basis)
def parse(s):
    inner=s[1:-1]
    if inner.startswith('+/- '):
        r=F(inner[4:]);return -r,r
    c,r=map(F,inner.split(' +/- '));return c-r,c+r
def pi_bounds():
    def atan(q):
        lo=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(80)),F(0))
        return lo,lo+F(1,161*q**161)
    a,b=atan(5),atan(239);return 16*a[0]-4*b[1],16*a[1]-4*b[0]
def heat_rational(Z,simplices,x,t):
    dots=[sum((F(a)*b for a,b in zip(x,z)),F(0)) for z in Z]
    scale=math.lcm((144*t).denominator,*(d.denominator for d in dots))
    a=int(scale*144*t);assert a>0
    bs=[int(scale*d) for d in dots];real=imag=F(0)
    for s in simplices:
        P,Q=1,0
        for j in s['vertices']:P,Q=P*a+Q*bs[j],Q*a-P*bs[j]
        denominator=P*P+Q*Q;assert denominator>0
        weight=s['det_numerator']*scale**5
        real+=F(weight*P,denominator);imag-=F(weight*Q,denominator)
    assert imag==0 and real>0
    return real,imag

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',required=True,type=Path);parser.add_argument('--model-root',required=True,type=Path)
    parser.add_argument('--run',action='store_true');args=parser.parse_args()
    source,model=project(args.source_root),project(args.model_root);author=source/'astra_vertex017'
    manifest=author/'manifest.json';assert sha(manifest)==MANIFEST;info=json.loads(manifest.read_text());locked={'astra_vertex017/manifest.json':MANIFEST}
    for row in info['white_list']:
        p=author/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes'];locked['astra_vertex017/'+row['path']]=sha(p)
    assert sha(author/'boundary_triangulation.json.gz')==TABLE
    anchors=json.loads((author/'anchors_result.json').read_text())
    for legacy,expected in anchors['source_hashes'].items():
        suffix=Path(legacy).parts[-2:];choices=[r.joinpath(*suffix) for r in (source,model)]
        p=next((p for p in choices if p.exists() and sha(p)==expected),None);assert p is not None;locked['/'.join(suffix)]=sha(p)
    contacts=json.loads((author/'contacts_result.json').read_text());classification={}
    for name in ('005_vertex_rays_exact.py','005_vertex_rays_exact_results.json'):
        expected=contacts['source_hashes'][name];choices=[r/'audit'/name for r in (source,model)]
        p=next((p for p in choices if p.exists() and sha(p)==expected),None);assert p is not None
        locked['audit/'+name]=sha(p);classification[name]=p
    if not args.run:
        print(json.dumps({'run_handle':RUN,'status':'locked_inputs_only','calculation_executed':False}));return
    pid=os.getpid();print(json.dumps({'run_handle':RUN,'pid':pid,'status':'started'}),flush=True);began=time.monotonic()
    table=json.loads(gzip.decompress((author/'boundary_triangulation.json.gz').read_bytes()))
    payload={k:table[k] for k in ('denominator','vertices','facets','simplices')}
    payload_hash=hashlib.sha256(canonical(payload)).hexdigest();assert payload_hash==info['geometry_payload_sha256']
    Z=[tuple(z) for z in table['vertices']];facets=table['facets'];simplices=table['simplices']
    assert table['denominator']==144 and len(Z)==280 and len(facets)==1882 and len(simplices)==3504
    assert Z==sorted(set(Z)) and all(len(z)==5 and psi(z)==144 for z in Z)
    assert all(tuple(-x for x in z) in set(Z) for z in Z) and rank(Z)==5
    original=json.loads(classification['005_vertex_rays_exact_results.json'].read_text())
    true_records=[r for r in original['records'] if r['normalized_ray_is_true_vertex']]
    assert [(tuple(r['primitive_ray']),r['q_equals_psi_ray']) for r in true_records]==[((1,1,0,0,0),16),((1,1,1,1,0),24),((2,1,1,1,1),36)]
    expanded=set();orbit_counts=[]
    for record in true_records:
        ray=tuple(record['primitive_ray']);q=record['q_equals_psi_ray'];assert psi(ray)==q and 144%q==0
        images={tuple(144//q*a*s for a,s in zip(p,signs)) for p in set(permutations(ray)) for signs in product((-1,1),repeat=5)}
        orbit_counts.append(len(images));expanded|=images
    assert expanded==set(Z) and orbit_counts==[40,80,160]
    incidence=[[] for _ in Z];owner=Counter();normals=set();contact_sets=set();witness=None
    for fi,facet in enumerate(facets):
        a=tuple(facet['normal']);assert len(a)==5 and a not in normals;normals.add(a)
        values=[dot(a,z) for z in Z];assert max(values)==144
        contact=tuple(i for i,v in enumerate(values) if v==144)
        assert contact==tuple(facet['vertices']) and contact not in contact_sets;contact_sets.add(contact)
        assert rank([Z[i] for i in contact])==5
        center_sum=tuple(sum(Z[i][j] for i in contact) for j in range(5))
        dots=[dot(v,center_sum) for v in V];assert all(dots)
        signed_sum=tuple(sum((1 if d>0 else -1)*v[j] for v,d in zip(V,dots)) for j in range(5))
        assert signed_sum==a # exact certificate that the halfspace contains {psi<=1}
        for i in contact:incidence[i].append(a)
        for si in facet['simplices']:
            assert 0<=si<len(simplices);owner[si]+=1
            assert set(simplices[si]['vertices'])<=set(contact)
        if witness is None and len(contact)==5 and len(facet['simplices'])==1:
            si=facet['simplices'][0];assert tuple(simplices[si]['vertices'])==contact
            witness={'facet_index':fi,'normal':a,'vertices':contact,'simplex_index':si}
        if time.monotonic()-began>120:raise TimeoutError('finite audit geometry checkpoint')
    assert owner==Counter({i:1 for i in range(3504)}) and witness is not None
    assert all(rank(rows)==5 for rows in incidence)
    ridges=Counter();oriented=Counter();simplex_sets=set();total_det=0
    for s in simplices:
        ids=tuple(s['vertices']);assert len(ids)==5 and ids==tuple(sorted(set(ids))) and ids not in simplex_sets;simplex_sets.add(ids)
        D=det5([Z[i] for i in ids]);assert D!=0 and abs(D)==s['det_numerator'];total_det+=abs(D)
        for j in range(5):
            ridge=ids[:j]+ids[j+1:];ridges[ridge]+=1;oriented[ridge]+=(-1)**j*(1 if D>0 else -1)
    assert len(ridges)==8760 and set(ridges.values())=={2} and all(v==0 for v in oriented.values())
    mass=F(total_det,144**5);assert mass==F(35509,23887872)
    # Parse immutable old BASE_CONES only; do not import or expand 26880 terms.
    old_source=model/'numerics/body_fan_geometry.py';assert sha(old_source)==locked['numerics/body_fan_geometry.py']
    node=next(n for n in ast.parse(old_source.read_text()).body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='BASE_CONES' for t in n.targets))
    base=ast.literal_eval(node.value);old_sum=3840*sum((F(abs(det5(rows)),math.prod(psi(r) for r in rows)) for rows in base),F(0));assert old_sum==mass
    ray_types={tuple(x//math.gcd(*r) for x in r) for rows in base for r in rows}
    assert ray_types=={tuple(row['primitive_ray']) for row in original['records']}
    pil,pih=pi_bounds();assert F(314159,100000)<pil<pih<F(22,7)
    anchor_results=[]
    for row in anchors['records']:
        x=tuple(map(F,row['point']));t=F(row['time']);real,imag=heat_rational(Z,simplices,x,t)
        lo,hi=real/(32*pih**5),real/(32*pil**5);scale=10**95
        lo,hi=F((lo*scale).__floor__(),scale),F((hi*scale).__ceil__(),scale)
        assert 0<lo<hi
        for label in ('new3504_real','old26880_real'):
            lower,upper=parse(row[label]);assert lower<=lo<hi<=upper
        assert max(parse(row['new3504_real'])[0],parse(row['old26880_real'])[0])<=min(parse(row['new3504_real'])[1],parse(row['old26880_real'])[1])
        assert parse(row['new3504_imag'])[0]<=0<=parse(row['new3504_imag'])[1] and parse(row['new_minus_old'])[0]<=0<=parse(row['new_minus_old'])[1]
        anchor_results.append({'point':[str(v) for v in x],'time':str(t),'exact_rational_imaginary':str(imag),
                               'rational_pre_pi_sum_sha256':hashlib.sha256(str(real).encode()).hexdigest(),
                               'rational_pre_pi_numerator_bits':real.numerator.bit_length(),
                               'heat_rational_interval':[str(lo),str(hi)],'inside_saved_new_and_old_balls':True})
    assert len(anchor_results)==2
    out={'run_handle':RUN,'pid':pid,'status':'pass_original_vertex_fan_coverage_and_two_exact_rational_anchors',
         'model':'alpha=1 Cauchy psi=sum_16 |v.xi|; not the alpha=3/2 Green model',
         'geometry':{'vertices':280,'facet_support_and_sign_certificates':1882,'nondegenerate_boundary_simplices':3504,
                      'orbits':orbit_counts,'all_280_active_normal_rank_five':True,'ridges':len(ridges),'ridge_incidence':2,
                      'all_outward_boundary_coefficients_zero':True,'degree_one_simplicial_facet_witness':witness,
                      'coverage_proof':'positive closed boundary chain on convex K has constant multiplicity; simplicial facet witness fixes it to one; facet sign certificates give K={psi<=1}',
                      'normalized_determinant_sum':str(mass),'volume_P':str(mass/120),'old_seven_base_sum_matches':True,
                      'p1_zero_times_pi5':str(mass/32),'geometry_payload_sha256':payload_hash},
         'anchors':anchor_results,
         'scope_limit':'Independent full coverage and relative-interior nonoverlap for cone integration. Full face lattice/f-vector and stronger face-to-face pulling restrictions are locked production evidence, not regenerated. No global g/I/C5 claim.',
         'arithmetic':'Python Fraction/integers only, 80-term Machin pi intervals; no Arb/SymPy or old heat evaluator executed',
         'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'registration_sha256':sha(HERE/'registration.md'),
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began,'contacts_or_faces_or_pulling_replayed':False}}
    (HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('run_handle','pid','status','geometry','runtime')},indent=2))

if __name__=='__main__':main()
