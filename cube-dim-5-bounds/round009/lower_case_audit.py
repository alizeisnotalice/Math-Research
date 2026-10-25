"""Integer witness checker for the one new low-resource planar qualification."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import hashlib,json,os,time
P=Path(__file__).resolve().parents[1];start=time.monotonic()
source=P/'lower_round018/K8_n2_results.json'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(source)=='378c4ca2620674940bb9bc0e0a5ff00b3163be7ae4664f5c9da2b7de56432c77'
a=json.loads(source.read_text());assert a['K']==8 and a['n']==2 and a['selected_mass']==16
atoms={};indices=set()
for atom in a['normalized_lattice_atoms']:
    m,n=atom['m'],atom['n'];x,y=5*F(atom['x']),5*F(atom['y'])
    assert x.denominator==y.denominator==1 and x==4*m-3*n and y==3*m+4*n
    assert (m,n) not in indices;indices.add((m,n));atoms[atom['id']]=(int(x),int(y))
assert len(atoms)==225
reps={(x%25,y%25) for x,y in atoms.values()};assert len(reps)==25
assert all((4*x+3*y)%25==0 and (-3*x+4*y)%25==0 for x,y in reps)
rectangles=set();xs={0,25};ys={0,25}
for cell,ids in a['cell_witnesses']:
    c=[5*F(z) for z in cell];assert all(z.denominator==1 for z in c)
    l,r,b,t=map(int,c);assert 0<=l<r<=25 and 0<=b<t<=25
    assert len(ids)==len(set(ids))==2
    for id in ids:
        x,y=atoms[id];assert x-4<=l<r<=x+4 and y-4<=b<t<=y+4
    assert (l,r,b,t) not in rectangles;rectangles.add((l,r,b,t));xs|={l,r};ys|={b,t}
xs=sorted(xs);ys=sorted(ys)
expected={(l,r,b,t) for (l,r),(b,t) in product(zip(xs,xs[1:]),zip(ys,ys[1:]))}
assert rectangles==expected and len(rectangles)==625
assert sum((r-l)*(t-b) for l,r,b,t in rectangles)==625
sl,sh=map(F,a['s_interval']);rl,rh=map(F,a['r_interval'])
assert 0<2*sl-1<2*sh-1 and (2*sl-1)**5<2<(2*sh-1)**5
assert (2*rl)**5<16<(2*rh)**5 and rl-F(4,5)*sh>0
assert rl-F(4,5)*sh==F(a['strict_physical_radius_margin_lower'])
assert os.getpriority(os.PRIO_PROCESS,0)>=10 and time.monotonic()-start<1
out={'status':'PASS_independent_integer_witness_cover_one_case_only','case':[8,2],
     'distinct_atoms':225,'period_representatives':25,'closed_cells':625,
     'scaled_cell_area':625,'scaled_radius':4,'scaling':'x_scaled=5*x_normalized',
     'normalised_period_area':'25','physical_plane_capture_area':'25*s^2',
     'root_endpoint_and_strict_radius_margin_verified':True,
     'remaining_cases':16,'C5_improvement':False,
     'independence':'Checks integer certificate witnesses and Cartesian tiling; no producer bitmask rescan, no union or event optimization.',
     'sources_sha256':{str(p):sha(p) for p in [source,Path(__file__),P/'round009/lower_case_audit_registration.md']},
     'runtime':{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'seconds':time.monotonic()-start}}
(P/'round009/lower_case_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
