"""Independent fixed CI4B geometry/residual/integral/cost verifier, no optimizer.

Uses the locked global N/D data whose derivation is independently checked by
035_representation_check.py. No author geometry/engine/checker is executed.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import argparse,hashlib,json,time
import sympy as sy
from flint import fmpq_mpoly_ctx,fmpq,arb,ctx as actx
HERE=Path(__file__).resolve().parent
def aq(x):
    x=Q(x);return arb(x.numerator)/x.denominator
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def det3(v):
    a,b,c=v
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def volume(tet):return abs(det3([tuple(tet[j][i]-tet[0][i] for i in range(1,4)) for j in (1,2,3)]))/6
def refine(tet,depth):
    level=[tet]
    for _ in range(depth):
        children=[]
        for t in level:
            lengths={ij:sum((t[ij[0]][k]-t[ij[1]][k])**2 for k in range(1,4)) for ij in combinations(range(4),2)}
            top=max(lengths.values());i,j=min(ij for ij,d in lengths.items() if d==top)
            midpoint=tuple((a+b)/2 for a,b in zip(t[i],t[j]));left=list(t);right=list(t);left[i]=midpoint;right[j]=midpoint
            assert volume(tuple(left))+volume(tuple(right))==volume(t)
            children.extend((tuple(left),tuple(right)))
        level=children
    return level

def run(root):
    start=time.monotonic();sources=json.loads((HERE/'035_source_snapshot.json').read_text())
    for name,h in sources.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    raw=json.loads((root/'upper/u11_minimum_global_results.json').read_text())['physical_cases']
    cert=json.loads((root/'upper/continuous_integral_n4/ci4b_certificate.json').read_text())
    actx.prec=256
    R=[[(1,0,0,0),(1,1,0,0),(2,1,1,0),(3,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,1),(3,1,1,1)]]
    roots=[tuple(tuple(Q(v,r[0]) for v in r) for r in rs) for rs in R]
    assert [volume(t) for t in roots]==[Q(1,36),Q(1,12),Q(1,18)]
    # Exact two-way shell barycentric identities, not only volume equality.
    x,y,z=sy.symbols('x y z');t=x+y+z
    ref=[[(1,0,0),(0,1,0),(0,0,1),(Q(1,2),0,0)],[(0,1,0),(0,0,1),(Q(1,2),0,0),(0,Q(1,2),0)],[(0,0,1),(Q(1,2),0,0),(0,Q(1,2),0),(0,0,Q(1,2))]]
    bary=[(x+t-1,y,z,2*(1-t)),(t+x+y-1,z,2*x,2*(1-t-x)),(2*t-1,2*x,2*y,2*(1-t-x-y))]
    for vertices,b in zip(ref,bary):
        assert sy.expand(sum(b)-1)==0
        assert all(sy.expand(sum(b[j]*vertices[j][k] for j in range(4))-v)==0 for k,v in enumerate((x,y,z)))
    cells=[(0,t) for t in refine(roots[0],8)]
    corner=(Q(1),)*4
    for case in (1,2):
        assert corner in roots[case]
        abc=[v for v in roots[case] if v!=corner]
        assert max(sum((v[i]-1)**2 for i in range(1,4)) for v in abc)<=2
        for k in range(28):
            scale=Q(1,2**k)
            outer=[tuple(q+scale*(v-q) for q,v in zip(corner,a)) for a in abc]
            inner=[tuple(q+scale*(v-q)/2 for q,v in zip(corner,a)) for a in abc]
            a,b,c=outer;aa,bb,cc=inner
            for tet in [(a,b,c,aa),(b,c,aa,bb),(c,aa,bb,cc)]:
                cells.extend((case,z) for z in refine(tet,max(0,8-2*k)))
    assert len(cells)==2440 and len(cert['cells'])==2440
    covered=sum((volume(t) for _,t in cells),Q(0));omitted=Q(5,36*2**84)
    assert covered==Q(1,6)-omitted
    assert Q(cert['ordered_covered_volume'])==covered
    CTX=fmpq_mpoly_ctx.get(('v0','v1','v2','v3'),ordering='lex');v=CTX.gens();zero=CTX.constant(0);sigma=sum(v,zero)
    def constant(q):q=Q(q);return CTX.constant(fmpq(q.numerator,q.denominator))
    def poly(dic):return CTX.from_dict({tuple(map(int,e.split(','))):fmpq(c) for e,c in dic.items()})
    models=[(v[0]**3*poly(row['numerator_dict']),poly(row['denominator_dict'])) for row in raw]
    total=arb(0);constraints=0;minimum=None;cell_records=[]
    print(json.dumps({'stage':'035_integral','cells':len(cells),'geometry':'independently regenerated;all dyadic bisections exact','covered_volume':str(covered),'precision_bits':256,'LP_calls':0}),flush=True)
    for i,((case,tet),row) in enumerate(zip(cells,cert['cells'])):
        assert row['case']==case and tuple(tuple(map(Q,p)) for p in row['vertices'])==tet
        assert volume(tet)>0 and all(p[0]==1 for p in tet)
        pp=list(map(Q,row['P']));assert len(pp)==4 and all(p>=0 for p in pp)
        sub=[sum((constant(tet[j][k])*v[j] for j in range(4)),zero) for k in range(4)]
        numerator,denominator=models[case]
        n=numerator.compose(*sub);d=denominator.compose(*sub)
        dc=d.to_dict();assert dc and all(c>=0 for c in dc.values())
        affine=sum((constant(p)*u for p,u in zip(pp,v)),zero)
        residual=affine*d-n*sigma;coeff=residual.to_dict();assert all(c>=0 for c in coeff.values()),(i,'negative residual')
        constraints+=len(coeff)
        if coeff:
            small=min(coeff.values());minimum=small if minimum is None or small<minimum else minimum
        contribution=aq(volume(tet))*sum((aq(p).root(3)**4 for p in pp),arb(0))/(4*arb.pi()**4)
        total+=contribution
        cell_records.append({'cell':i,'case':case,'volume':str(volume(tet)),'nonzero_residual_coefficients':len(coeff),'all_nonnegative':True})
        if (i+1)%400==0:print(json.dumps({'stage':'035_integral','checked':i+1,'of':2440,'nonnegative_coefficients':constraints,'elapsed':time.monotonic()-start}),flush=True)
    assert constraints==4932474
    # The nearest true line intersection in this small corner ball is q itself.
    # 030 gives the whole3D radial cap; positive face occupies exactly one octant.
    pi=arb.pi();D=1/(48*pi**4);B=arb(180224)/pi**3;eps=aq(Q(1,2**27))
    cap=arb(2).root(3)*4*pi*(D.root(3)**4*eps**3/3+aq(Q(3,5))*B.root(3)**4*eps.root(3)**5)/8
    Ibound=6*total+cap
    Iupper=Q(str(Ibound.upper().fmpq()));I0=Q(92072351,500000000000);g0=Q(21,25000);candidate=Q(14753,1000)
    assert Iupper<=I0<Q(19,100000)
    costcube=256*I0**3/g0**4;gap=candidate**3-costcube
    expected=json.loads((root/'results/round011_n4_actual_candidate_cost.json').read_text())
    assert costcube==Q(expected['cost_cube_exact']) and gap==Q(expected['strict_cube_gap']) and gap>0 and candidate<16
    cost=aq(costcube).root(3);assert cost<aq(candidate)
    for name,h in sources.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source changed during check')
    return {'status':'PASS independently checked2440cell exact residuals, continuous I4 cap and fixed candidate cost; representation/g acceptance is companion035 check',
            'source_sha256':sources,'cells':2440,'nonzero_nonnegative_residual_coefficients':constraints,'minimum_nonzero_residual':str(minimum),
            'ordered_covered_volume':str(covered),'omitted_ordered_volume':str(omitted),'independent_geometry_matches_certificate':True,'shell_barycentric_identities':3,'cells_summary':cell_records,
            'cap_interval':ends(cap),'cap_ball':cap.str(55),'I4_upper_interval':ends(Ibound),'I4_upper_ball':Ibound.str(55),'I4_upper_exact':str(Iupper),
            'I_outward_rational':str(I0),'g_fixed_rational':str(g0),'cost_cube_exact':str(costcube),'candidate_U4':str(candidate),'strict_cube_gap':str(gap),'cost_interval':ends(cost),'cost_ball':cost.str(55),
            'precision_bits':256,'LP_calls':0,'author_engine_or_checker_calls':0,'shared_libraries':'flint fmpq_mpoly/Arb256bit,Fraction;locked N/D data independently validated in035representation;new geometry/substitution/checker implementation',
            'seconds':time.monotonic()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');a=ap.parse_args();out=run(a.project);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if a.save:(HERE/'035_integral_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','cells_summary')},indent=2))
