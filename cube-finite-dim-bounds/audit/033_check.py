"""Fresh full136-patch minimum audit at384bit, no disk node cache.

New exact-rational real-Hermite point evaluator and polynomial construction.
Only locked cb11 prepared ray/rate data and the Arb arithmetic engine are shared.
No author checker/point evaluator or cached node values are executed/read.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import argparse,hashlib,json,runpy,time,platform
import sympy as sy
from flint import arb,arb_mat,ctx

HERE=Path(__file__).resolve().parent
N=21;DEGREE=20;BITS=128;PRECISION=384
EPS=F(1,50);HALFWIDTH=F(1,32);RADIUS=F(1,50);RHO=F(3)
STRONG=F(131,25000)
def aq(x):
    x=F(x);return arb(x.numerator)/x.denominator
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]

def run(root):
    start=time.monotonic();snapshot=json.loads((HERE/'033_source_snapshot.json').read_text())
    for name,h in snapshot.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    author=json.loads((root/'upper/continuous_certificate/cc11_cheb_all_result.json').read_text())
    # Read-only prepared rates from the independently audited032 source.
    module=runpy.run_path(str(root/'upper/continuous_bridge/cb11_frame_green.py'))
    prepared=module['_prepared'](F(1,3));assert len(prepared)==48
    ctx.prec=PRECISION
    rays={}
    for cone,rates,w in prepared:
        assert w==2/(rates[0]*rates[1]*rates[2])
        for ray,k in zip(cone,rates):
            assert k>=3*(abs(ray[1])+abs(ray[2]))
            rays[ray,k]=None
    min_ratio=min(k/F(abs(v[1])+abs(v[2])) for v,k in rays if abs(v[1])+abs(v[2]))
    assert min_ratio==3
    EPS2=EPS*EPS
    point_calls=0
    def point(y,z):
        nonlocal point_calls
        point_calls+=1
        projection={(ray,k):(F(ray[0])+ray[1]*y+ray[2]*z)/k for ray,k in rays}
        logs={};divdiff={}
        def value(b):
            b=abs(b)
            if b not in logs:logs[b]=aq(EPS2+b*b).log()/2
            return logs[b]
        total=arb(0)
        for cone,rates,w in prepared:
            bs=tuple(sorted(projection[ray,k] for ray,k in zip(cone,rates)))
            if bs not in divdiff:
                lo,mid,hi=bs
                if lo==hi:
                    dd=aq((EPS2-lo*lo)/(2*(EPS2+lo*lo)**2))
                elif lo==mid or mid==hi:
                    a,b=(lo,hi) if lo==mid else (hi,lo)
                    dd=(value(b)-value(a)-aq(a/(EPS2+a*a))*aq(b-a))*aq(1/(b-a)**2)
                else:
                    dd=sum((value(b)*aq(1/((b-bs[(j+1)%3])*(b-bs[(j+2)%3]))) for j,b in enumerate(bs)),arb(0))
                divdiff[bs]=dd
            total+=aq(w)*divdiff[bs]
        return total/(8*arb.pi()**3)
    # Global analytic constants, independently simplified as rational factors.
    strip_delta=EPS-HALFWIDTH*(RHO-1/RHO)/6
    assert strip_delta==F(11,1800)>0
    derivative_delta=EPS-RADIUS/3;assert derivative_delta==F(1,75)
    M=aq(F(15,1088)/strip_delta**2)/arb.pi()**3
    D=aq(F(15,1088)/(derivative_delta**2*RADIUS))/arb.pi()**3
    displacement=2*D*aq(HALFWIDTH/F(2**BITS))
    error=aq(F(84,3**20))*M
    assert displacement<aq(F(2296921,10**44)) #2.296921e-38
    assert error<aq(F(2868313237155,10**19)) #2.868313237155e-7
    # Exact symmetric dyadic replacement of true Chebyshev roots.
    roots=[None]*N;root_records=[]
    for j in range(DEGREE//2):
        exact=(arb.pi()*(2*j+1)/(2*N)).cos()
        floor=(exact*2**BITS).floor().unique_fmpz();assert floor is not None
        q=F(int(floor),2**BITS);assert abs(exact-aq(q))<aq(F(1,2**BITS))
        roots[j]=q;roots[DEGREE-j]=-q
        root_records.append({'index':j,'dyadic':str(q),'true_root_interval':ends(exact)})
    roots[DEGREE//2]=F(0)
    assert len(set(roots))==N and all(abs(q)<1 for q in roots)
    A=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos()*aq(F(1 if k==0 else 2,N)) for j in range(N)] for k in range(N)])
    Q=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos() for k in range(N)] for j in range(N)])
    identity=A*Q
    for i in range(N):
        for j in range(N):assert identity[i,j].contains(int(i==j)) and identity[i,j].rad()<arb(2)**-300
    # SymPy polynomial coefficients rather than author's custom recurrence.
    u=sy.symbols('u');polys=[sy.Poly(sy.chebyshevt(k,2*u-1),u) for k in range(N)]
    power=[[F(polys[k].nth(l)) for l in range(k+1)] for k in range(N)]
    bf=[[sum((power[k][l]*F(comb(j,l),comb(DEGREE,l)) for l in range(min(j,k)+1)),F(0)) for k in range(N)] for j in range(N)]
    for k in range(N):
        back=sum((sy.Rational(bf[j][k].numerator,bf[j][k].denominator)*comb(DEGREE,j)*u**j*(1-u)**(DEGREE-j) for j in range(N)),sy.Integer(0))
        assert sy.Poly(sy.expand(back)-polys[k].as_expr(),u).is_zero
    B=arb_mat([[aq(v) for v in row] for row in bf]);T=B*A
    print(json.dumps({'audit':'033','fresh_nodes':136*N*N,'precision_bits':ctx.prec,'cache':'disabled; no cache files opened or written','target':str(STRONG),'algebra_checks':'PASS all21 exact polynomial identities and441 DCT entries'}),flush=True)
    digest=hashlib.sha256();records=[];max_node_radius=arb(0);max_bernstein_radius=arb(0);minimum=None;minimum_patch=None
    patches=[(i,j) for i in range(16) for j in range(i+1)]
    assert len(patches)==len(set(patches))==136
    source_patches={tuple(r['patch']) for r in author['patches']};assert source_patches==set(patches)
    for count,(i,j) in enumerate(patches,1):
        patch_start=time.monotonic();ys=[F(2*i+1,32)+t/32 for t in roots];zs=[F(2*j+1,32)+t/32 for t in roots]
        values=[]
        for y in ys:
            row=[]
            for z in zs:
                v=point(y,z);assert v.is_finite()
                v+=arb(0,displacement.upper())
                max_node_radius=max(max_node_radius,v.rad())
                digest.update((str(y)+'|'+str(z)+'|'+str(v.lower().fmpq())+'|'+str(v.upper().fmpq())+'\n').encode())
                row.append(v)
            values.append(row)
        C=T*arb_mat(values)*T.transpose()
        poly_lower=min(F(str(C[a,b].lower().fmpq())) for a in range(N) for b in range(N))
        bound=aq(poly_lower)-error
        low=F(str(bound.lower().fmpq()))
        assert low>STRONG,(i,j,'root lower failed strong target; stop, no split/precision changes',str(low))
        max_bernstein_radius=max(max_bernstein_radius,max(C[a,b].rad() for a in range(N) for b in range(N)))
        if minimum is None or low<minimum:minimum=low;minimum_patch=[i,j]
        records.append({'patch':[i,j],'closed_box':[str(F(i,16)),str(F(i+1,16)),str(F(j,16)),str(F(j+1,16))],
                        'fresh_nodes':N*N,'H_lower_exact':str(low),'H_lower_display':str(aq(low)),
                        'strong_margin_exact':str(low-STRONG),'degree':DEGREE,'subdivisions':0,'seconds':time.monotonic()-patch_start})
        if count%8==0 or count==136:print(json.dumps({'done':count,'total':136,'fresh_nodes':point_calls,'minimum_lower':str(aq(minimum)),'minimum_patch':minimum_patch,'elapsed':time.monotonic()-start}),flush=True)
    assert point_calls==59976 and all(r['subdivisions']==0 for r in records)
    # Recheck all source hashes after the long run, before publishing.
    for name,h in snapshot.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source changed while fresh run')
    return {'status':'PASS fresh full-domain strong g>=131/25000 certificate, no J or Cn upper',
            'source_sha256':snapshot,'precision_bits':ctx.prec,'configuration':{'degree':20,'nodes_per_axis':21,'node_bits':128,'epsilon':'1/50','rho':'3','patch_grid':16,'strong_target':str(STRONG)},
            'fresh_point_calls':point_calls,'cache_reads':0,'cache_writes':0,'all_node_enclosures_stream_sha256':digest.hexdigest(),
            'tube_ray_min_ratio':str(min_ratio),'ellipse_delta_exact':str(strip_delta),'derivative_delta_exact':str(derivative_delta),
            'ellipse_M_interval':ends(M),'Cauchy_derivative_interval':ends(D),'node_displacement_value_error_interval':ends(displacement),'interpolation_error_interval':ends(error),
            'positive_root_records':root_records,'exact_polynomial_identities':21,'DCT_identity_entries':441,
            'largest_node_radius':str(max_node_radius),'largest_Bernstein_radius':str(max_bernstein_radius),
            'minimum_H_lower_exact':str(minimum),'minimum_H_lower_display':str(aq(minimum)),'minimum_patch':minimum_patch,
            'strong_target_margin_exact':str(minimum-STRONG),'patch_count':136,'patches':records,'closed_coverage_complete':True,
            'implementation':'fresh own real-Hermite384bit evaluator and SymPy polynomial-derived Bernstein conversion; only locked032 module prepared rays/rates shared; no author point/checker/cache calls',
            'shared_libraries':{'SymPy':sy.__version__,'Fraction':'stdlib','Arb':'python-flint384bit'},'python':platform.python_version(),
            'scope':'fixed tau1/3 full cube pointwise G>=131/25000 via H1/50 lower; not J/cost/newU3; no adaptive branch/parameter/grid search','seconds':time.monotonic()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    out=run(args.project);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:(HERE/'033_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('patches','source_sha256','positive_root_records')},indent=2),flush=True)
