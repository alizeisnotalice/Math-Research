"""Selected019 point only: time-resolvent partial fractions, no DD/Arb."""
import argparse,gzip,hashlib,json,math,os,sys,time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent
RUN='audit_zero041_resolvent_001';N=96;GRID=2**384
MANIFEST='54d220d4627afac14b2f4c705904c8564356ba8c99a834e1cacd22add3dee59d'
TABLE='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
CONTINUITY='298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81'
CHECKPOINT=Path('/tmp/audit_zero041_checkpoint.json')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def project(p):
    p=p.resolve();return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def phase(q,e):
    e%=4
    return ((q,F(0)),(F(0),q),(-q,F(0)),(F(0),-q))[e]
def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def add(a,b):return a[0]+b[0],a[1]+b[1]
def reciprocal(a):
    d=a[0]*a[0]+a[1]*a[1];return a[0]/d,-a[1]/d
def gp(a,j):
    out=(F(1),F(0))
    for _ in range(j):out=mul(out,a)
    return out
def parse(s):
    a,b=map(F,s[1:-1].split(' +/- '));return a-b,a+b
def log_bound(u):
    assert 0<=u<=F(1,3)
    power=u;total=F(0);u2=u*u
    for j in range(N):total+=2*power/(2*j+1);power*=u2
    high=total+2*power/((2*N+1)*(1-u2))
    return (total*GRID).__floor__(),(high*GRID).__ceil__()
def pi_bounds():
    def atan(q):
        low=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(80)),F(0))
        return low,low+F(1,161*q**161)
    a,b=atan(5),atan(239);return 16*a[0]-4*b[1],16*a[1]-4*b[0]
def local_resolvent(nodes):
    """Pole Taylor coefficients: B_(k,j)=i^(j-5)*a_(k,j)."""
    counts=Counter(nodes);out={}
    for k,m in counts.items():
        poly=[F(1)]+[F(0)]*(m-1)
        for l,d in counts.items():
            if l==k:continue
            q=[F((-1)**h*math.comb(d+h-1,h),(k-l)**(d+h)) for h in range(m)]
            poly=[sum((poly[h]*q[t-h] for h in range(t+1)),F(0)) for t in range(m)]
        for j in range(1,m+1):out[k,j]=poly[m-j]
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-root',type=Path,required=True);ap.add_argument('--model-root',type=Path,required=True)
    ap.add_argument('--run',action='store_true');ap.add_argument('--arm',action='store_true');args=ap.parse_args()
    source,model=project(args.source_root),project(args.model_root);author=source/'astra_green019'
    manifest=author/'manifest.json';assert sha(manifest)==MANIFEST
    locked={'astra_green019/manifest.json':MANIFEST}
    for row in json.loads(manifest.read_text())['white_list']:
        p=author/row['path'];assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes'];locked['astra_green019/'+row['path']]=sha(p)
    tablepath=source/'astra_vertex017/boundary_triangulation.json.gz';assert sha(tablepath)==TABLE
    locked['astra_vertex017/boundary_triangulation.json.gz']=TABLE
    contpath=model/'upper_alpha1_029/report.md';assert sha(contpath)==CONTINUITY
    locked['upper_alpha1_029/report.md']=CONTINUITY
    selected=json.loads((author/'point_0.json').read_text())
    x=tuple(map(F,selected['point']));assert x==(F(1),F(1),F(1),F(2,3),F(1,16))
    if not args.run:
        print(json.dumps({'run_handle':RUN,'status':'locked_inputs_only','calculation_executed':False}));return
    pid=os.getpid();print(json.dumps({'run_handle':RUN,'pid':pid,'status':'armed' if args.arm else 'started'}),flush=True)
    if args.arm:input()
    began=time.monotonic()
    def checkpoint(stage,partial):
        if time.monotonic()-began>120:
            save(CHECKPOINT,{'run_handle':RUN,'pid':pid,'status':'STOP_120_second_checkpoint','stage':stage,
                            'partial':partial,'locked_inputs_sha256':locked,'seconds':time.monotonic()-began})
            raise TimeoutError(str(CHECKPOINT))
    table=json.loads(gzip.decompress(tablepath.read_bytes()));assert table['denominator']==144
    Z=table['vertices'];simplices=table['simplices'];assert len(Z)==280 and len(simplices)==3504
    L=math.lcm(*(v.denominator for v in x));n=[int(v*L) for v in x]
    k=[sum(a*b for a,b in zip(n,z)) for z in Z]
    assert L==48 and n==[48,48,48,32,3] and k.count(0)==6 and len(set(k))==95
    assert len(set(map(abs,x)))==3 # >2 absolute coordinate values excludes all body pair spans
    signatures=Counter();total_D=0
    for s in simplices:
        signatures[tuple(sorted(k[v] for v in s['vertices']))]+=s['det_numerator'];total_D+=s['det_numerator']
    A={}
    for index,(nodes,D) in enumerate(sorted(signatures.items())):
        for key,value in local_resolvent(nodes).items():A[key]=A.get(key,F(0))+D*value
        checkpoint('partial_fraction_signatures',{'completed':index+1,'total':len(signatures),
                    'coefficients':[[k,j,str(a)] for (k,j),a in sorted(A.items())]})
    zero=[phase(A.get((0,j),F(0)),j-5) for j in range(1,6)]
    assert all(re==0 for re,im in zero)
    laurent=[]
    for p in range(1,6):
        q=sum((a*math.comb(p-1,j-1)*k**(p-j) for (k,j),a in A.items() if j<=p),F(0))
        assert q==(total_D if p==5 else 0);laurent.append(str(q))
    # Different Gaussian-rational code path, evaluate rational identity at s=1.
    direct=(F(0),F(0))
    for nodes,D in signatures.items():
        term=(F(D),F(0))
        for kv in nodes:term=mul(term,reciprocal((F(1),F(-kv))))
        direct=add(direct,term)
    expansion=(F(0),F(0))
    for (kv,j),a in A.items():
        expansion=add(expansion,mul(phase(a,j-5),gp(reciprocal((F(1),F(-kv))),j)))
    assert direct==expansion
    logs={};rational=F(0)
    for (kv,j),a in A.items():
        if not kv:continue # only after ALL real zero-pole coefficients verified zero
        if j==1:logs[abs(kv)]=logs.get(abs(kv),F(0))-a
        else:rational+=F((-1)**(j+1),(j-1)*kv**(j-1))*a
    logs={k:c for k,c in logs.items() if c};assert sum(logs.values(),F(0))==0
    prefactor=F(L**4,144*32);logs={k:prefactor*c for k,c in logs.items()};rational*=prefactor
    ref=max(logs);terms={};eweight=F(0)
    for integer,c in logs.items():
        if integer==ref:continue
        y=F(integer,ref);e=0
        while y<1:y*=2;e-=1
        while y>=2:y/=2;e+=1
        u=(y-1)/(y+1)
        if u:terms[u]=terms.get(u,F(0))+c
        eweight+=e*c
    terms[F(1,3)]=terms.get(F(1,3),F(0))+eweight
    terms={u:c for u,c in terms.items() if c};assert len(terms)<=140
    coeff={'point':[str(v) for v in x],'L':L,'n':n,'integer_projections':k,
           'phase_convention':'B_(k,j)=i^(j-5) a_(k,j), resolvent=sum B/(s-ik)^j',
           'all_resolvent_a_coefficients':[[kv,j,str(a)] for (kv,j),a in sorted(A.items())],
           'zero_pole_Gaussian_rational_coefficients':[[j,str(z[0]),str(z[1])] for j,z in enumerate(zero,1)],
           'prefactor_for_G_times_pi5':str(prefactor),'G_pi5_rational_remainder':str(rational),
           'G_pi5_log_coefficients':{str(k):str(c) for k,c in sorted(logs.items())},
           'sum_log_coefficients':'0','log_reference':ref,
           'rescaled_atanh_terms':{str(u):str(c) for u,c in sorted(terms.items())}}
    cp=HERE/'coefficients.json';save(cp,coeff)
    common=math.lcm(rational.denominator,*(c.denominator for c in terms.values()))
    low_num=high_num=rational.numerator*(common//rational.denominator)*GRID
    for index,(u,c) in enumerate(sorted(terms.items())):
        low,high=log_bound(u);weight=c.numerator*(common//c.denominator)
        low_num+=weight*(low if weight>=0 else high);high_num+=weight*(high if weight>=0 else low)
        checkpoint('log_series',{'completed':index+1,'total':len(terms),'coefficients_sha256':sha(cp)})
    low,high=F(low_num,common*GRID),F(high_num,common*GRID);assert 0<low<high
    pl,ph=pi_bounds();gl,gh=low/ph**5,high/pl**5;scale=10**85
    gl,gh=F((gl*scale).__floor__(),scale),F((gh*scale).__ceil__(),scale)
    printed=parse(selected['G_ball']);dyadic=F(selected['G_lower_dyadic']),F(selected['G_upper_dyadic'])
    assert 0<printed[0]<gl<gh<printed[1] and dyadic[0]<gl<gh<dyadic[1]
    out={'run_handle':RUN,'pid':pid,'status':'PASS_selected019_first_point_time_resolvent_fullGreen',
         'point':[str(v) for v in x],'same_alpha1_model_and_Fourier_factor':True,
         'partial_fractions':{'simplices':3504,'distinct_integer_projections':len(set(k)),
             'zero_projection_vertices':k.count(0),'signatures':len(signatures),
             'maximum_pole_order':max(j for k,j in A),'all_real_zero_principal_coefficients_zero':True,
             'zero_Gaussian_rational_coefficients':[[str(a),str(b)] for a,b in zero],
             'Laurent_s_minus_1_through_5_rational_coefficients':laurent,
             'direct_3504_resolvents_equal_partial_fractions_at_s1':True},
         'integration':{'G_pi5_rational_remainder':str(rational),'log_count':len(logs),
             'log_coefficients_sum_zero':True,'atanh_series_count':len(terms),'N':N,'log_quantization_bits':384,
             'coefficient_common_denominator_bits':common.bit_length(),'Machin_terms_per_arctan':80},
         'pure_Fraction_G_interval':[str(gl),str(gh)],'inside_source_first_G_ball_and_dyadic':True,
         'essential_obstacle':{'outside_S2_by_three_distinct_absolute_coordinates':True,
             'continuity_uses_locked_029_five_dimensional_Cauchy_proof':True,
             'not_using_nonzero_projection_DCT':True,'not_a_global_g_lower_or_minimizer_or_C5_bound':True},
         'coefficients_sha256':sha(cp),'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'registration_sha256':sha(HERE/'registration.md'),
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began,'point_evaluations':1,
                    'transcendental_library_calls':0,'author_modules_imported':0},
         'scope':'First019 point only; geometry and five-dimensional continuity remain prior locked inputs. No DD/Hermite functional/Arb, other019 points,017/036 anchors,oldfan,search,gate or Gamma baseline executed.'}
    save(HERE/'result.json',out)
    print(json.dumps({k:out[k] for k in ('run_handle','pid','status','partial_fractions','integration','pure_Fraction_G_interval','runtime')},indent=2))

if __name__=='__main__':main()
