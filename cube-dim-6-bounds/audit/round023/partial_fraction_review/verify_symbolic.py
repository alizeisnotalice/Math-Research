#!/usr/bin/env python3
"""One independent full Hermite-DD5 Laurent Green audit, plus eight formal cases."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import gzip,hashlib,json,math,os,platform,resource,signal,sys,time
import sympy as sp
import flint
from flint import fmpq as Q

HERE=Path(__file__).resolve().parent
TOP=Path('/Users/zhengzhihao/Desktop/T')
B=TOP/'Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round023/partial_fraction_kernel'
C=TOP/'Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round022/reduced_point_preflight'
INPUTS={}
CASES=[(1,2,3,4,5,6),(-3,-2,-1,1,2,3),(0,1,2,3,4,5),
       (0,0,-2,-1,1,3),(1,1,-1,-1,2,3),(-2,-2,-2,1,3,4),
       (1,1,1,1,-1,2),(0,0,0,0,0,0)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    INPUTS[str(p)]=sha(p);return json.loads(p.read_text())
def pin(p,h):
    INPUTS[str(p)]=sha(p);assert INPUTS[str(p)]==h,str(p)
def write(name,data):(HERE/name).write_text(json.dumps(data,indent=2)+'\n')
def jet(b,k):
    """[pi/2 constant coefficient, eps^-1, eps^-3, eps^-5]."""
    out=[Q(0) for _ in range(4)]
    if k==0:out[0]=Q(1 if b>0 else -1 if b<0 else 0)
    elif not b and k%2:out[(k+1)//2]=Q((-1)**((k-1)//2),k)
    return out
def hermite_laurent(nodes):
    nodes=sorted(nodes);row=[jet(b,0) for b in nodes]
    for k in range(1,len(nodes)):
        nxt=[]
        for i in range(len(nodes)-k):
            delta=nodes[i+k]-nodes[i]
            nxt.append([(row[i+1][j]-row[i][j])/delta for j in range(4)]
                       if delta else jet(nodes[i],k))
        row=nxt
    return row[0]

def scalar_checks():
    # Independent interpolation matrix, not the full-run divided-difference recurrence.
    x,t=sp.symbols('x t',real=True);g=x/(t*t+x*x)
    derivatives=[sp.diff(g,x,k)/sp.factorial(k) for k in range(6)]
    output=[]
    for case in CASES:
        counts=Counter(case);rows=[(b,k) for b in sorted(counts) for k in range(counts[b])]
        matrix=sp.Matrix([[sp.binomial(d,k)*sp.Integer(b)**(d-k) if d>=k else 0
                           for d in range(6)] for b,k in rows])
        leading=matrix.inv()[5,:]
        interpolation=sp.zeros(6,4)
        for i,(b,k) in enumerate(rows):
            if k==0:interpolation[i,0]=sp.sign(b)
            elif b==0 and k%2:interpolation[i,(k+1)//2]=sp.Rational((-1)**((k-1)//2),k)
        expected=leading*interpolation
        actual=hermite_laurent(list(map(Q,case)))
        assert list(map(str,actual))==[str(expected[j]) for j in range(4)]
        dd_density=sum(leading[i]*derivatives[k].subs(x,b) for i,(b,k) in enumerate(rows))
        numerator=sp.re(sp.expand(sp.prod(t+sp.I*b for b in case)))
        target=numerator/sp.prod(t*t+b*b for b in case)
        assert sp.cancel(dd_density-target)==0
        output.append({'nodes':list(case),'Laurent_jet_interpolation_matches':True,
                       'formal_density_identity':True,'Laurent_vector':list(map(str,actual))})
    assert len(output)==8
    return output

def saved_table_checks(result,tables):
    # This reads saved tables only and is NOT a second full tuple traversal.
    signed={(F(r['b']),r['order']):F(r['coefficient']) for r in tables['signed_A']}
    assert len(signed)==len(tables['signed_A']) and all(signed.values())
    fold={}
    for b in sorted({abs(b) for b,j in signed if b}):
        for j in range(1,7):
            c=signed.get((b,j),F(0))+(-1)**j*signed.get((-b,j),F(0))
            if c:fold[b,j]=c
    assert fold=={(F(r['positive_b']),r['order']):F(r['coefficient']) for r in tables['positive_fold_C']}
    zero={r['order']:F(r['A']) for r in tables['zero_pole_channels']}
    assert all(zero[j]==signed.get((F(0),j),F(0)) for j in range(1,7))
    assert all(zero[j]==0 for j in (2,4,6))
    real=defaultdict(F);powers={}
    for (b,j),c in fold.items():
        q=b*b;powers[q]=max(j,powers.get(q,0));atpole=F(0)
        for r in range(j//2+1):
            n=(-1 if (j-r-3)%2 else 1)*math.comb(j,2*r)*b**(j-2*r)
            atpole+=n*(-q)**r
            for s in range(r+1):real[q,j-s]+=c*n*math.comb(r,s)*(-q)**(r-s)
        assert atpole==(-1)**(j+1)*2**(j-1)*b**j and atpole
    real={k:v for k,v in real.items() if v}
    assert real=={(F(r['b_squared']),r['order']):F(r['coefficient']) for r in tables['real_z_principal']}
    assert all(real[q,j]!=0 for q,j in powers.items())
    assert powers=={F(r['b_squared']):r['exponent'] for r in tables['denominator_factors']}
    assert len(powers)==359 and sum(powers.values())==417
    assert Counter(powers.values())=={1:311,2:38,3:10}
    gamma=sum((c for (b,j),c in fold.items() if j==1),F(0))
    integrated=F(0)
    for (q,j),c in real.items():
        n,d=math.isqrt(q.numerator),math.isqrt(q.denominator)
        assert n*n==q.numerator and d*d==q.denominator
        b=F(n,d);integrated+=c*math.comb(2*j-2,j-1)/(4**(j-1)*b**(2*j-1))
    assert integrated==gamma==F(result['Gamma_G'])
    assert F(result['G_times_pi5_exact'])==gamma/128
    return {'signed_slots':len(signed),'folded_slots':len(fold),'real_principal_slots':len(real),
            'positive_poles':359,'true_degree_z':417,'exponent_histogram':dict(Counter(powers.values())),
            'even_zero_channels_exact_zero':True,'highest_pole_witnesses_checked':359,
            'Gamma_fold_and_real_integral_equal':True}

def main():
    assert not (HERE/'run_receipt.json').exists(),'one complete symbolic pass only'
    started=time.perf_counter()
    pin(B/'manifest.json','8b5eacb4e8bb09ca95ea1abe01edabe140c79e07048a3a76c5e15efe9b77285c')
    manifest=read(B/'manifest.json')
    for row in manifest['files']:pin(B/row['path'],row['sha256'])
    meta=read(B/'input_manifest.json')
    for row in meta['inputs']:pin(Path(row['path']),row['sha256'])
    bresult=read(B/'results.json');tables=read(B/'partial_tables.json')
    compressed=(C/'merged_tuples.json.gz').read_bytes()
    content=gzip.decompress(compressed)
    assert hashlib.sha256(content).hexdigest()=='7deb058449e11fa74d7e1a4658311cb738a38fbbb57ab03118b494138cd39390'
    data=json.loads(content)
    assert data['point_numerators']==[50,50,50,45,27,0] and data['denominator']==50
    assert len(data['tuples'])==40320 and data['logical_terms']==967680
    write('input_manifest.json',{'status':'REGISTERED_BEFORE_SCALARS_AND_ONE_FULL_SYMBOLIC_AUDIT',
          'inputs':INPUTS,'script_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.md'),
          'full_symbolic_cap':1,'scalar_cases':8,'native_cap':0,'positive_time_cap':0,
          'method':'Hermite DD5 atan Laurent recurrence','decompressed_sha256':hashlib.sha256(content).hexdigest()})
    signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(TimeoutError('900-second audit cap')))
    signal.alarm(900)
    formal=scalar_checks()
    write('scalar_checks.json',{'status':'PASS_EIGHT_INDEPENDENT_FORMAL_HERMITE_CASES','cases':formal,
                              'full_symbolic_passes':0,'numeric_times':0})
    full_start=time.perf_counter();total=[Q(0) for _ in range(4)];weight_sum=Q(0)
    record={'status':'ONE_INDEPENDENT_SYMBOLIC_G_AUDIT_RUNNING','pid':os.getpid(),
            'full_symbolic_passes':1,'method':'Hermite-DD5 Laurent of atan(b/epsilon)',
            'tuples_processed':0,'logical_terms':0,'point_orbit_cone_terms':0,
            'zero_tuple_logical_terms':0,'repeated_tuple_logical_terms':0,
            'native_calls':0,'positive_time_values':0,'Sturm_calls':0,'MC_calls':0,'MATLAB_calls':0}
    write('run_receipt.json',record)
    try:
        for row in data['tuples']:
            nodes=tuple(Q(x) for x in row['nodes']);w=Q(row['weight'])
            assert len(nodes)==6 and list(nodes)==sorted(nodes) and w>0
            assert row['global_negation_orbit_size']==2
            assert row['logical_terms']==12*row['point_orbit_cone_terms']
            coefficients=hermite_laurent(nodes)
            for j in range(4):total[j]+=w*coefficients[j]
            weight_sum+=w;record['tuples_processed']+=1
            record['logical_terms']+=row['logical_terms']
            record['point_orbit_cone_terms']+=row['point_orbit_cone_terms']
            counts=Counter(nodes)
            record['zero_tuple_logical_terms']+=row['logical_terms']*int(Q(0) in counts)
            record['repeated_tuple_logical_terms']+=row['logical_terms']*int(max(counts.values())>1)
            if record['tuples_processed']%2048==0:
                rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
                assert rss<768*1024**2
                record['sampled_peak_rss_bytes']=rss
        assert record['tuples_processed']==40320 and record['logical_terms']==967680
        assert record['point_orbit_cone_terms']==80640
        assert record['zero_tuple_logical_terms']==88320 and record['repeated_tuple_logical_terms']==201408
        assert weight_sum==Q(520066141,48769479475200)
        assert total[0]==Q(bresult['Gamma_G']) and all(c==0 for c in total[1:])
        zero={r['order']:Q(r['A']) for r in tables['zero_pole_channels']}
        assert total[1]==zero[2] and total[2]==-zero[4]/3 and total[3]==zero[6]/5
        full_seconds=time.perf_counter()-full_start
        saved=saved_table_checks(bresult,tables)
        result={'status':'PASS_INDEPENDENT_FULL_SYMBOLIC_G_AND_PARTIAL_SOURCE_AUDIT',
                'Gamma_G':str(total[0]),'G_times_pi5_exact':str(total[0]/128),
                'Laurent_negative_coefficients':{'epsilon^-1':str(total[1]),'epsilon^-3':str(total[2]),'epsilon^-5':str(total[3])},
                'Gamma_matches_B023_exactly':True,'all_negative_powers_cancel':True,
                'full_weight_checksum':str(weight_sum),'saved_table_checks':saved,
                'scalar_cases':8,'full_symbolic_passes':1,'additional_density_tuple_passes':0,
                'native_calls':0,'positive_time_values':0,'Sturm_calls':0,'MC_calls':0,'MATLAB_calls':0,
                'no_global_numerator':True,'new_C6_bound':False,'all_time_certificate':False,
                'full_symbolic_seconds':full_seconds,'elapsed_seconds':time.perf_counter()-started,
                'pid':os.getpid(),'python':platform.python_version(),'python_flint':flint.__version__,'sympy':sp.__version__,
                'script_sha256':sha(Path(__file__)),'input_sha256':INPUTS}
        write('results.json',result);record['status']='INDEPENDENT_SYMBOLIC_AUDIT_COMPLETE_STOPPED'
        record['full_symbolic_seconds']=full_seconds
    except BaseException as error:
        record['status']='FAILED_STOPPED';record['exception']=repr(error);raise
    finally:
        record['elapsed_seconds']=time.perf_counter()-started;write('run_receipt.json',record);signal.alarm(0)
    print(json.dumps({k:result[k] for k in ('status','pid','full_symbolic_seconds','elapsed_seconds','full_symbolic_passes','scalar_cases','Gamma_matches_B023_exactly')}))
if __name__=='__main__':main()
