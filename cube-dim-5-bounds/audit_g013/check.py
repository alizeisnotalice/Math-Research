"""Portable finite audit; default hashes only, --run executes the registered batch."""
import argparse,hashlib,json,math,subprocess,sys,tempfile,time
from collections import Counter
from fractions import Fraction as F
from itertools import permutations,product
from pathlib import Path

HERE=Path(__file__).resolve().parent
MANIFEST='05bda95ed6fdd4a10a676c0b3a9dbd558130050a96969b9390ee418668b5523c'
HEADER='7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def project_path(p):
    p=p.resolve()
    return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def poly(p):return {tuple(map(int,k.split(','))):F(v) for k,v in p.items()}
def pmul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
def ppow(a,n):
    out=[F(1)]
    for _ in range(n):out=pmul(out,a)
    return out
def corner(p):
    return sum((v*F(math.factorial(5),math.prod(math.factorial(n) for n in Counter(k).values())) for k,v in p.items()),F(0))
def pi_upper():
    def atan(q):
        lo=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(20)),F(0))
        return lo,lo+F(1,41*q**41)
    a,b=atan(5),atan(239)
    return 16*a[1]-4*b[0]

def independent_heldout():
    # Four active characters, ordinary reciprocal Taylor coefficients.
    edges=((0,1),(0,3),(1,2),(2,3));shape=(5,)*4;beta=F(3,4)
    V=[(1,)+v for v in product((-1,1),repeat=4)]
    choose=[F(1)]
    for n in range(1,21):choose.append(choose[-1]*(beta-n+1)/n)
    zero=(0,)*4;inv={zero:F(1)};den={}
    for k in sorted(product(range(6),repeat=4),key=sum)[1:]:
        avg=F(sum(math.prod((v[i]*v[j])**p for (i,j),p in zip(edges,k)) for v in V),16)
        if not avg:inv[k]=F(0);continue
        n=sum(k);multi=math.factorial(n)//math.prod(math.factorial(p) for p in k)
        den[k]=choose[n]*(-1)**n*multi*avg
        inv[k]=-sum((value*inv[tuple(x-y for x,y in zip(k,j))] for j,value in den.items()
                    if all(y<=x for x,y in zip(k,j))),F(0))
    return inv[shape]*math.factorial(5)**4/math.prod((beta+j for j in range(20)),start=F(1))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,required=True)
    parser.add_argument('--backend-root',type=Path,required=True)
    parser.add_argument('--gmp-prefix',type=Path,default=Path('/opt/homebrew/opt/gmp'))
    parser.add_argument('--run',action='store_true');args=parser.parse_args()
    project=project_path(args.project_root);backend=project_path(args.backend_root)/'upper'
    author=project/'astra_upper005';manifest=author/'manifest.json'
    assert sha(manifest)==MANIFEST
    locked={'astra_upper005/manifest.json':sha(manifest)}
    for name,info in json.loads(manifest.read_text())['files'].items():
        p=author/name;assert sha(p)==info['sha256'] and p.stat().st_size==info['bytes']
        locked['astra_upper005/'+name]=sha(p)
    assert sha(backend/'round012_graph_moments_gmp.hpp')==HEADER
    locked['backend/round012_graph_moments_gmp.hpp']=HEADER
    checkpoint=json.loads((author/'checkpoint.json').read_text())
    verification=json.loads((author/'verification.json').read_text())
    for receipt in (checkpoint,verification):
        for legacy,expected in receipt['source_sha256'].items():
            p=project.joinpath(*Path(legacy).parts[-2:]);assert sha(p)==expected
            locked[str(p.relative_to(project))]=sha(p)
    if not args.run:
        print(json.dumps({'status':'locked_inputs_only','calculation_executed':False,'locked_count':len(locked)}));return
    began=time.monotonic()
    families={key:[poly(p) for p in values] for key,values in json.loads((author/'space_moments.json').read_text())['families'].items()}
    keys={'d1':'E_d_W_power','d2':'E_d_squared_W_power','b':'E_b_corner_W_power'}
    for fam,values in families.items():
        for n,p in enumerate(values):
            assert all(sum(k)==n+(fam=='b') for k in p)
            assert corner(p)==F(checkpoint['moments_to_order_22'][keys[fam]][n])
    A,B,c=(F(checkpoint['constants'][k]) for k in ('A_lower','B_upper','c_lower'))
    tau=F(15800608,26277615)
    assert A>B>0 and A**11<F(17,11)**11/tau**6 and B**11>F(6,11)**11/tau**17
    assert c*c*2**15*pi_upper()**5<9
    kappa,a=F(1,256),F(27,40);f=[F(1),F(1,8),-kappa]
    q=lambda w:sum(v*w**i for i,v in enumerate(f))
    assert q(F(-2))==F(47,64)>0 and q(F(10))==F(119,64)>0
    f7,f11=ppow(f,7),ppow(f,11);power={}
    for n,coef in enumerate(f7):
        for fam,scale in (('d1',A),('d2',-B)):
            factor=F(11,4)*a**7*coef*scale
            for key,val in families[fam][n].items():power[key]=power.get(key,F(0))+factor*val
    for n,coef in enumerate(f11):
        factor=-F(7,4)*a**11*coef
        for key,val in families['b'][n].items():power[key]=power.get(key,F(0))+factor*val
    power={k:v for k,v in power.items() if v}
    target=next(row for row in verification['candidates'] if row['kappa']=='1/256')
    assert power==poly(target['symmetric_power_coefficients']) and target['a']=='27/40'
    face={}
    for partition,value in power.items():
        assert partition==tuple(sorted(partition,reverse=True))
        for exponent in set(permutations(partition)):
            key=exponent[1:];face[key]=face.get(key,F(0))+value
    face={k:v for k,v in face.items() if v}
    degrees=tuple(max(k[i] for k in face) for i in range(4));assert degrees==(12,)*4
    n=13;N=n**4;D=math.lcm(*(v.denominator for v in face.values()))
    flatten=lambda k:sum(ki*n**(3-i) for i,ki in enumerate(k))
    tensor=[0]*N
    for k,v in face.items():tensor[flatten(k)]=v.numerator*(D//v.denominator)
    L=math.lcm(*(math.comb(12,j) for j in range(13)))
    weights=[[math.comb(k,j)*(L//math.comb(12,j)) for j in range(k+1)] for k in range(13)]
    # Four exact integer binomial transforms; common denominator stays outside.
    for axis in range(4):
        stride=n**(3-axis);out=[0]*N
        for outer in range(0,N,n*stride):
            for inner in range(stride):
                base=outer+inner
                for k in range(13):out[base+k*stride]=sum(tensor[base+j*stride]*w for j,w in enumerate(weights[k]))
        tensor=out
    denominator=D*L**4;minimum_integer=min(tensor);minimum=F(minimum_integer,denominator)
    minimum_indices=[k for k in product(range(13),repeat=4) if tensor[flatten(k)]==minimum_integer]
    assert minimum==F(target['minimum_Bernstein_value']) and minimum>0
    assert c*minimum>F(4583,10**8)
    coeff={','.join(map(str,k)):str(F(tensor[flatten(k)],denominator)) for k in product(range(13),repeat=4)}
    coefficient_hash=hashlib.sha256(json.dumps(coeff,sort_keys=True).encode()).hexdigest()
    assert coefficient_hash==target['Bernstein_coefficients_sha256']
    with tempfile.TemporaryDirectory(prefix='audit-g013-') as temporary:
        temporary=Path(temporary);binary=temporary/'heldout'
        command=['clang++','-O2','-std=c++17','-I',str(backend),'-isystem',str(args.gmp_prefix/'include'),
                 str(HERE/'heldout.cpp'),'-L',str(args.gmp_prefix/'lib'),'-lgmpxx','-lgmp','-o',str(binary)]
        subprocess.run(command,check=True,capture_output=True,text=True,timeout=30)
        result=subprocess.run([str(binary)],check=True,capture_output=True,text=True,timeout=15)
        heldout=F(result.stdout.strip());independent=independent_heldout();assert heldout==independent
    out={'status':'pass_quadratic_Z_global_g_audit_with_shared_spatial_moments',
         'candidate':{'kappa':'1/256','a':'27/40','g':'4583/100000000'},
         'q_endpoint_values':['47/64','119/64'],'Z_uniform_positive_lower':str(a*F(47,64)),
         'Bernstein':{'degrees':degrees,'coefficient_count':N,'minimum_indices':minimum_indices,
                      'minimum':str(minimum),'c_lower_times_minimum':str(c*minimum),'coefficients_sha256':coefficient_hash},
         'heldout':{'graph':'four-cycle all four edges multiplicity 5','total_degree':20,
                    'exact_moment':str(independent),'fresh_GMP_no_journal_matches_Python_Fraction_Stieltjes':True},
         'dependency_limit':'All saved spatial moments remain locked shared input from author production; no 217-second production replay or complete independent 143289-moment reconstruction.',
         'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'heldout_source_sha256':sha(HERE/'heldout.cpp'),'registration_sha256':sha(HERE/'registration.md'),
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began}}
    (HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ('status','Bernstein','heldout','runtime')},indent=2))

if __name__=='__main__':main()
