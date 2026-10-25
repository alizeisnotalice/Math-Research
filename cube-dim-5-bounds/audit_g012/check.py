"""Portable independent audit. Default locks inputs; --run executes fixed batch."""
import argparse,hashlib,json,math,subprocess,sys,tempfile,time
from fractions import Fraction as F
from itertools import combinations,permutations,product
from pathlib import Path

HERE=Path(__file__).resolve().parent
MANIFEST='92398b5b1d5b5e26552fbaadbab4e3627d34a02fa26d9f6ce1f65e5a5d832da4'
EDGES=list(combinations(range(5),2))
V=[(1,)+v for v in product((-1,1),repeat=4)]
TESTS=[(2,0,0,0,0,0,0,0,0,0),(4,0,0,0,0,0,0,0,0,0),(12,0,0,0,0,0,0,0,0,0),
       (1,1,0,0,1,0,0,0,0,0),(5,5,0,0,5,0,0,0,0,0),
       (3,0,3,0,3,0,0,3,0,0),(3,0,0,3,3,0,0,3,0,3),
       (1,1,1,1,1,1,1,1,1,1),(8,8,0,0,0,0,0,0,0,0)]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def load_poly(poly):return {tuple(map(int,k.split(','))):F(v) for k,v in poly.items()}

def direct_stieltjes(shape):
    # Direct multivariate power-series reciprocal, without canonical graphs.
    beta=F(3,4);zero=(0,)*10;n=sum(shape)
    indices=sorted(product(*(range(k+1) for k in shape)),key=sum)
    choose=[F(1)]
    for j in range(1,n+1):choose.append(choose[-1]*(beta-j+1)/j)
    den={};inv={zero:F(1)}
    for k in indices[1:]:
        avg=F(sum(math.prod((v[i]*v[j])**power for (i,j),power in zip(EDGES,k)) for v in V),16)
        if avg:
            degree=sum(k)
            multi=math.factorial(degree)//math.prod(math.factorial(p) for p in k)
            den[k]=choose[degree]*(-1)**degree*multi*avg
        inv[k]=-sum((value*inv[tuple(a-b for a,b in zip(k,j))]
                     for j,value in den.items() if all(b<=a for a,b in zip(k,j))),F(0))
    return inv[shape]*math.prod(math.factorial(k) for k in shape)/math.prod((beta+j for j in range(n)),start=F(1))

def atan_bounds(q):
    lower=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(20)),F(0))
    return lower,lower+F(1,41*q**41)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,required=True)
    parser.add_argument('--run',action='store_true')
    args=parser.parse_args();project=args.project_root.resolve()
    if (project/'cube-dim-5-bounds').is_dir():project=project/'cube-dim-5-bounds'
    author=project/'astra_upper004';manifest=author/'manifest.json'
    assert sha(manifest)==MANIFEST
    locked={str(manifest.relative_to(project)):sha(manifest)}
    for name,info in json.loads(manifest.read_text()).items():
        p=author/name;assert sha(p)==info['sha256'] and p.stat().st_size==info['bytes']
        locked[str(p.relative_to(project))]=sha(p)
    checkpoint=json.loads((author/'corner_checkpoint.json').read_text())
    for legacy,expected in checkpoint['source_sha256'].items():
        suffix=Path(legacy).parts[-2:];p=project.joinpath(*suffix)
        assert sha(p)==expected;locked[str(p.relative_to(project))]=sha(p)
    if not args.run:
        print(json.dumps({'status':'locked_inputs_verified_only','run':False,'locked_count':len(locked)}));return
    began=time.monotonic()
    verification=json.loads((author/'verification.json').read_text())
    data=json.loads((author/'space_moments.json').read_text())
    families={name:[load_poly(poly) for poly in fam] for name,fam in data['families'].items()}
    A,B,c=(F(checkpoint['constants'][key]) for key in ('A_lower','B_upper','c_lower'))
    tau=F(15800608,26277615)
    # Endpoints verified without Arb using positive eleventh and second powers.
    assert A>B>0 and A**11<F(17,11)**11/tau**6
    assert B**11>F(6,11)**11/tau**17
    a5,a239=atan_bounds(5),atan_bounds(239)
    pi_upper=16*a5[1]-4*a239[0]
    assert c>0 and c*c*2**15*pi_upper**5<9
    # Check exact W_v vertex extrema, proving W>=-2 by multiaffinity and convexity.
    wvertices=[sum(v[i]*v[j]*x[i]*x[j] for i,j in EDGES) for v in V for x in product((-1,1),repeat=5)]
    assert min(wvertices)==-2
    eta,a=F(1,8),F(333,500)
    assert a*(1-2*eta)==F(999,2000)>0
    power={}
    for n in range(12):
        if n<=7:
            for fam,coef in (('d1',A),('d2',-B)):
                factor=F(11,4)*a**7*math.comb(7,n)*eta**n*coef
                for key,val in families[fam][n].items():power[key]=power.get(key,F(0))+factor*val
        factor=-F(7,4)*a**11*math.comb(11,n)*eta**n
        for key,val in families['b'][n].items():power[key]=power.get(key,F(0))+factor*val
    power={k:v for k,v in power.items() if v}
    target=next(row for row in verification['candidates'] if row['eta']=='1/8')
    assert power==load_poly(target['symmetric_power_coefficients'])
    face={}
    for partition,value in power.items():
        assert tuple(sorted(partition,reverse=True))==partition
        for monomial in set(permutations(partition)):
            key=monomial[1:];face[key]=face.get(key,F(0))+value
    face={k:v for k,v in face.items() if v}
    degrees=tuple(max(k[i] for k in face) for i in range(4));assert degrees==(6,6,6,6)
    D=math.lcm(*(value.denominator for value in face.values()))
    indices=list(product(range(7),repeat=4));tensor={k:0 for k in indices}
    # 60 is lcm of binomial(6,j). All inner computations are arbitrary integers.
    weights=[[math.comb(k,j)*(60//math.comb(6,j)) if j<=k else 0 for j in range(7)] for k in range(7)]
    for j,value in face.items():
        coeff=value.numerator*(D//value.denominator)
        for k in product(*(range(exponent,7) for exponent in j)):
            tensor[k]+=coeff*math.prod(weights[ki][ji] for ki,ji in zip(k,j))
    tensor={k:F(v,D*60**4) for k,v in tensor.items()}
    minimum=min(tensor.values());minimum_keys=[k for k,v in tensor.items() if v==minimum]
    assert len(tensor)==2401 and minimum==F(target['minimum_Bernstein_value'])
    assert minimum>0 and c*minimum>F(7,156250)
    tensor_hash=hashlib.sha256(json.dumps({','.join(map(str,k)):str(v) for k,v in sorted(tensor.items())},sort_keys=True).encode()).hexdigest()
    assert tensor_hash==target['Bernstein_coefficients_sha256']
    with tempfile.TemporaryDirectory(prefix='audit-g012-') as temp:
        temp=Path(temp);binary=temp/'spots';rawpath=temp/'raw.json';spotpath=temp/'spots.tsv'
        # Renaming main for the probe removes C++'s special implicit return 0.
        # Add that same explicit return in a temporary copy; no author file edit.
        text=(author/'graph_moments.cpp').read_text()
        assert text.endswith('}\n')
        instrumented=temp/'graph_source.cpp'
        instrumented.write_text(text[:-2]+'return 0;\n}\n')
        compiler=['clang++','-O2','-std=c++17','-DAUTHOR_SOURCE="'+str(instrumented)+'"',str(HERE/'moment_spots.cpp'),'-o',str(binary)]
        subprocess.run(compiler,check=True,capture_output=True,text=True,timeout=60)
        compiled=subprocess.run([str(binary),str(rawpath),str(spotpath)],check=True,capture_output=True,text=True,timeout=60)
        assert sha(rawpath)==verification['raw_graph_moment_input_sha256']
        raw=json.loads(rawpath.read_text());den=list(map(int,raw['denominators']))
        assert den==[math.prod(4*j+3 for j in range(n)) for n in range(18)]
        for fam,polys in families.items():
            decoded=[{tuple(map(int,k.split(','))):sum((F(int(num),den[i]) for i,num in enumerate(row)),F(0)) for k,row in poly.items()} for poly in raw[fam]]
            decoded=[{k:v for k,v in poly.items() if v} for poly in decoded]
            assert decoded==polys
        expected={int(code):F(int(num),int(denom)) for code,num,denom in (line.split('\t') for line in spotpath.read_text().splitlines())}
        spot_receipts=[]
        for shape in TESTS:
            code=sum(k<<(4*i) for i,k in enumerate(shape));value=direct_stieltjes(shape)
            assert value==expected[code]
            spot_receipts.append({'edge_multiplicities':shape,'degree':sum(shape),'independent_exact_moment':str(value)})
        rawhash=sha(rawpath)
    out={'status':'pass_independent_Young_Bernstein_and_nine_moment_checks',
         'audited_candidate':{'eta':str(eta),'a':str(a),'g':'7/156250'},
         'exact_endpoints':{'A_lower':str(A),'B_upper':str(B),'c_lower':str(c),'rational_power_checks':True},
         'Z_uniform_lower':str(a*(1-2*eta)),
         'Bernstein':{'degrees':degrees,'coefficients':len(tensor),'minimum_keys':minimum_keys,
                      'minimum':str(minimum),'c_lower_times_minimum':str(c*minimum),'coefficient_sha256':tensor_hash},
         'independent_spot_moments':spot_receipts,'regenerated_raw_graph_sha256':rawhash,
         'canonical_moments':raw['canonical_moment_count'],'canonical_cache':raw['canonical_cache_count'],
         'dependency_limit':'All spatial moments use the locked author C++ generator; nine moments independently use Python Fraction multivariate series. No independent rederivation of all 3389 numeric moments is claimed.',
         'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'driver_sha256':sha(HERE/'moment_spots.cpp'),'registration_sha256':sha(HERE/'registration.md'),
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began}}
    (HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','Bernstein','canonical_moments','runtime')},indent=2))

if __name__=='__main__':main()
