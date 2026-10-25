"""One portable fixed K3 audit; hashes only by default, --run calculates."""
import argparse,hashlib,json,math,sys,time
from fractions import Fraction as F
from itertools import permutations,product
from pathlib import Path

HERE=Path(__file__).resolve().parent
M7='b8a7b09d6327428facdc8e6bc7524a87e226d42f2bf859e41240a42392672b5a'
M6='398f33d8565a30a68b9d4afbb5d29e1047a2a3a37783ab91a56127437bece136'
HEADER='7ae71c1a9ddc0d4a1b6a8bd94393988adc55bf4154a8ee3afa692337f97391a8'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def project(p):
    p=p.resolve();return p/'cube-dim-5-bounds' if (p/'cube-dim-5-bounds').is_dir() else p
def poly(p):return {tuple(map(int,k.split(','))):F(v) for k,v in p.items()}
def corner(p):return sum((value*len(set(permutations(k))) for k,value in p.items()),F(0))
def pi_upper():
    def atan(q):
        lo=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(20)),F(0))
        return lo,lo+F(1,41*q**41)
    a,b=atan(5),atan(239);return 16*a[1]-4*b[0]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',required=True,type=Path)
    parser.add_argument('--backend-root',required=True,type=Path)
    parser.add_argument('--run',action='store_true');args=parser.parse_args()
    root=project(args.project_root);backend=project(args.backend_root)
    locked={}
    for directory,expected in (('astra_upper006',M6),('astra_upper007',M7)):
        manifest=root/directory/'manifest.json';assert sha(manifest)==expected
        locked[str(manifest.relative_to(root))]=expected
        for name,info in json.loads(manifest.read_text())['files'].items():
            p=root/directory/name;assert sha(p)==info['sha256'] and p.stat().st_size==info['bytes']
            locked[str(p.relative_to(root))]=sha(p)
    header=backend/'upper/round012_graph_moments_gmp.hpp';assert sha(header)==HEADER
    locked['backend/upper/round012_graph_moments_gmp.hpp']=HEADER
    receipt6=json.loads((root/'astra_upper006/verification.json').read_text())
    receipt7=json.loads((root/'astra_upper007/verification.json').read_text())
    for receipt in (receipt6,receipt7):
        for legacy,expected in receipt['source_sha256'].items():
            p=root.joinpath(*Path(legacy).parts[-2:]);assert sha(p)==expected
            locked[str(p.relative_to(root))]=sha(p)
    if not args.run:
        print(json.dumps({'status':'locked_inputs_only','calculation_executed':False,'locked_count':len(locked)}));return
    began=time.monotonic()
    tau,s=F(receipt6['tau']),F(5,11);rho=1-s
    assert tau==F(15800608,26277615) and F(1,2)<tau<1
    taylor=[F(1)]
    for n in range(1,4):taylor.append(taylor[-1]*(rho+n-1)/n)
    assert taylor==[F(v) for v in receipt6['Taylor_u_coefficients']]
    assert taylor==[F(1),F(6,11),F(51,121),F(476,1331)]
    assert rho*(rho+1)*(rho+2)*(rho+3)>0
    assert 1-1/tau>-1 and taylor[3]/taylor[2]==F(28,33)<1
    TL,TU=F(receipt6['T_lower']),F(receipt6['T_upper'])
    assert 0<TL<TU and TL**11*tau**6<1<TU**11*tau**6
    normalized=[F(0)]*5
    for j,coef in enumerate(taylor):
        for k in range(j+1):normalized[k+1]+=coef*math.comb(j,k)*(-1/tau)**k
    assert normalized==[F(v) for v in receipt6['normalized_d_coefficients']]
    a,eta,c=(F(receipt6[k]) for k in ('a_selected','eta','c_lower'))
    assert a==F(679,1000) and eta==F(1,8) and a*(1-2*eta)>0
    assert c>0 and c*c*2**15*pi_upper()**5<9
    families={name:[poly(p) for p in values] for name,values in json.loads((root/'astra_upper007/space_moments.json').read_text())['families'].items()}
    assert set(families)=={'d1','d2','d3','d4','b'}
    for degree in range(1,5):
        assert len(families['d'+str(degree)])==8
        for n,p in enumerate(families['d'+str(degree)]):
            assert all(sum(k)==n and k==tuple(sorted(k,reverse=True)) for k in p)
            assert corner(p)==F(receipt6['moments'][str(degree)][n])
    assert len(families['b'])==12
    assert all(sum(k)==n+1 for n,p in enumerate(families['b']) for k in p)
    power={}
    for n in range(8):
        factor=F(11,4)*a**7*math.comb(7,n)*eta**n*TL
        for degree in range(1,5):
            for key,value in families['d'+str(degree)][n].items():power[key]=power.get(key,F(0))+factor*normalized[degree]*value
    for n in range(12):
        factor=-F(7,4)*a**11*math.comb(11,n)*eta**n
        for key,value in families['b'][n].items():power[key]=power.get(key,F(0))+factor*value
    power={k:v for k,v in power.items() if v}
    assert power==poly(receipt7['symmetric_power_coefficients']) and len(power)==127
    assert corner(power)==F(receipt6['corner_polynomial_value'])
    face={}
    for partition,value in power.items():
        for exponent in set(permutations(partition)):
            key=exponent[1:];face[key]=face.get(key,F(0))+value
    face={k:v for k,v in face.items() if v};assert len(face)==1316
    assert tuple(max(k[i] for k in face) for i in range(4))==(6,)*4
    D=math.lcm(*(v.denominator for v in face.values()));indices=list(product(range(7),repeat=4))
    tensor={k:0 for k in indices};L=math.lcm(*(math.comb(6,j) for j in range(7)))
    weights=[[math.comb(k,j)*(L//math.comb(6,j)) if j<=k else 0 for j in range(7)] for k in range(7)]
    for j,value in face.items():
        coeff=value.numerator*(D//value.denominator)
        for k in product(*(range(exponent,7) for exponent in j)):
            tensor[k]+=coeff*math.prod(weights[x][y] for x,y in zip(k,j))
    tensor={k:F(value,D*L**4) for k,value in tensor.items()}
    minimum=min(tensor.values());minimum_keys=[k for k,v in tensor.items() if v==minimum]
    assert minimum==F(receipt7['minimum_Bernstein_value']) and minimum>0
    assert c*minimum==F(receipt7['global_g_exact_rational'])==F(receipt6['exact_corner_g_lower'])
    assert c*minimum>F(5563,10**8)>F(4583,10**8)
    coefficient_hash=hashlib.sha256(json.dumps({','.join(map(str,k)):str(v) for k,v in sorted(tensor.items())},sort_keys=True).encode()).hexdigest()
    assert coefficient_hash==receipt7['Bernstein_coefficients_sha256']
    output={'status':'pass_K3_global_g_with_shared_d3_d4_space_moments',
            'g':'5563/100000000','eta':str(eta),'a':str(a),
            'minorant':{'tau':str(tau),'Taylor_coefficients':[str(v) for v in taylor],
                        'T_lower':str(TL),'T_upper':str(TU),'normalized_d_coefficients':[str(v) for v in normalized],
                        'common_factor_power_bracket_verified':True,'positive_remainder_and_difference_factor_checked':True},
            'Bernstein':{'degrees':[6]*4,'coefficient_count':len(tensor),'minimum_indices':minimum_keys,
                         'minimum':str(minimum),'c_lower_times_minimum':str(c*minimum),'coefficients_sha256':coefficient_hash},
            'dependency_limit':'All degree27 graph moments and d3/d4 spatial coefficients remain locked author inputs. This batch independently verifies constants, assembly and all 2401 Bernstein coefficients, not 147258 graph production.',
            'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
            'registration_sha256':sha(HERE/'registration.md'),
            'runtime':{'python':sys.executable,'seconds':time.monotonic()-began}}
    (HERE/'result.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ('status','Bernstein','runtime')},indent=2))

if __name__=='__main__':main()
