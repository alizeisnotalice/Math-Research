"""Fixed portable g/optimization eligibility audit; --run is explicit."""
import argparse,hashlib,json,math,sys,time
from fractions import Fraction as F
from itertools import permutations,product
from pathlib import Path

HERE=Path(__file__).resolve().parent
RUN='audit_g015_fixed_001'
LOCKS={'astra_upper009':'09f09e048c0d3b4f8a1c40e7341159b5a1baf7e7c3e479af47a86136da675896',
       'astra_upper006':'398f33d8565a30a68b9d4afbb5d29e1047a2a3a37783ab91a56127437bece136',
       'astra_upper007':'b8a7b09d6327428facdc8e6bc7524a87e226d42f2bf859e41240a42392672b5a'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):return {tuple(map(int,k.split(','))):F(v) for k,v in p.items()}
def trim(p):
    while len(p)>1 and not p[-1]:p.pop()
    return p
def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)
def der(p):return [i*p[i] for i in range(1,len(p))] or [F(0)]
def evaluate(p,x):
    out=F(0)
    for a in reversed(p):out=out*x+a
    return out
def primitive(p):
    p=trim(list(p));den=math.lcm(*(q.denominator for q in p))
    ints=[q.numerator*(den//q.denominator) for q in p]
    gcd=math.gcd(*ints);return [F(q,gcd) for q in ints]
def sign_variations(values):
    s=[1 if x>0 else -1 for x in values if x]
    return sum(a!=b for a,b in zip(s,s[1:]))
def descartes(p,a,b):
    # (1+t)^N p((a+b*t)/(1+t)): t>0 corresponds to x in (a,b).
    n=len(p)-1;out=[F(0)]*(n+1)
    for j,c in enumerate(p):
        for h in range(j+1):
            factor=c*math.comb(j,h)*a**(j-h)*b**h
            for k in range(n-j+1):out[h+k]+=factor*math.comb(n-j,k)
    return sign_variations(out)
def sturm(p):
    seq=[primitive(p),primitive(der(p))]
    while len(seq[-1])>1:
        rem=list(seq[-2]);div=seq[-1]
        while len(rem)>=len(div) and any(rem):
            degree=len(rem)-len(div);scale=rem[-1]/div[-1]
            for j,c in enumerate(div):rem[j+degree]-=scale*c
            trim(rem)
        if not any(rem):break
        seq.append(primitive([-x for x in rem]))
    return seq
def bernstein_range(p,a,b):
    if a==b:
        value=evaluate(p,a);return value,value
    n=len(p)-1
    translated=[sum((p[k]*math.comb(k,j)*a**(k-j)*(b-a)**j for k in range(j,n+1)),F(0)) for j in range(n+1)]
    coeff=[sum((translated[j]*F(math.comb(k,j),math.comb(n,j)) for j in range(k+1)),F(0)) for k in range(n+1)]
    return min(coeff),max(coeff)
def score_range(U,V,a,b):
    ul,uh=bernstein_range(U,a,b);vl,vh=bernstein_range(V,a,b)
    assert ul>0 and vl>0
    return ul**11/vh**7,uh**11/vl**7
def pi_bounds():
    def atan(q):
        lo=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(20)),F(0))
        return lo,lo+F(1,41*q**41)
    a,b=atan(5),atan(239);return 16*a[0]-4*b[1],16*a[1]-4*b[0]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',required=True,type=Path)
    parser.add_argument('--run',action='store_true');args=parser.parse_args()
    root=args.project_root.resolve()
    if (root/'cube-dim-5-bounds').is_dir():root=root/'cube-dim-5-bounds'
    locked={}
    for directory,expected in LOCKS.items():
        manifest=root/directory/'manifest.json';assert sha(manifest)==expected
        locked[str(manifest.relative_to(root))]=expected
        for name,info in json.loads(manifest.read_text())['files'].items():
            p=root/directory/name;assert sha(p)==info['sha256'] and p.stat().st_size==info['bytes']
            locked[str(p.relative_to(root))]=sha(p)
    receipt=json.loads((root/'astra_upper009/verification.json').read_text())
    for legacy,expected in receipt['source_sha256'].items():
        p=root.joinpath(*Path(legacy).parts[-2:]);assert sha(p)==expected
        locked[str(p.relative_to(root))]=sha(p)
    if not args.run:
        print(json.dumps({'run_handle':RUN,'status':'locked_inputs_only','calculation_executed':False}));return
    began=time.monotonic()
    corner=json.loads((root/'astra_upper006/verification.json').read_text())
    old=json.loads((root/'astra_upper004/corner_checkpoint.json').read_text())
    constants=json.loads((root/'astra_upper003/verification.json').read_text())['constants']
    normalized=list(map(F,corner['normalized_d_coefficients']))
    dm={k:list(map(F,corner['moments'][str(k)])) for k in range(1,5)}
    U=[math.comb(7,n)*sum((normalized[k]*dm[k][n] for k in range(1,5)),F(0)) for n in range(8)]
    V=[math.comb(11,n)*F(old['moments']['E_b_corner_W_power'][n]) for n in range(12)]
    assert U==list(map(F,receipt['normalized_U_power_coefficients'])) and V==list(map(F,receipt['V_power_coefficients']))
    left,right=mul(der(U),V),mul(U,der(V))
    D=primitive([11*a-7*b for a,b in zip(left,right)])
    assert len(D)==17 and D==list(map(F,receipt['primitive_derivative_coefficients']))
    lo,hi=map(F,receipt['stationary_root_intervals'][0]);zero,half=F(0),F(1,2)
    assert zero<lo<hi<half and all(evaluate(D,x)!=0 for x in (zero,lo,hi,half))
    variations=[descartes(D,a,b) for a,b in ((zero,half),(zero,lo),(lo,hi),(hi,half))]
    assert evaluate(D,lo)*evaluate(D,hi)<0
    fallback=None
    if variations[0]==1 or variations[1:]==[0,1,0]:
        root_method='Mobius_Descartes_plus_exact_endpoint_sign_change'
    else:
        seq=sturm(D)
        counts=[sign_variations([evaluate(p,a) for p in seq])-sign_variations([evaluate(p,b) for p in seq]) for a,b in ((zero,half),(lo,hi))]
        assert counts==[1,1];fallback={'degrees':[len(p)-1 for p in seq],'counts':counts}
        root_method='exact_Fraction_Sturm_fallback'
    critical=[(zero,zero),(lo,hi),(half,half)]
    scores=[score_range(U,V,a,b) for a,b in critical]
    assert scores[1][0]>max(scores[0][1],scores[2][1])
    TL,TU=F(corner['T_lower']),F(corner['T_upper']);tau=F(corner['tau'])
    assert TL>0 and TL**11*tau**6<1<TU**11*tau**6
    cL,cU=F(constants['c_lower']),F(constants['c_upper']);pil,pih=pi_bounds()
    assert cL>0 and cL*cL*2**15*pih**5<9<cU*cU*2**15*pil**5
    upper4=cU**4*TU**11*scores[1][1]
    assert upper4<F(5643,10**8)**4
    eta,a=F(27,250),F(69,100)
    assert lo*1000//1==hi*1000//1==108 and eta==F(108,1000)
    assert a**4<=TL*evaluate(U,eta)/evaluate(V,eta)<(a+F(1,1000))**4
    target=next(row for row in receipt['candidates'] if row['eta']=='27/250')
    families={name:[poly(p) for p in values] for name,values in json.loads((root/'astra_upper007/space_moments.json').read_text())['families'].items()}
    power={}
    for n in range(8):
        factor=F(11,4)*a**7*math.comb(7,n)*eta**n*TL
        for k in range(1,5):
            for key,value in families['d'+str(k)][n].items():power[key]=power.get(key,F(0))+factor*normalized[k]*value
    for n in range(12):
        factor=-F(7,4)*a**11*math.comb(11,n)*eta**n
        for key,value in families['b'][n].items():power[key]=power.get(key,F(0))+factor*value
    power={k:v for k,v in power.items() if v}
    assert power==poly(target['symmetric_power_coefficients']) and len(power)==127
    face={}
    for partition,value in power.items():
        for exponent in set(permutations(partition)):
            key=exponent[1:];face[key]=face.get(key,F(0))+value
    face={k:v for k,v in face.items() if v};assert tuple(max(k[i] for k in face) for i in range(4))==(6,)*4
    den=math.lcm(*(v.denominator for v in face.values()));indices=list(product(range(7),repeat=4));tensor={k:0 for k in indices}
    weights=[[math.comb(k,j)*(60//math.comb(6,j)) if j<=k else 0 for j in range(7)] for k in range(7)]
    for j,value in face.items():
        numerator=value.numerator*(den//value.denominator)
        for k in product(*(range(x,7) for x in j)):tensor[k]+=numerator*math.prod(weights[x][y] for x,y in zip(k,j))
    tensor={k:F(v,den*60**4) for k,v in tensor.items()};minimum=min(tensor.values())
    minimum_keys=[k for k,v in tensor.items() if v==minimum]
    assert minimum==F(target['minimum_Bernstein_value']) and minimum>0
    g=cL*minimum;assert g==F(target['global_g_exact_rational']) and g>F(2821,50000000)
    tensor_hash=hashlib.sha256(json.dumps({','.join(map(str,k)):str(v) for k,v in sorted(tensor.items())},sort_keys=True).encode()).hexdigest()
    assert tensor_hash==target['Bernstein_coefficients_sha256']
    out={'run_handle':RUN,'status':'pass_new_global_g_and_fixed_K3_linear_family_ceiling',
         'global_certificate':{'eta':str(eta),'a':str(a),'g':'2821/50000000','coefficient_count':len(tensor),
                               'minimum_indices':minimum_keys,'minimum':str(minimum),'exact_g':str(g),'coefficient_sha256':tensor_hash},
         'family_optimization':{'domain':'eta in [0,1/2], a>=0, fixed ideal K3 minorant only',
                                'root_method':root_method,'Descartes_variations_full_left_root_right':variations,
                                'Sturm_fallback':fallback,'root_interval':[str(lo),str(hi)],
                                'U_and_V_bounds_method':'exact Bernstein range on each critical interval; no Horner or Arb',
                                'score_intervals':[[str(x),str(y)] for x,y in scores],
                                'cU4_TU11_score_upper':str(upper4),'strict_family_ceiling':'5643/100000000',
                                'unique_global_maximum_root':True,'not_an_upper_bound_for_G':True},
         'shared_dependency':'All mixed graph/spatial moments remain author locked inputs; no 007/147258-moment replay. Only strongest rational candidate full-face certificate checked.',
         'locked_inputs_sha256':locked,'source_sha256':sha(Path(__file__)),
         'registration_sha256':sha(HERE/'registration.md'),
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began}}
    (HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'run_handle':RUN,'status':out['status'],'root_method':root_method,'Descartes_variations':variations,
                      'global_coefficient_sha256':tensor_hash,'runtime':out['runtime']},indent=2))

if __name__=='__main__':main()
