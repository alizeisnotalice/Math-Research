"""Independent exact fmpq Hermite DD6 reference; never uses barycentric AD."""
from pathlib import Path
from fractions import Fraction as Q
from flint import fmpq
import itertools,json,hashlib,time,sys
HERE=Path(__file__).resolve().parent
sys.set_int_max_str_digits(0)

def sign_dd(nodes):
    z=sorted(nodes)
    if any(not v for v in z):raise ValueError('unsupported exact zero node')
    b=[fmpq(1 if v>0 else -1)for v in z]
    for order in range(1,len(z)):
        b=[fmpq(0)if z[j]==z[j+order]else(b[j+1]-b[j])/(z[j+order]-z[j])for j in range(len(z)-order)]
    return b[0]

def run_gradient():
    started=time.perf_counter()
    source=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel')
    fan=source/'fan_integer_input.txt';old=source/'local_minimum/exact_local_reference_results.json'
    ref=json.loads(old.read_text());x=ref['point_numerators'];den=ref['denominator'];tokens=iter(fan.read_text().split());cones=[]
    for _ in range(int(next(tokens))):
        det=int(next(tokens));rays=[];alpha=fmpq(det)
        for j in range(6):
            q=int(next(tokens));r=tuple(int(next(tokens))for _ in range(6));rays.append((q,r));alpha/=q
        cones.append((alpha,rays))
    total=fmpq(0);terms=0;constant_terms=0;copied_DD6=0
    for perm in itertools.permutations(range(6)):
        pos=perm.index(1)
        for mask in range(64):
            signs=tuple(-1 if mask>>i&1 else 1 for i in range(6));y=tuple(signs[i]*x[perm[i]]for i in range(6))
            for alpha,rays in cones:
                nodes=[fmpq(sum(y[i]*r[i]for i in range(6)),den*q)for q,r in rays]
                if any(not b for b in nodes):raise ValueError('entire reference unsupported: zero node')
                terms+=1
                if all(b>0 for b in nodes)or all(b<0 for b in nodes):constant_terms+=1;continue
                local=fmpq(0)
                for j,(q,r)in enumerate(rays):
                    if r[pos]:
                        local+=fmpq(signs[pos]*r[pos],q)*sign_dd(nodes+[nodes[j]])
                        copied_DD6+=1
                total+=alpha*local
            if time.perf_counter()-started>600:raise TimeoutError('independent reference600s limit')
    assert terms==967680
    q=total/128
    result={'status':'PASS_EXACT_FULL_GROUP_DD6_GRADIENT_X2','point_numerators':x,'denominator':den,
            'component_index':1,'gradient_Q_over_pi5_numerator':str(q.numerator),'gradient_Q_over_pi5_denominator':str(q.denominator),
            'full_group_terms':terms,'all_same_nonzero_sign_terms_exact_zero':constant_terms,'copied_node_DD6_calls':copied_DD6,
            'global_negative_pairing_used':False,'elapsed_seconds':time.perf_counter()-started,
            'fan_sha256':hashlib.sha256(fan.read_bytes()).hexdigest(),'old_G_certificate_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'gradient_exact_reference.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if not k.startswith('gradient_Q')},indent=2))

if __name__=='__main__':run_gradient()
