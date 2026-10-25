"""Exact global zero-node principal parts before the fullGreen epsilon limit."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
from collections import Counter
from math import lcm,gcd,comb
import gzip,json,hashlib,time,os
from flint import arb,fmpq,ctx
P=Path(__file__).resolve().parent
TABLE=P.parent/'astra_vertex017'/'boundary_triangulation.json.gz'
TABLE_SHA='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def ar(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))
def canonical(obj):return json.dumps(obj,sort_keys=True,separators=(',',':')).encode()
def save(path,obj):
    path=Path(path);tmp=Path(str(path)+'.partial');tmp.write_bytes(canonical(obj));os.replace(tmp,path)
def rank(rows):
    a=[list(map(int,r)) for r in rows];r=0
    for col in range(5):
        p=next((j for j in range(r,len(a)) if a[j][col]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];pivot=a[r][col]
        for j in range(r+1,len(a)):
            if a[j][col]:
                z=[pivot*x-a[j][col]*y for x,y in zip(a[j],a[r])];g=0
                for x in z:g=gcd(g,abs(x))
                a[j]=[x//g for x in z] if g else z
        r+=1
        if r==len(a):break
    return r
def functional(nodes):
    """DD4 as rational coefficients of independent Hermite data (node,j)."""
    nodes=sorted(nodes);table=[]
    for right,node in enumerate(nodes):
        row=[{(node,0):F(1)}]
        for order in range(1,right+1):
            left=nodes[right-order]
            if node==left:out={(node,order):F(1)}
            else:
                out={}
                for key,value in row[order-1].items():out[key]=out.get(key,F(0))+value/F(node-left)
                for key,value in table[right-1][order-1].items():out[key]=out.get(key,F(0))-value/F(node-left)
                out={key:value for key,value in out.items() if value}
            row.append(out)
        table.append(row)
    return table[-1][-1]

class ZeroGreenKernel:
    def __init__(self):
        assert sha(TABLE)==TABLE_SHA
        self.table=json.loads(gzip.decompress(TABLE.read_bytes()))
        assert self.table['denominator']==144 and len(self.table['vertices'])==280 and len(self.table['simplices'])==3504
    def evaluate(self,point,checkpoint,seconds=120):
        began=time.monotonic();point=tuple(map(F,point));L=lcm(*(x.denominator for x in point));n=[int(x*L) for x in point]
        V=[(1,)+s for s in product((-1,1),repeat=4)]
        pair_ranks=[rank([v,w,n]) for v,w in combinations(V,2)];assert len(pair_ranks)==120
        common={'point':[str(x) for x in point],'common_denominator':L,'integer_numerator':n,
                'body_pair_count':120,'pair_rank_counts':dict(sorted(Counter(pair_ranks).items())),
                'source_hashes':{str(x):sha(x) for x in [Path(__file__),TABLE]},'pid':os.getpid()}
        if min(pair_ranks)<3:return dict(common,status='STOP_point_in_S2',seconds=time.monotonic()-began)
        projections=[sum(a*b for a,b in zip(n,z)) for z in self.table['vertices']]
        patterns=Counter();weights=Counter();total_det=0
        for simplex in self.table['simplices']:
            nodes=tuple(sorted(projections[j] for j in simplex['vertices']));D=simplex['det_numerator'];weights[nodes]+=D;total_det+=D
            patterns[tuple(sorted(Counter(nodes).values(),reverse=True))]+=1
        assert 0 in projections # this batch is specifically the zero-ray extension
        total={};done=0
        for nodes,D in sorted(weights.items()):
            for key,value in functional(nodes).items():total[key]=total.get(key,F(0))+D*value
            done+=1
            if time.monotonic()-began>seconds:
                state=dict(common,status='STOP_120_second_checkpoint',completed_signatures=done,total_signatures=len(weights),
                           coefficients=[[k,j,str(v)] for (k,j),v in sorted(total.items())],seconds=time.monotonic()-began)
                save(checkpoint,state);return state
        # Exact polynomial reproduction, before using even-function symmetry.
        reproduction=[]
        for power in range(5):
            value=sum((c*comb(power,j)*k**(power-j) for (k,j),c in total.items() if j<=power),F(0))
            assert value==(total_det if power==4 else 0);reproduction.append(str(value))
        zero=[total.get((0,j),F(0)) for j in range(5)];C=144*L
        principal={'log_epsilon':zero[0],'epsilon_minus2':zero[2]/(2*C**2),'epsilon_minus4':-zero[4]/(4*C**4)}
        common.update(zero_projection_vertices=projections.count(0),distinct_projection_count=len(set(projections)),
                      signature_count=len(weights),multiplicity_patterns={','.join(map(str,k)):v for k,v in patterns.items()},
                      polynomial_reproduction_degrees_0_to_4=reproduction,
                      zero_Hermite_coefficients_orders_0_to_4=[str(v) for v in zero],
                      principal_coefficients_inside_common_negative_factor={k:str(v) for k,v in principal.items()})
        if any(principal.values()):
            out=dict(common,status='STOP_nonzero_global_principal_part',raw_coefficients=[[k,j,str(v)] for (k,j),v in sorted(total.items()) if v],seconds=time.monotonic()-began)
            save(checkpoint,out);return out
        # For every epsilon>0 h_epsilon is even; odd derivatives have odd parity.
        folded={}
        for (k,j),c in total.items():
            if k:
                key=(abs(k),j);folded[key]=folded.get(key,F(0))+c*((-1)**j if k<0 else 1)
        folded={key:value for key,value in folded.items() if value}
        assert sum((c for (k,j),c in folded.items() if j==0),F(0))==0
        rational=sum((c*F((-1)**(j-1),j*k**j) for (k,j),c in folded.items() if j),F(0))
        remainder=ar(rational)
        for (k,j),c in folded.items():
            if j==0:remainder+=ar(c)*arb(k).log()
        factor=arb(L)**4/(144*(2*arb.pi())**5);G=-factor*remainder
        assert G.is_finite() and G>0
        coefficients=[[k,j,str(v)] for (k,j),v in sorted(folded.items())]
        out=dict(common,status='PASS_exact_zero_principal_parts_and_full_Green',precision_bits=ctx.prec,
                 nonzero_folded_coefficients=coefficients,nonzero_coefficients_sha256=hashlib.sha256(canonical(coefficients)).hexdigest(),
                 rational_remainder=str(rational),nonzero_remainder_ball=str(remainder),G_ball=str(G),
                 G_lower_dyadic=str(G.lower().fmpq()),G_upper_dyadic=str(G.upper().fmpq()),
                 odd_zero_data_identically_zero_for_all_epsilon=True,seconds=time.monotonic()-began)
        save(checkpoint,out);return out
