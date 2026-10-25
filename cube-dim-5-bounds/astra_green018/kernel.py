"""Full generic Green DD4 kernel on the locked3504 original-vertex simplices."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from math import lcm
import gzip,json,hashlib,time
from flint import arb,fmpq,ctx
P=Path(__file__).resolve().parent
TABLE=P.parent/'astra_vertex017'/'boundary_triangulation.json.gz'
TABLE_SHA='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def ar(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))
def endpoint(z,upper):return F(str((z.upper() if upper else z.lower()).fmpq()))
def dd4(nodes,logs=None):
    """Sorted confluent divided difference, nonzero integer nodes only."""
    if any(n==0 for n in nodes):raise ValueError('zero log node is unsupported')
    nodes=sorted(nodes);table=[]
    for right,node in enumerate(nodes):
        row=[arb(abs(node)).log() if logs is None else logs[abs(node)]]
        for order in range(1,right+1):
            left=nodes[right-order]
            row.append(ar(F((-1)**(order-1),order*node**order)) if node==left
                       else (row[order-1]-table[right-1][order-1])/(node-left))
        table.append(row)
    return table[-1][-1]
class GreenKernel:
    def __init__(self):
        assert sha(TABLE)==TABLE_SHA
        data=json.loads(gzip.decompress(TABLE.read_bytes()));assert data['denominator']==144
        self.vertices=data['vertices'];self.simplices=data['simplices']
        assert len(self.vertices)==280 and len(self.simplices)==3504
    def projections(self,point,distinct=False):
        point=tuple(map(F,point));D=lcm(*(x.denominator for x in point))
        n=[int(x*D) for x in point]
        k=[sum(a*b for a,b in zip(n,z)) for z in self.vertices]
        if 0 in k:raise ValueError('zero projection on a true vertex; no extension attempted')
        if distinct and len(set(k))!=280:raise ValueError('repeated true-vertex projections in strict-generic batch')
        return D,k
    def evaluate(self,point,distinct=False):
        began=time.monotonic();D,k=self.projections(point,distinct)
        logs={abs(n):arb(abs(n)).log() for n in k};total=arb(0);absolute_sum=arb(0);patterns=Counter()
        for s in self.simplices:
            nodes=tuple(k[i] for i in s['vertices']);patterns[tuple(sorted(Counter(nodes).values(),reverse=True))]+=1
            term=s['det_numerator']*dd4(nodes,logs);total-=term;absolute_sum+=abs(term)
        factor=arb(D)**4/(144*(2*arb.pi())**5);value=total*factor
        assert value.is_finite() and value>0
        ordered=sorted(set(k));lower,upper=endpoint(value,False),endpoint(value,True)
        return {'point':[str(F(x)) for x in point],'common_denominator':D,'precision_bits':ctx.prec,
                'G_ball':str(value),'G_lower_dyadic':str(lower),'G_upper_dyadic':str(upper),
                'projection_count':len(k),'distinct_projection_count':len(set(k)),
                'minimum_abs_projection':str(F(min(map(abs,k)),144*D)),
                'minimum_distinct_projection_gap':str(F(min(b-a for a,b in zip(ordered,ordered[1:])),144*D)),
                'multiplicity_patterns':{','.join(map(str,key)):v for key,v in patterns.items()},
                'absolute_term_sum':str(absolute_sum*factor),'cancellation_ratio':str(absolute_sum*factor/value),
                'seconds':time.monotonic()-began}
