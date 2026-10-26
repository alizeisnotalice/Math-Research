"""C115 same512 native physics; explicit general0<a<T window schema."""
import math
from fractions import Fraction as Q
def require(x,s):
    if not x:raise ValueError(s)

def validate_native(m,expected):
    # Same accepted pointjet physics; Euler now uses the actual positive window.
    w=expected['window']
    require(w['kind']=='B' and 0<Q(w['a'])<Q(w['T']),'General positive band window')
    require(m.get('status')=='all_real_nodes_supported'and m.get('precision_bits')==512,'Native status/precision')
    require(m.get('full_group_terms')==967680 and m.get('computed_sign_pairs')==483840,'Incomplete body32 group')
    for k in ['confluent_terms','zero_node_terms']:
        v=m.get(k);require(type(v)is int and 0<=v<=967680 and v%2==0,'Native count '+k)
    require(type(m.get('unique_node_values'))is int and 1<=m['unique_node_values']<=2903040,'Distinct-node count')
    require(isinstance(m.get('kernel_seconds'),(int,float))and not isinstance(m['kernel_seconds'],bool)
        and math.isfinite(m['kernel_seconds'])and m['kernel_seconds']>=0,'Native elapsed metadata')
    require(isinstance(m.get('gradient'),list)and len(m['gradient'])==6,'Six physical derivatives')
    F=M.dyadic(m['band']);g=[M.dyadic(d)for d in m['gradient']]
    pa=M.dyadic(m['p_a']);pt=M.dyadic(m['p_T']);ne=M.dyadic(m['Euler_xgrad_plus5B_plusTpT_minusapa'])
    require(F[0]>0 and pa[0]>0 and pt[0]>0,'Positive band/endpoints')
    require(ne[0]<=0<=ne[1],'Native Euler')
    def scale(i,q):return tuple(sorted(v*q for v in i))
    terms=[scale(F,Q(5)),scale(pt,Q(expected['window']['T'])),scale(pa,-Q(expected['window']['a']))]
    terms += [scale(i,Q(x))for i,x in zip(g,expected['center'])]
    e=tuple(sum(t[k]for t in terms)for k in [0,1]);require(e[0]<=0<=e[1],'Exact-rational Euler')
    return F,g,pa,pt,ne,e
