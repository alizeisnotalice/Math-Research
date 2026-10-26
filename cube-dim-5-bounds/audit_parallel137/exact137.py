"""New certificate-wide integer clearing path, independent of audit130 helpers."""
from math import gcd


def f(s):
    a = str(s).split('/')
    assert len(a) in (1,2)
    n,d = int(a[0]),int(a[1]) if len(a)==2 else 1
    assert d > 0
    return n,d


def eq(a,b):
    n,d=f(a);m,e=f(b)
    return n*e == m*d


def le(a,b):
    n,d=f(a);m,e=f(b)
    return n*e <= m*d


def reduced(n,d):
    k=gcd(n,d);n//=k;d//=k
    return str(n) if d==1 else f'{n}/{d}'


def key(node):
    rid,path=node
    assert type(rid) is int and 0<=rid<52360
    assert isinstance(path,str) and not set(path)-{'0','1'}
    return rid,path


def canonical_roots():
    return [[a,b,c,d] for a in range(32) for b in range(a+1)
            for c in range(b+1) for d in range(c+1)]


def box(node,roots):
    rid,path=key(node)
    low=list(roots[rid]);den=[32]*4
    for i,bit in enumerate(path):
        j=i%4;low[j]=2*low[j]+int(bit);den[j]*=2
    return [(low[j],low[j]+1,den[j]) for j in range(4)]


def in_box(point,node,roots):
    assert len(point)==5 and eq(point[0],1)
    return all(lo*pd<=pn*d<=hi*pd
               for (lo,hi,d),(pn,pd) in zip(box(node,roots),map(f,point[1:])))


def check_child(child,roots):
    bounds=box(child['node'],roots)
    assert len(child['box_lower'])==len(child['box_upper'])==4
    for (lo,hi,d),a,b in zip(bounds,child['box_lower'],child['box_upper']):
        an,ad=f(a);bn,bd=f(b)
        assert lo*ad==an*d and hi*bd==bn*d


def new_q_check(source):
    c=source['core'];qn,qd=f(c['penalty_coefficient'])
    assert c['target']=='11/100000' and 0<qn<=1024*qd
    hn,hd=f(source['H_interval'][0]);un,ud=f(source['H_interval'][1])
    assert hn>0 and hn*ud<=un*hd
    if c['epsilon']=='1/64':
        assert c['curvature_rule']=='075_H64_Q256' and qn==256*qd
        assert source['H64_g105_ancestry'] is None
    else:
        assert c['epsilon']=='1/128'
        if c['curvature_rule']=='129_anchor_Jensen_g105_U128':
            a=source['H64_g105_ancestry'];assert a
            assert a['H64_component_epsilon']=='1/64' and eq(a['H64_lower'],'21/200000')
            assert eq(a['H128_upper'],source['H_interval'][1])
            assert 21*ud<=200000*un and 256*qd<=qn<1024*qd
            assert qn*(200000*un)==qd*(1024*200000*un-768*21*ud)
        else:
            assert c['curvature_rule']=='075_H128_Q1024_no_H64g105_ancestor'
            assert qn==1024*qd and source['H64_g105_ancestry'] is None


def positive_certificate(source,child,certificate,roots):
    """One common denominator for ALL terms of ALL16 corners; no pair algebra."""
    c=source['core'];check_child(child,roots);new_q_check(source)
    for name in ('epsilon','point','target','Delta_lower','loggradient_directed',
                 'curvature_rule','penalty_coefficient'):
        assert certificate[name]==c[name]
    assert certificate['child_node']==child['node']
    assert certificate['source_id']==source['source_id']
    assert certificate['checked_points']==16 and certificate['passed'] is True
    point=list(map(f,c['point'][1:]));bounds=box(child['node'],roots)
    W=1
    for _,d in point:W=W//gcd(W,d)*d
    for lo,hi,d in bounds:W=W//gcd(W,d)*d
    integer_anchor=[n*(W//d) for n,d in point]
    integer_bounds=[(lo*(W//d),hi*(W//d)) for lo,hi,d in bounds]
    dn,dd=f(c['Delta_lower']);qn,qd=f(c['penalty_coefficient'])
    gradient=[tuple(map(f,pair)) for pair in c['loggradient_directed']]
    assert len(gradient)==4
    T=dd
    for (an,ad),(bn,bd) in gradient:
        assert an*bd<=bn*ad
        T=T//gcd(T,ad*W)*(ad*W);T=T//gcd(T,bd*W)*(bd*W)
    T=T//gcd(T,qd*W*W)*(qd*W*W)
    delta_int=dn*(T//dd)
    gradient_int=[(an*(T//(ad*W)),bn*(T//(bd*W))) for (an,ad),(bn,bd) in gradient]
    quadratic_int=qn*(T//(qd*W*W))
    all16=[]
    for mask in range(16):
        bits=[(mask>>(3-j))&1 for j in range(4)]
        u=[integer_bounds[j][bits[j]]-integer_anchor[j] for j in range(4)]
        m=delta_int+sum(min(a*z,b*z) for (a,b),z in zip(gradient_int,u))-quadratic_int*sum(z*z for z in u)
        all16.append((m,bits))
    minimum,bits=min(all16)
    assert minimum>0
    mn,md=f(certificate['minimum_log_margin'])
    assert minimum*md==mn*T
    assert certificate['worst_corner_bits']==bits
    assert certificate['anchor_inside_child']==in_box(c['point'],child['node'],roots)
    return {'minimum_margin':reduced(minimum,T),'worst_corner_bits':bits,
            'all16_integer_numerators_positive':True,'clearing_denominator_bits':T.bit_length()}


def full_tree(leaves,splits):
    """Heap-coded full binary forest; independent of130 interval-tiling route."""
    def heap(n):
        r,p=key(n);return r,(1<<len(p))+(int(p,2) if p else 0)
    L=set(map(heap,leaves));I=set(map(heap,splits))
    assert len(L)==len(leaves) and len(I)==len(splits) and not L&I
    allnodes=L|I
    assert all((r,1) in allnodes for r in range(52360))
    for r,k in allnodes:
        if k>1:assert (r,k//2) in I
    for r,k in I:assert (r,2*k) in allnodes and (r,2*k+1) in allnodes
    by=[[] for _ in range(52360)]
    for r,k in L:by[r].append(k.bit_length()-1)
    for depths in by:
        assert depths
        D=max(depths)
        assert sum(1<<(D-d) for d in depths)==1<<D
    assert len(L)==52360+len(I)
    return {'roots':52360,'leaves':len(L),'splits':len(I),
            'heap_parent_and_both_child_closure':True,'antichain':True,'every_root_Kraft_exact_one':True}
