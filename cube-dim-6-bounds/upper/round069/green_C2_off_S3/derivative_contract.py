"""Exact off-true-S3 qualification and formal tuple tools; no body32 fan engine."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import comb, lcm
import hashlib
import json
import sympy as s

def rational(v):
    if isinstance(v,bool) or not isinstance(v,(str,int,F)):
        raise ValueError('Exact rational strings/integers/Fractions required')
    return F(v)

def sym(q):return s.Rational(q.numerator,q.denominator)
def frac(q):return F(str(q))
def canonical_hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def qualify_off_S3(columns,point):
    cols=[list(map(rational,col)) for col in columns]
    x=list(map(rational,point))
    if len(x)!=6 or len(cols)<6 or any(len(col)!=6 for col in cols):
        raise ValueError('A spanning six-dimensional frame and point are required')
    if any(not any(col) for col in cols):raise ValueError('Remove zero columns with a separate symbol proof')
    A=s.Matrix.hstack(*(s.Matrix(list(map(sym,col))) for col in cols))
    if A.rank()!=6:raise ValueError('Frame rank must be six')
    selected=[]
    for j in range(len(cols)):
        trial=selected+[j]
        if A[:,trial].rank()>len(selected):selected.append(j)
        if len(selected)==6:break
    N=A[:,selected];inverse=N.inv();xx=s.Matrix(list(map(sym,x)));norm2=(xx.T*xx)[0]
    rows=[]
    for triple in combinations(range(len(cols)),3):
        independent=[]
        for j in triple:
            if A[:,independent+[j]].rank()>len(independent):independent.append(j)
        U=A[:,independent];gram=U.T*U;v=U.T*xx
        residual=s.cancel(norm2-(v.T*gram.inv()*v)[0])
        rho=frac(residual)
        if rho<0:raise ValueError('Negative exact squared residual')
        rows.append({'triple':list(triple),'rank':len(independent),'independent_indices':independent,'squared_residual':str(rho)})
    assert len(rows)==comb(len(cols),3)
    rho=min(F(row['squared_residual']) for row in rows)
    canonical={'columns':[[str(v) for v in col] for col in cols]}
    C=frac(sum(v*v for v in inverse));c=1/(1+C)
    denom=lcm(*(v.denominator for col in cols for v in col))
    out={'status':'OFF_TRUE_S3' if rho>0 else 'ON_TRUE_S3_UNSUPPORTED',
         'frame_sha256':canonical_hash(canonical),'point':list(map(str,x)),
         'columns':canonical['columns'],'frame_rank':6,'basis_indices':selected,
         'basis_inverse':[[str(inverse[i,j]) for j in range(6)] for i in range(6)],
         'triple_count':len(rows),'expected_triple_count':comb(len(cols),3),
         'triples':rows,'min_squared_distance':str(rho),
         'delta_point':str(min(rho,F(1))/2) if rho>0 else None,
         'coercivity_c':str(c),'column_l1_sum':str(sum(abs(v) for col in cols for v in col)),
         'common_column_denominator':denom,'C2_on_S3_claim':False}
    return out

def local_tail_contract(qualification,ball_radius,epsilon):
    if qualification['status']!='OFF_TRUE_S3':raise ValueError('Full-G derivative contract requires strict off-true-S3 qualification')
    radius,eps=map(rational,(ball_radius,epsilon));delta=F(qualification['delta_point'])-radius
    if radius<0 or delta<=0 or eps<=0:raise ValueError('Require nonnegative radius, positive residual distance and epsilon')
    c=F(qualification['coercivity_c']);S=F(qualification['column_l1_sum']);R=delta/(2*S)
    m=len(qualification['columns']);d=qualification['common_column_denominator']
    M=F(comb(m,4)*d**4)/(64*c*c*R**8)
    inv=[list(map(F,row)) for row in qualification['basis_inverse']]
    col=[[inv[i][j] for i in range(6)] for j in range(6)]
    c1=[sum(map(abs,v)) for v in col]
    c2=[[c1[j]*c1[k]+2*sum(abs(a*b) for a,b in zip(col[j],col[k])) for k in range(6)] for j in range(6)]
    return {'status':'QUALIFIED_LOCAL_C2_TAIL_CONTRACT','frame_sha256':qualification['frame_sha256'],
            'center':qualification['point'],'ball_radius':str(radius),'uniform_distance_delta':str(delta),
            'epsilon':str(eps),'Mbar_delta':str(M),'C1_coordinate':list(map(str,c1)),
            'C2_coordinate':[[str(v) for v in row] for row in c2],
            'value_tail_upper':str(M*eps**3/3),
            'gradient_tail_absolute':list(map(str,[M*v*eps**2/2 for v in c1])),
            'Hessian_tail_absolute':[[str(M*v*eps) for v in row] for row in c2],
            'tail_object':'G-H_epsilon and matching derivatives; no fixed-cutoff homogeneity substitution',
            'complete_fan_proof_required_for_finite_part':True,'global_cover_certified':False}

def bind_derivative_request(qualification,fan_metadata,fields):
    if qualification['status']!='OFF_TRUE_S3':raise ValueError('Strict off-S3 qualification required')
    if fan_metadata.get('frame_sha256')!=qualification['frame_sha256']:
        raise ValueError('Fan/frame normalization mismatch')
    if fan_metadata.get('complete_fan_proved') is not True or fan_metadata.get('multiplicities_included_once') is not True:
        raise ValueError('Complete fan and multiplicity proof required')
    if not isinstance(fan_metadata.get('proof_reference'),str) or not fan_metadata['proof_reference']:
        raise ValueError('External fan proof reference required')
    if any(len(field)>2 or any(not isinstance(k,int) or isinstance(k,bool) or k<0 or k>=6 for k in field) for field in fields):
        raise ValueError('Only physical coordinate fields of orders0,1,2 are supported')
    return {'status':'QUALIFIED_OFF_S3_CONDITIONAL_ON_EXTERNAL_COMPLETE_FAN_PROOF',
            'frame_sha256':qualification['frame_sha256'],'point':qualification['point'],
            'fields':[list(field) for field in fields],'fan_proof_reference':fan_metadata['proof_reference'],
            'fan_proof_independently_authenticated_by_this_function':False,
            'constant_normalization':'1/(128*pi^5)','negative_channels_must_cancel_in_complete_physical_sum':[1,3,5,7],
            'zero_repeated_nodes':'exact signs/equality and confluent tuples retained',
            'G_C2_on_S3_claim':False,'global_cover_certified':False}

CHANNELS=('Gamma','L1','L3','L5','L7')
def formal_dd_channels(nodes):
    """Formal coefficient extraction only: no assertion of a Green derivative."""
    nodes=tuple(map(rational,nodes))
    if len(nodes) not in (6,7,8):raise ValueError('Only DD5,DD6,DD7 tuple lengths')
    counts=Counter(nodes);answer={key:F(0) for key in CHANNELS};weights=[]
    for b,m in sorted(counts.items()):
        coefficients=[F(1)]+[F(0)]*(m-1)
        for other,n in counts.items():
            if other==b:continue
            factor=[F((-1)**l*comb(n+l-1,l))/(b-other)**(n+l) for l in range(m)]
            coefficients=[sum(coefficients[j]*factor[l-j] for j in range(l+1)) for l in range(m)]
        for k in range(m):
            weight=coefficients[m-1-k];weights.append((b,k,weight))
            if k==0 and b:answer['Gamma']+=weight*(1 if b>0 else -1)
            if b==0 and k%2==1:answer['L'+str(k)]+=weight*F((-1)**((k-1)//2),k)
    return answer,weights

def formal_cone_field(base,jacobian,field):
    """One formal cone field; caller must separately prove qualification/full fan."""
    base=list(map(rational,base));J=[list(map(rational,row)) for row in jacobian]
    if len(base)!=6 or len(J)!=6 or any(len(row)!=6 for row in J):raise ValueError('Six labelled nodes/J rows required')
    total={key:F(0) for key in CHANNELS}
    def add(nodes,coefficient):
        if coefficient:
            values,_=formal_dd_channels(nodes)
            for key in CHANNELS:total[key]+=coefficient*values[key]
    if len(field)==0:add(base,F(1))
    elif len(field)==1:
        k=field[0]
        for i in range(6):add(base+[base[i]],J[i][k])
    elif len(field)==2:
        k,l=field
        for i in range(6):add(base+[base[i],base[i]],2*J[i][k]*J[i][l])
        for i in range(6):
            for j in range(i+1,6):add(base+[base[i],base[j]],J[i][k]*J[j][l]+J[j][k]*J[i][l])
    else:raise ValueError('Order above2 unsupported')
    return total
