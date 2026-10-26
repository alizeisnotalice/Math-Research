"""Independent integer closed-corner check for frozen EPS128 target110 jets."""
from fractions import Fraction as F
from itertools import product
from math import lcm
from interfaces import common as C,arith as A


def certificate(node,cert,record,roots):
    C.require(record['epsilon']==cert['leaf_epsilon']=='1/128' and cert['target']=='11/100000' and cert['lower_bound_object']=='G' and cert['component_function']=='H_epsilon','certificate epsilon/model binding differs')
    C.require(cert['mechanism']=='075_full_frame_16' and cert['quadratic']=='I4/(16epsilon^2)' and F(cert['penalty_coefficient'])==1024 and cert['checked_points']==16,'EPS128 curvature differs')
    C.require(cert['point_key']==','.join(record['point']),'certificate anchor differs')
    values=[list(map(F,pair)) for pair in record['values']];C.require(len(values)==5 and all(a<=b for a,b in values) and values[0][0]>F(11,100000),'H128/target interval fails')
    for lo,hi in values:C.require((hi-lo)*2**40<=max(abs(lo),abs(hi),F(1,10**6)),'actual production width condition fails')
    gradients=[list(map(F,pair)) for pair in cert['loggradient_directed']];C.require(len(gradients)==4 and all(lo<=hi for lo,hi in gradients),'gradient interval order differs')
    for (gl,gu),(dl,du) in zip(gradients,values[1:]):
        for h in values[0]:
            for d in (dl,du):C.require(gl*h<=d<=gu*h,'same-epsilon physical derivative/H containment fails')
    delta=F(cert['Delta_lower']);ratio=values[0][0]/F(11,100000);lower,upper=A.log_bounds(ratio.numerator,ratio.denominator)
    C.require(0<delta<=F(lower,A.S),'fresh H128 Delta not below independent integer log')
    lo,hi=C.physical_box(roots,node);point=list(map(F,record['point']))
    C.require(point[0]==1 and all(a<=x<=b for a,x,b in zip(lo,point[1:],hi)),'anchor outside independently constructed closed box')
    offsets=[(a-x,b-x) for a,x,b in zip(lo,point[1:],hi)];vd=lcm(*(x.denominator for pair in offsets for x in pair));gd=lcm(delta.denominator,*(x.denominator for pair in gradients for x in pair))
    zz=[tuple(int(x*vd) for x in pair) for pair in offsets];gg=[tuple(int(x*gd) for x in pair) for pair in gradients];di=int(delta*gd);rows=[]
    for bits in product((0,1),repeat=4):
        v=[zz[j][bits[j]] for j in range(4)];num=di*vd*vd+vd*sum(min(a*x,b*x) for (a,b),x in zip(gg,v))-1024*gd*sum(x*x for x in v);rows.append((num,bits))
    minimum,bits=min(rows);margin=F(minimum,gd*vd*vd)
    C.require(cert['passed']==(minimum>0) and F(cert['minimum_log_margin'])==margin and tuple(cert['worst_corner_bits'])==bits,'all16 exact H128 corners do not strictly pass')
    return {'node':list(node),'margin':str(margin),'certificate_sha256':C.digest(cert),'record_sha256':C.digest(record),'closed_corners':16,'penalty':1024,'passed':minimum>0}

