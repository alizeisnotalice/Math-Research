"""Independent rational replay: stationary signs and Bernstein polynomial bounds."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json
import sys

def stationary(d, q):
    return 6*q**3-3*(3+d)*q*q+3*(1+d)*q-d

def moment(d, q):
    a=1-d/(2*q); b=1+d/(2*(1-q))
    return q*a**3+(1-q)*b**3-1

def bernstein_nonnegative(c, left, right):
    # Substitute q=left+(right-left)t, then convert power basis to degree-4 Bernstein.
    e=[sum((c[j]*comb(j,k)*left**(j-k)*(right-left)**k
            for j in range(k,5)), Q(0)) for k in range(5)]
    b=[sum((e[k]*Q(comb(i,k),comb(4,k)) for k in range(i+1)), Q(0))
       for i in range(5)]
    return all(z>=0 for z in b)

def verify(row):
    d=Q(row['d']); c=row['certificate']; assert 0<d<2
    if c['branch']=='outer':
        assert d>=Q(2,3) and Q(c['q'])==d/2
        assert Q(c['H_exact'])==d*(4-d)/(2-d)**2
    elif c['branch']=='inner_exact_root':
        q=Q(c['q_exact']); assert 0<d<Q(2,3) and Q(1,3)<q<Q(1,2)
        assert stationary(d,q)==0 and Q(c['H_exact'])==moment(d,q)
    else:
        assert c['branch']=='inner' and 0<d<Q(2,3)
        l,u=Q(c['q_lo']),Q(c['q_hi']);lo,hi=Q(c['H_lo']),Q(c['H_hi'])
        assert Q(1,3)<=l<u<=Q(1,2) and stationary(d,l)>0>stationary(d,u)
        assert 0<lo<=hi
        # 4q²(1-q)²*(G(q)-lo) and *(hi-G(q)) are nonnegative quartics.
        low=[-d**3/2,d*d*(3+d),-3*d*d-4*lo,8*lo,-4*lo]
        high=[d**3/2,-d*d*(3+d),3*d*d+4*hi,-8*hi,4*hi]
        assert bernstein_nonnegative(low,l,u) and bernstein_nonnegative(high,l,u)

if __name__=='__main__':
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('certificates_minimal.json')
    rows=json.loads(path.read_text())['records']
    for row in rows: verify(row)
    print('independent_certificate_replay=%d PASS' % len(rows))
