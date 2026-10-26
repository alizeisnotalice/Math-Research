#!/usr/bin/env python3
"""Independent fixed H_1/128 first-jet kernel for mixed G leaf bounds; no driver."""
from __future__ import annotations
import argparse
from collections import deque
from fractions import Fraction as F
from itertools import product
from math import comb,lcm
from pathlib import Path
import gzip,hashlib,json,os,resource,sys,time
from flint import arb,fmpq,ctx
import flint

if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ab(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))
def endpoint(x,upper=False):return F(str((x.upper() if upper else x.lower()).fmpq()))
def pack(x):
    if not x.is_finite():return str(x)
    return [str(endpoint(x)),str(endpoint(x,True))]
def unpack(x):
    lo,hi=map(F,x);return arb(fmpq((lo+hi).numerator,2*(lo+hi).denominator),fmpq((hi-lo).numerator,2*(hi-lo).denominator))
def atomic(p,x):
    q=p.with_suffix(p.suffix+'.tmp')
    with q.open('wb') as f:f.write(canonical(x)+b'\n');f.flush();os.fsync(f.fileno())
    os.replace(q,p)
def pointkey(x):return ','.join(map(str,x))
def numeric_quality(values):
    if not all(v.is_finite() for v in values):return False
    for v in values:
        lo,hi=map(F,pack(v))
        if hi-lo>max(abs(lo),abs(hi),F(1,10**6))/2**40:return False
    return values[0]>0


class FirstJetKernel:
    """H_epsilon and four derivatives. No frozen producer is imported.

    Integer projection eligibility is exact. All subsequent real arithmetic is
    outward Arb; partial sums are saved using exact dyadic interval endpoints.
    """
    def __init__(self,table,epsilon=F(1,128)):
        self.vertices=table['vertices'];self.simplices=table['simplices'];self.epsilon=F(epsilon)
        assert self.epsilon==F(1,128) and table['denominator']==144 and len(self.vertices)==280 and len(self.simplices)==3504

    def projections(self,point):
        point=tuple(map(F,point));assert len(point)==5
        D=lcm(*(p.denominator for p in point));n=[int(D*p) for p in point]
        k=[sum(a*b for a,b in zip(n,z)) for z in self.vertices]
        return D,k

    def eligible(self,point):
        D,k=self.projections(point)
        return {'eligible':0 not in k and len(set(k))==280,'zero_count':k.count(0),
                'distinct_count':len(set(k)),'common_denominator':str(D)}

    def evaluate(self,point,bits=256,partial=None,deadline=None,stride=128):
        """Return {complete, values/partial, ...}; a deadline never omits a cone."""
        began=time.monotonic();ctx.prec=bits;point=tuple(map(F,point));D,k=self.projections(point)
        assert 0 not in k and len(set(k))==280,'generic-center-only DD API'
        scaled=F(144*D,128)
        assert scaled==144*D*self.epsilon
        A,B=scaled.numerator,scaled.denominator
        slopes=[tuple(D*z for z in v[1:]) for v in self.vertices]
        jets=[]
        for node,l in zip(k,slopes):
            w=B*B*node*node+A*A
            jets.append([arb(w).log()/2,*[arb(s*B*B*node)/w for s in l]])
        total=[arb(0) for _ in range(5)];absolute=arb(0);start=0
        if partial:
            assert partial['point']==list(map(str,point)) and partial['bits']==bits
            start=partial['next_cone'];total=list(map(unpack,partial['total']));absolute=unpack(partial['absolute'])
        for ci in range(start,3504):
            s=self.simplices[ci];ids=sorted(s['vertices'],key=lambda i:k[i]);table=[]
            for right,idx in enumerate(ids):
                row=[jets[idx]]
                for order in range(1,right+1):
                    left=ids[right-order];d=k[idx]-k[left]
                    delta=[a-b for a,b in zip(slopes[idx],slopes[left])]
                    P=row[-1];Q=table[right-1][order-1]
                    r0=(P[0]-Q[0])/d
                    row.append([r0,*[(P[j+1]-Q[j+1]-delta[j]*r0)/d for j in range(4)]])
                table.append(row)
            jet=table[-1][-1];det=s['det_numerator']
            for j in range(5):total[j]-=det*jet[j]
            absolute+=abs(det*jet[0]);done=ci+1
            if deadline is not None and done%stride==0 and time.monotonic()>=deadline:
                return {'complete':False,'partial':{'point':list(map(str,point)),'bits':bits,'next_cone':done,
                  'total':list(map(pack,total)),'absolute':pack(absolute)},'seconds':time.monotonic()-began}
        factor=arb(D)**4/(144*(2*arb.pi())**5);values=[factor*v for v in total]
        return {'complete':True,'values':values,'seconds':time.monotonic()-began,
          'cancellation_ratio':str(absolute*factor/abs(values[0])) if not values[0].contains(0) else None,
          'minimum_integer_projection_gap':str(min(b-a for a,b in zip(sorted(k),sorted(k)[1:]))),
          'common_denominator':str(D),'all_cones':3504,'bits':bits}

