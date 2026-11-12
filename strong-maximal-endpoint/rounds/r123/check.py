"""Exact time-fiber balance for two fixed transverse box profiles."""
from fractions import Fraction as F
from itertools import product
from concurrent.futures import ProcessPoolExecutor

def lengths(seq,p):
    balance=0
    result=[]
    for g in seq:
        slope=g-p
        result.append(1 if slope>0 or (slope==0 and balance>0) else
                      min(1,balance/(-slope)) if slope<0 else 0)
        balance=max(0,balance+slope)
    return result

def verify(first):
    p,a=F(1,2),F(1,4)
    for second in product([0,p,1],repeat=3):
        # Long zero buffers also account for exterior receptors.
        profiles=[[0]*4+list(row)+[0]*4 for row in (first,second)]
        l=list(map(max,zip(*(lengths(row,p) for row in profiles))))
        r=list(map(max,zip(*(lengths(row[::-1],p)[::-1] for row in profiles))))
        pl=pr=delta=upper=volume=0
        for g,left,right in zip(map(max,zip(*profiles)),l,r):
            o=max(0,left+right-1)
            u=left+right-o
            pl+=left*(g-p)
            pr+=right*(g-p)
            delta+=(1-u)*max(g-a,0)+u*max(a-g,0)+(1-o)*max(g-p,0)+o*max(p-g,0)
            upper+=max(g-a,0)+max(g-p,0)
            volume+=u
        assert min(pl,pr,delta)>=0 and (p-a)*volume==upper-delta-pl-pr
    return 27

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=12) as pool:
        print(sum(pool.map(verify,product([0,F(1,2),1],repeat=3))))
