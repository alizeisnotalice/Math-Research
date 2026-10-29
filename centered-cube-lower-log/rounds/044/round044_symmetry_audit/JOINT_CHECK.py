#!/usr/bin/env python3
"""Exact Krawtchouk and score-retaining multi-table tests; no new C_n estimate."""
from fractions import Fraction as F
from itertools import product,permutations
from math import comb,prod
from pathlib import Path
import json,time
import numpy as np
import checks as source
HERE=Path(__file__).resolve().parent

def K(n,d,k):
    return sum((-1)**(d-j)*comb(k,j)*comb(n-k,d-j)
        for j in range(max(0,d-(n-k)),min(d,k)+1))
def main():
    start=time.time();rng=np.random.default_rng(440045)
    orth=0;cauchy=0;spectral=0;joint=0;physical=0;largest=F(0)
    for n in range(1,9):
        binlaw=[F(comb(n,k),2**n)for k in range(n+1)]
        for d,e in product(range(n+1),repeat=2):
            v=sum(binlaw[k]*K(n,d,k)*K(n,e,k)for k in range(n+1))
            assert v==(comb(n,d)if d==e else 0);orth+=1
        w=[[F(0)]*(n+1)for b in range(4)]
        for k in range(n+1):
            nums=[int(v)for v in rng.integers(1,10,size=4)]
            for b in range(4):w[b][k]=F(nums[b],sum(nums))
        pi=[sum(binlaw[k]*row[k]for k in range(n+1))for row in w]
        coeff2=[]
        for b in range(4):
            coeff2.append([sum(binlaw[k]*w[b][k]*K(n,d,k)for k in range(n+1))**2/F(comb(n,d))for d in range(n+1)])
        for d in range(n+1):
            assert sum(coeff2[b][d]/pi[b]for b in range(4))<=1;cauchy+=1
        if n>4:continue
        for q in (1,2,3):
            cap=F(1,(2*q+1)**2);v=F(1,10*(2*q+1))
            scores=list(permutations([F(1,2),F(1),F(3,2)]))[:4]
            pp=[];laws=[];lam=[]
            for a in scores:
                p=[F(1,2)+v/a[0],F(1,2)-v/a[1],F(1,2)]
                law=[aa/3 for aa in a]
                assert sum(law)==1 and sum(z*p0 for z,p0 in zip(law,p))==F(1,2)
                ll=4*sum(z*(p0-F(1,2))**2 for z,p0 in zip(law,p));assert ll<=cap
                pp.append(p);laws.append(law);lam.append(ll)
            ET=[F(0)]*4;whole=F(0);separated=F(0)
            for inds in product(range(3),repeat=n):
                A=[];T=[]
                for b in range(4):
                    rows=[(pp[b][i],1-pp[b][i])for i in inds]
                    distribution=source.poly(rows)
                    tb=sum(distribution[k]*w[b][k]for k in range(n+1))
                    ab=prod(scores[b][i]for i in inds)
                    A.append(ab);T.append(tb)
                    ET[b]+=prod(laws[b][i]for i in inds)*tb*tb
                FV=sum(a*t for a,t in zip(A,T));GV=sum(a*p for a,p in zip(A,pi))
                sep=sum(A[b]*T[b]**2/pi[b]for b in range(4))
                assert FV*FV/GV<=sep
                whole+=FV*FV/GV/F(3**n);separated+=sep/F(3**n)
            for b in range(4):
                assert ET[b]==sum(lam[b]**d*coeff2[b][d]for d in range(n+1));spectral+=1
            assert separated==sum(ET[b]/pi[b]for b in range(4))
            assert whole<=separated<=sum(cap**d for d in range(n+1))
            assert whole-1<=F(1,4*q*(q+1));joint+=1;largest=max(largest,whole)
        # Real physical sources, independent exact pointwise Cauchy checks.
        for trial in range(7):
            x=[F(2*i+trial+1,19)for i in range(n)];L=F(5,4)
            A=[];T=[]
            for b in range(4):
                rows=[source.class_window(b,xx,L)for xx in x]
                ab=prod(sum(row)/L for row in rows)
                comp=sum(vv*w[b][k]for k,vv in enumerate(source.poly(rows)))/L**n
                A.append(ab);T.append(comp/ab)
            FV=sum(a*t for a,t in zip(A,T));GV=sum(a*p for a,p in zip(A,pi))
            assert FV*FV/GV<=sum(A[b]*T[b]**2/pi[b]for b in range(4));physical+=1
    out=dict(status='PASS',exact_Krawtchouk_orthogonality_checks=orth,
        exact_cross_table_weight_Cauchy_checks=cauchy,exact_spectral_expansion_checks=spectral,
        exact_joint_abstract_channel_integrals=joint,exact_real_source_pointwise_CS_checks=physical,
        largest_exact_abstract_joint_value=str(largest),
        scope='Exact small finite-channel identities and real-source pointwise inequalities; general integral bound is proved separately; no C_n estimate',seconds=time.time()-start)
    (HERE/'JOINT_CHECK.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
