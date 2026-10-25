#!/usr/bin/env python3
"""One 120-pair rank3-quotient/three-linear-factor flag batch."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import sys
import time

import flint
from flint import arb, ctx

HERE = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def ar(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def det_exact(rows):
    a=[list(map(F,row)) for row in rows]
    ans=F(1)
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]),None)
        if j is None:
            return F(0)
        if j!=i:
            a[i],a[j]=a[j],a[i]
            ans=-ans
        p=a[i][i]
        ans*=p
        for j in range(i+1,len(a)):
            f=a[j][i]/p
            for k in range(i+1,len(a)):
                a[j][k]-=f*a[i][k]
            a[j][i]=0
    return ans


def classes(indices,cross):
    unseen=set(indices)
    answer=[]
    while unseen:
        i=min(unseen)
        line={j for j in unseen if cross[i][j]==0}
        assert all(cross[u][v]==0 for u in line for v in line)
        assert all(cross[u][v]>0 for u in line for v in set(indices)-line)
        answer.append(line)
        unseen-=line
    return answer


def histogram(coeff):
    total=sum(coeff)
    assert total>0
    return tuple(sorted(Counter(F(c,total) for c in coeff if c).items()))


def parse(s):
    c,r=map(F,s[1:-1].split(' +/- '))
    return c-r,c+r


def main():
    started=time.monotonic()
    reg=json.loads((HERE/'registration.json').read_text())
    for row in reg['inputs']:
        assert sha(Path(row['source']))==row['sha256']
    oldpath=next(Path(r['source']) for r in reg['inputs'] if r['source'].endswith('upper_I021/results.json'))
    old=json.loads(oldpath.read_text())
    oldJ=parse(old['arb_384bit']['J_upper'])
    oldI=parse(old['arb_384bit']['complete_I_upper'])
    V=[(1,)+s for s in product((-1,1),repeat=4)]
    VV=[[dot(v,w) for w in V] for v in V]
    W,denominator=16**5,60*16**5
    types={}
    pair_types=Counter()
    occupied=[0]*16
    total_C=0
    pair_C_sum=0
    flag_count=0
    example=None
    stream=hashlib.sha256()

    for J in combinations(range(16),2):
        h=VV[J[0]][J[1]]
        D=25-h*h
        assert D in (16,24)
        adj=((5,-h),(-h,5))
        z=[(VV[J[0]][i],VV[J[1]][i]) for i in range(16)]
        za=[(z[i][0]*adj[0][0]+z[i][1]*adj[1][0],
             z[i][0]*adj[0][1]+z[i][1]*adj[1][1]) for i in range(16)]
        A=[[D*VV[i][j]-dot(za[i],z[j]) for j in range(16)] for i in range(16)]
        nonzero={i for i in range(16) if A[i][i]}
        assert set(range(16))-nonzero==set(J)
        cross=[[A[i][i]*A[j][j]-A[i][j]**2 for j in range(16)] for i in range(16)]
        assert all(x>=0 for row in cross for x in row)
        Xclasses=classes(nonzero,cross)
        pair_types[(D,len(Xclasses),tuple(sorted(len(x) for x in Xclasses)))]+=1
        local_C=0
        for X in Xclasses:
            r=min(X)
            xr=A[r][r]
            xsum=F(sum(A[i][i] for i in X),xr)
            PX=histogram([A[i][i] for i in X])
            B=[[xr*A[i][j]-A[i][r]*A[r][j] for j in range(16)] for i in range(16)]
            remaining={i for i in range(16) if B[i][i]}
            assert remaining==nonzero-X
            divisor=D*D*xr
            w=[[0]*16 for _ in range(16)]
            for i,k in combinations(range(16),2):
                numerator=B[i][i]*B[k][k]-B[i][k]**2
                assert numerator>=0 and numerator%divisor==0
                w[i][k]=w[k][i]=numerator//divisor
            assert sum(w[i][k] for i,k in combinations(range(16),2))==xr*16**2
            Yclasses=classes(remaining,w)
            C_for_X=0
            for Y in Yclasses:
                s=min(Y)
                ysum=F(sum(B[i][i] for i in Y),B[s][s])
                outside=remaining-Y
                L=[w[s][k] for k in outside]
                lsum=sum(L)
                assert lsum>0
                assert all(w[i][k]*B[s][s]==B[i][i]*w[s][k] for i in Y for k in outside)
                Cflag=xsum*ysum*lsum
                assert Cflag.denominator==1
                Cflag=int(Cflag)
                support=[set(J),X,Y,outside]
                assert set.union(*support)==set(range(16))
                assert sum(map(len,support))==16
                PY=histogram([B[i][i] for i in Y])
                PL=histogram(L)
                key=tuple(sorted((PX,PY,PL)))
                if key not in types:
                    types[key]={'sum_C':0,'count':0,'histograms':key,
                                'example_pair':J,'example_X':sorted(X),'example_Y':sorted(Y)}
                row=types[key]
                row['sum_C']+=Cflag
                row['count']+=1
                total_C+=Cflag
                local_C+=Cflag
                C_for_X+=Cflag
                flag_count+=1
                for i in range(16):
                    occupied[i]+=Cflag
                stream.update(json.dumps([J,sorted(X),sorted(Y),Cflag,str(key)],separators=(',',':')).encode())
                if example is None:
                    example=(J,X,Y,outside,A,B,w,r,s)
            assert C_for_X==xsum*2*xr*16**2
        assert local_C==6*D*16**3
        pair_C_sum+=D*16**3

    assert pair_C_sum==10*W and total_C==denominator
    assert occupied==[denominator]*16
    assert sum(F(row['sum_C'],denominator) for row in types.values())==1
    # One flag's full separated-coefficient product, checked with independent 5x5 arithmetic.
    J,X,Y,outside,A,B,w,r,s=example
    example_count=0
    for i,j,k in product(sorted(X),sorted(Y),sorted(outside)):
        predicted=F(A[i][i],A[r][r])*F(B[j][j],B[s][s])*w[s][k]
        actual=det_exact([V[a] for a in J]+[V[i],V[j],V[k]])**2
        assert actual==predicted
        example_count+=1

    beta,b,q=F(3,4),F(3,20),F(10,7)
    ctx.prec=384
    gain=arb(1)
    rows=[]
    for key,row in sorted(types.items()):
        theta=F(row['sum_C'],denominator)
        R=arb(1)
        for hist in row['histograms']:
            R*=sum((count*ar(p)**ar(beta) for p,count in hist),arb(0))
        gain*=R**ar(-theta/5)
        rows.append({'count':row['count'],'sum_C':row['sum_C'],'theta':str(theta),
                     'probability_histograms':[{str(p):count for p,count in hist} for hist in key],
                     'example_pair':row['example_pair'],'example_X':row['example_X'],'example_Y':row['example_Y'],
                     'R_product_ball':R.str(85)})
    moment=ar(1+b/beta).gamma()/ar(1+b).gamma()
    Jnew=2*ar(F(7,4)).gamma()*moment**5*gain
    c=3/(ar(2)**ar(F(15,2))*arb.pi()**ar(F(5,2)))
    Inew=c**ar(q)*arb.pi()**2/60*Jnew**ar(q)
    upper=F(str(Inew.upper().fmpq()))
    cap=F((upper*10**13).__ceil__(),10**13)
    assert Inew<ar(cap)
    checks={'J_improves_locked021':bool(Jnew<ar(oldJ[0])),
            'I_improves_locked021':bool(Inew<ar(oldI[0])),
            'I_below_target_11e_minus6':bool(Inew<ar(F(11,10**6))),
            'I_below_118e_minus7':bool(Inew<ar(F(118,10**7)))}
    out={'status':'finite_two_subset_three_linear_factor_batch_complete','new_C5_bound':False,
         'script_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.json'),
         'source_inputs':reg['inputs'],
         'geometry':{'pair_count':120,'flag_count':flag_count,'type_count':len(types),
                     'sum_Cflag':total_C,'expected_60W':denominator,'coordinate_scope_rho':'1',
                     'pair_types':[{'Gram_det':D,'line_class_count':n,'line_sizes':sizes,'count':count}
                                   for (D,n,sizes),count in sorted(pair_types.items())],
                     'probability_types':rows,'geometry_stream_sha256':stream.hexdigest(),
                     'independent_first_flag_checks':example_count},
         'parameters':{'beta':str(beta),'b':str(b),'q':str(q),'rho':'1'},
         'arb_384bit':{'three_linear_factor_gain':gain.str(85),'ES_minus_b':moment.str(85),
                      'J_upper':Jnew.str(85),'complete_I_upper':Inew.str(85),
                      'new_over_locked021_I_lower':(Inew/ar(oldI[0])).str(85)},
         'strict_rational_I_upper':str(cap),'checks':checks,
         'runtime':{'python':sys.executable,'flint':flint.__version__,'precision_bits':ctx.prec,
                    'seconds':time.monotonic()-started,'five_subset_scan':False,'spatial_or_probability_samples':0},
         'scope':'Full-time/full-face I candidate by complete two-subset flag structure; no type-weight optimization, no new g/C5'}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'flags':flag_count,'types':len(types),'rho':'1',
                      'J':Jnew.str(28),'I':Inew.str(28),'rational_I_upper':str(cap),
                      'checks':checks,'seconds':out['runtime']['seconds']},indent=2))


if __name__=='__main__':
    main()
