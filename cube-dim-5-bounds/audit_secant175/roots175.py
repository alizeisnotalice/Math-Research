"""Independent rational fourth-root binary search; no isqrt or producer import."""
from math import gcd
from division175 import gain
N=1<<160
B=1<<64

def root_ratio(n,d):
    assert type(n) is int and n>=0 and type(d) is int and d>0
    lo=0;hi=1<<max(0,(n.bit_length()-d.bit_length()+4)//4);checks=0
    assert lo**4*d<=n<hi**4*d
    while hi-lo>1:
        mid=(lo+hi)>>1;checks+=1
        if mid**4*d<=n:lo=mid
        else:hi=mid
    left=n-lo**4*d;right=(lo+1)**4*d-n
    assert left>=0 and right>0
    return lo,left,right,checks

def witness(s):
    L,U=s['range_L_i'],s['range_U_i'];assert 0<=L<U and U>=1<<144
    assert L<=s['mean162_i']<=s['mean112_i']<=U
    Y,yl,yr,ychecks=root_ratio(U<<256,N)
    Z,zl,zr,zchecks=root_ratio(L<<256,U)
    assert Y>0 and 0<=Z<=B
    # Homogeneous Horner sums: different algebra from producer B*D+Z^4.
    D=1
    for j in range(1,4):D=D*Z+B**j
    P=D*Z+B**4
    a,b=Y*P,B*B*D;g=gcd(a,b);a//=g;b//=g
    sn,sd=s['slope164_numerator'],s['slope164_denominator']
    assert sn>0 and sd>0 and 16*sn>=sd
    chosen=a*sd>sn*b
    if chosen:sn,sd=a,b
    g=gcd(sn,sd);sn//=g;sd//=g
    assert sn>0 and sd>0 and 16*sn>=sd
    delta=s['mean112_i']-s['mean162_i']
    q,A,C,r=gain(s['weight'],delta,sn,sd,s['level'])
    # The published remainder uses the uncancelled ratio. No second division.
    rawden=sd<<(4*s['level']);rawrem=s['weight']*delta*sn-q*rawden
    assert 0<=rawrem<rawden
    value=s['self112']-q;assert s['selfL']<=value<=s['self166']
    result={'key':s['key'],'Y':Y,'Z':Z,'Y_left_residual':yl,'Y_next_residual':yr,
        'Z_left_residual':zl,'Z_next_residual':zr,'new_secant_numerator':a,'new_secant_denominator':b,
        'selected_slope_numerator':sn,'selected_slope_denominator':sd,'selected_new_secant':chosen,
        'total_gain_relative112_units':q,'division_remainder':rawrem,'new_selfU':value,
        'baseline166_selfU':s['self166'],'strict_self_gain_units':s['self166']-value}
    proof={'root_comparisons':[ychecks,zchecks],'cancelled_gain_numerator':str(A),
        'cancelled_gain_denominator':str(C),'binary_gain_remainder':str(r)}
    return result,proof
