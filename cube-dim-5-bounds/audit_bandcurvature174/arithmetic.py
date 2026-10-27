"""174 prepared rational interface. No source load or real evaluation on import.

Caller must separately qualify EPS64 original jets,donor intervals,old155 q,
model/point/time identity and all source hashes before invoking this module.
"""
from fractions import Fraction as F
from math import comb
def exact(x):
    if isinstance(x,(float,bool)):raise TypeError('Exact rational required')
    return F(x)
def difference_weight(time_left,time_right):
    a=exact(time_left);b=None if time_right is None else exact(time_right)
    assert a>0 and (b is None or b>a)
    A,B=1/a,F(0) if b is None else 1/b
    return F(2,3)*(A**6-B**6)/(A**4-B**4)
def positive_weight(time_left,time_right):
    a=exact(time_left);b=None if time_right is None else exact(time_right)
    assert a>0 and (b is None or b>a)
    A,B=1/a,F(0) if b is None else 1/b;d=A-B
    P=lambda n:sum((F(comb(n,k),k+1)*B**(n-k)*d**k for k in range(n+1)),F(0))
    return P(5)/P(3)
def future_q_from_qualified_intervals(original_H64,donor73,time_keys,old155_q):
    """Future authorized case only;no gradient/Delta/log/corner computation.

    donor73 is the exact original list [(lower,upper),...],one per physical t.
    All source identities and raw SHA are a caller binding gate,not this API.
    """
    assert len(time_keys)==len(donor73)==73
    times=list(map(exact,time_keys));assert all(0<a<b for a,b in zip(times,times[1:]))
    assert times[0]==F(1,512) and times[-1]==1
    values=[tuple(map(exact,r)) for r in donor73];assert all(0<l<=u for l,u in values)
    gl,gu=map(exact,original_H64);assert 0<gl<=gu
    start=times.index(F(1,64));assert start==24 and len(times[start:])==49
    dl,du=values[start];overlap=(max(gl,dl),min(gu,du))
    if overlap[0]>overlap[1]:raise ValueError('Same-point H64 intervals do not overlap;stop,no clipping')
    L=max(gl,dl);assert L>0
    K=F(0);bands=[]
    for j in range(start,73):
        nxt=values[j+1] if j<72 else (F(0),F(0));massU=values[j][1]-nxt[0]
        assert massU>0
        right=times[j+1] if j<72 else None
        w=difference_weight(times[j],right);assert w==positive_weight(times[j],right) and w>0
        K+=massU*w;bands.append((massU,w))
    band_cap=K/(16*L);old=exact(old155_q);assert old>0
    # Genuine donor H32 lower and ORIGINAL g-source H64 upper,not L64.
    L32=values[times.index(F(1,32))][0];assert 0<=L32<gu
    A=F(1,16)/F(1,64)**2;split=A*(F(7,10)-F(8,15)*L32/gu)
    new=min(old,band_cap,split);assert 0<new<=old
    return {'K_upper':K,'H64_denominator_lower':L,'same_cutoff_overlap':overlap,
      'band_curvature_upper':band_cap,'split153_upper':split,'old155_q':old,'new_q':new,
      'strict_q_decrease':new<old,'band_count':len(bands),'H32_lower':L32,'original_H64_upper_for_split':gu}
