"""Atomic endpoint batches for finite interval-union one-sided slopes.

Only Fraction/integer inputs. An interval record is (lower,upper,lower_velocity,
upper_velocity), with endpoint derivatives in increasing physical path time.
All endpoints at an event coordinate are processed together, never sequentially.
"""
from fractions import Fraction as F
from normal_intervals import signed_length


def q(x):
    if isinstance(x,float):raise TypeError('floating event positions/velocities are not certificates')
    return F(x)


def union_slope(intervals,side='right'):
    """Return d|union|/dt from right or left, assuming actual endpoint expansions.

    side=left divides by negative time increment; it is not the coefficient
    with respect to a positive backwards-time parameter. Empty positive-gap
    intervals are stable locally; zero-width intervals may be born on one side.
    """
    if side not in ('left','right'):raise ValueError('side must be left/right')
    rows=[tuple(map(q,row)) for row in intervals]
    if any(len(row)!=4 for row in rows):raise ValueError('four values per interval required')
    direction=1 if side=='right' else -1
    positions=sorted(set(x for a,b,va,vb in rows for x in (a,b)))
    coefficient=F(0)
    for c in positions:
        crossing=any(a<c<b for a,b,va,vb in rows)
        left=[direction*vb for a,b,va,vb in rows if a<b and b==c]
        right=[direction*va for a,b,va,vb in rows if a<b and a==c]
        born=[(direction*va,direction*vb) for a,b,va,vb in rows if a==b==c]
        coefficient+=signed_length(left,right,born,crossing=crossing)
    return direction*coefficient


def union_length(intervals):
    """Independent exact endpoint-cell integration for finite interval fixtures."""
    rows=[tuple(map(q,row)) for row in intervals]
    cuts=sorted(set(x for row in rows for x in row))
    return sum((b-a for a,b in zip(cuts,cuts[1:])
                if any(lo<(a+b)/2<hi for lo,hi in rows)),F(0))


def _encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:_encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [_encode(v) for v in x]
    return x


def _check():
    from pathlib import Path
    import hashlib,json
    import sympy as sy
    p=Path(__file__).resolve().parent;path=p/'vg11b_input.json'
    assert hashlib.sha256(path.read_bytes()).hexdigest()=='7644648e0bc6d88a85c0bd684ec9039e451ce690eda39dad166707e152d7c947'
    data=json.loads(path.read_text());event=json.loads((p/'vg11b_event.json').read_text())['theta_certificate']
    caps={k:tuple(map(F,v)) for k,v in data['rational_caps'].items()}
    pl,ph=caps['P1'];Rl,Rh=caps['R36'];r0l,r0h=caps['rL_at_0'];r1l,r1h=caps['rL_at_1']
    assert (2*(pl-1))**2<2<(2*(ph-1))**2
    assert Rl**4<F(9,4)<Rh**4
    assert r0l**4<F(15,8)<r0h**4 and r1l**4<F(39,20)<r1h**4
    F0_upper=r0h-Rl;F1_lower=pl/25+r1l-Rh
    overlap_lower=Rl+r0l-ph*F(26,25)
    Hwidth_lower=2*Rl-ph*F(26,25);Lwidth_lower=2*r0l-ph
    assert F0_upper<0<F1_lower and min(overlap_lower,Hwidth_lower,Lwidth_lower,Rl-r1h)>0
    theta_lo,theta_hi=map(F,event['theta_interval'])
    assert 0<theta_lo<theta_hi<1 and theta_hi-theta_lo==F(1,2**40)
    assert F(event['F_at_lower'][1])<0<F(event['F_at_upper'][0])
    assert F(event['monotone_derivative_lower'])>0
    before_slope=(-ph/25,-pl/25)
    after_slope=(F(6,5)/(64*r1h**3),F(6,5)/(64*r0l**3))
    assert before_slope[1]<0<after_slope[0]
    # One rational event model in the actual path's certified sign/order class.
    a=F(683,20000);b=F(23,2000)
    rows=[(0,2,a,-a),(0,1,-a-b,-a+b)]
    right=union_slope(rows,'right');left=union_slope(rows,'left')
    assert right==b and left==-2*a
    for dt in (F(-1,100),F(1,100)):
        actual=[(q(lo)+q(va)*dt,q(hi)+q(vb)*dt) for lo,hi,va,vb in rows]
        quotient=(union_length(actual)-2)/dt
        assert quotient==(right if dt>0 else left)
    shielded=rows+[(-1,3,0,0)]
    assert union_slope(shielded,'right')==union_slope(shielded,'left')==0
    # No threshold solve: exact symbolic path formulas, supplied event label.
    t=sy.symbols('t',real=True);P1=1+sy.sqrt(2)/2;R=sy.sqrt(6)/2
    r=(30+sy.Rational(6,5)*t)**sy.Rational(1,4)/2
    Fexpr=P1*t/25+r-R;VH=2*R-P1-P1*t/25
    assert sy.simplify(sy.diff(r,t)-sy.Rational(6,5)/(64*r**3))==0
    assert sy.simplify(sy.diff(VH+Fexpr,t)-sy.diff(r,t))==0
    # Exact algebraic identity proving single-pair loss when F>0.
    x,A,B,d=sy.symbols('F A B d',real=True)
    assert sy.expand((-2*d*A+x*B)-(-(2*d-x)*A-x*(A-B)))==0
    # Actual tangential widths are positive, using coarse old-root caps only.
    U=F(6,5);V=F(4,3)
    assert U**4>2 and V**4>3
    P2_upper=(2+U**3*V)/4;P34_upper=(1+U)/2
    spans_upper=[2*P2_upper,P34_upper,2*P34_upper]
    assert all(2*r0l>s for s in spans_upper)
    out=dict(status='one fixed path-B local interval interface PASS; no full old-union derivative',
        event_batches=[dict(time='0',kind='right raw ties; permanent identical candidates grouped'),
                       dict(theta_interval=[theta_lo,theta_hi],kind='H.lower=L.lower processed simultaneously')],
        event_batch_count=2,event_limit=32,
        order_before=['H.lower','L.lower','L.upper','H.upper'],
        order_after=['L.lower','H.lower','L.upper','H.upper'],
        certified_fixed_order_subintervals=[[F(0),theta_lo],[theta_hi,F(1)]],
        guard_band=[theta_lo,theta_hi],guard_band_note='one unique event, no sequential endpoint update inside band',
        F0_upper=F0_upper,F1_lower=F1_lower,overlap_lower=overlap_lower,
        Hwidth_lower=Hwidth_lower,Lwidth_lower=Lwidth_lower,
        left_event_union_slope_interval=before_slope,right_event_union_slope_interval=after_slope,
        baseline_increment_derivative=0,finite_endpoint_increment_lower=F1_lower,
        rational_event_fixture=dict(a=a,b=b,left_slope=left,right_slope=right,shielded_slope=0),
        pair_obstruction='For all t>0, raw single H(t) union L(t) volume < raw H(0) volume; not the full old union',
        n4_threshold_scope='Imported B event certificate; no root isolation or path sampling repeated',
        source_sha256={name:hashlib.sha256((p/name).read_bytes()).hexdigest() for name in ['vg11b_interface.py','vg11b_input.json','vg11b_event.json','normal_intervals.py']})
    (p/'vg11b_results.json').write_text(json.dumps(_encode(out),indent=2)+'\n')
    print(json.dumps(_encode(out),indent=2))


if __name__=='__main__':_check()
