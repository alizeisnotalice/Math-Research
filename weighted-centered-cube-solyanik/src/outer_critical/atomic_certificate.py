"""Exact rational enclosures for finite atomic weak-type lower bounds.

Halo interval lists are canonical modulo finitely many endpoints. All weighted
power bounds use checked integer inequalities; no floating point accepts a proof.
"""
import json
from pathlib import Path
from fractions import Fraction as Q

ROOT=Path(__file__).resolve().parent
BITS=200


def integer_root(n,b):
    assert n>=0 and b>=1
    if n==0 or b==1:
        return n
    x=1<<((n.bit_length()+b-1)//b)
    while True:
        y=((b-1)*x+n//x**(b-1))//b
        if y>=x:
            break
        x=y
    while x**b>n:x-=1
    while (x+1)**b<=n:x+=1
    assert x**b<=n<(x+1)**b
    return x


def power_bounds(x,e):
    x,e=Q(x),Q(e)
    assert x>=0
    if e.denominator==1:
        assert x>0 or e>=0
        value=x**int(e)
        return value,value
    if x==0:
        assert e>0
        return Q(0),Q(0)
    if e<0:
        x,e=1/x,-e
    a,b=e.numerator,e.denominator
    scale=1<<BITS
    numerator=x.numerator**a*scale**b
    denominator=x.denominator**a
    k=integer_root(numerator//denominator,b)
    assert k**b*denominator<=numerator<(k+1)**b*denominator
    return Q(k,scale),Q(k+1,scale)


def canonical(intervals):
    out=[]
    for a,b in sorted(intervals):
        if a>=b:continue
        if out and a<=out[-1][1]:
            out[-1]=(out[-1][0],max(b,out[-1][1]))
        else:out.append((a,b))
    return out


def block_halo(atoms):
    xs=[x for x,m in atoms];prefix=[Q(0)]
    for x,m in atoms:prefix.append(prefix[-1]+m)
    return canonical((xs[k]-(prefix[k+1]-prefix[j])/2,
                      xs[j]+(prefix[k+1]-prefix[j])/2)
                     for j in range(len(atoms)) for k in range(j,len(atoms)))


def subset_halo(atoms):
    intervals=[]
    for mask in range(1,1<<len(atoms)):
        chosen=[atoms[j] for j in range(len(atoms)) if mask&(1<<j)]
        total=sum((m for x,m in chosen),Q(0))
        intervals.append((max(x for x,m in chosen)-total/2,
                          min(x for x,m in chosen)+total/2))
    return canonical(intervals)


def primitive_bounds(x,tau):
    low,high=power_bounds(abs(x),tau)
    return (low/tau,high/tau) if x>=0 else (-high/tau,-low/tau)


def certificate(record):
    tau=Q(str(record['tau']))
    combined={}
    for x,m in record['atoms']:
        x,m=Q(x),Q(m)
        assert x!=0 and m>0
        combined[x]=combined.get(x,Q(0))+m
    atoms=sorted(combined.items())
    assert 1<=len(atoms)<=8 and all(abs(x)<=1 for x,m in atoms)
    assert max(abs(x) for x,m in atoms)==1
    assert 0<tau<=1
    halo=block_halo(atoms)
    assert halo==subset_halo(atoms)
    if 'expected_halo' in record:
        assert halo==[tuple(map(Q,pair)) for pair in record['expected_halo']]
    mass_low=mass_high=Q(0)
    for a,b in halo:
        al,ah=primitive_bounds(a,tau);bl,bh=primitive_bounds(b,tau)
        mass_low+=bl-ah;mass_high+=bh-al
    cost_low=cost_high=Q(0)
    for x,m in atoms:
        lo,hi=power_bounds(abs(x),tau-1)
        cost_low+=m*lo;cost_high+=m*hi
    assert mass_low>0 and cost_low>0
    lo,hi=mass_low/cost_high,mass_high/cost_low
    assert 0<lo<=hi
    if tau==1:
        z=12*hi-11
        assert z<=0 or z*z<=61  # exact comparison to (11+sqrt(61))/12
    bl,bh=power_bounds(Q(5,3),tau)
    bl*=Q(3,4)/tau;bh*=Q(3,4)/tau
    relation='overlap'
    if hi<bl:relation='strictly_below_pair'
    if lo>bh:relation='strictly_above_pair'
    return dict(id=record['id'],tau=str(tau),lower=str(lo),upper=str(hi),
                halo=[[str(a),str(b)] for a,b in halo],pair_comparison=relation,
                approximate_midpoint=float((lo+hi)/2))


def a1_root_enclosure(tau):
    """Certify the unique root u^(1-tau)+(1-tau)u=tau in (0,1)."""
    tau=Q(tau);assert 0<tau<1
    lo,hi=Q(0),Q(1)
    for _ in range(100):
        mid=(lo+hi)/2
        pl,ph=power_bounds(mid,1-tau)
        gl,gh=pl+(1-tau)*mid-tau,ph+(1-tau)*mid-tau
        if gh<0:lo=mid
        elif gl>0:hi=mid
        elif gl==gh==0:
            lo=hi=mid;break
        else:
            raise ArithmeticError('Root sign unresolved: increase BITS')
    assert 0<lo<=hi<1
    # Negative exponent reverses the root bounds.
    al,_=power_bounds(hi,tau-1)
    _,ah=power_bounds(lo,tau-1)
    bl,bh=power_bounds(Q(5,3),tau);bl*=Q(3,4)/tau;bh*=Q(3,4)/tau
    ml,mh=power_bounds(Q(61),Q(1,2));ml=(11+ml)/12;mh=(11+mh)/12
    return dict(tau=str(tau),u_lower=str(lo),u_upper=str(hi),
                A1_lower=str(al),A1_upper=str(ah),
                operator_W_lower=str(max(bl,ml)),operator_W_upper=str(2*ah),
                approximate_W_upper=float(2*ah))


def main():
    inputs=json.loads((ROOT/'candidates.json').read_text())
    assert len(inputs)==175 and len({r['id'] for r in inputs})==175
    results=[certificate(r) for r in inputs]
    exact1=Q(135)-Q(17,5)**4
    exact2=Q(15,4)-Q(242,125)**2
    assert exact1==Q(854,625)>0 and exact2==Q(119,62500)>0
    summary=dict(status='passed',cases=len(results),bits=BITS,
                 root_method='integer nth-root plus exact rational inequalities',
                 pair_bounds={'W_1/4 > 17/5':str(exact1),'W_1/2 > 242/125':str(exact2)},
                 groups=[],operator_bounds=[a1_root_enclosure(t)
                     for t in (Q(1,4),Q(1,2),Q(3,4),Q(9,10))])
    for tau in sorted({r['tau'] for r in results},key=Q):
        group=[r for r in results if r['tau']==tau]
        best=max(group,key=lambda r:Q(r['lower']))
        summary['groups'].append(dict(tau=tau,cases=len(group),best_certified_lower=best,
            all_sample_values_below_melas=all(
                12*Q(r['upper'])-11<=0 or (12*Q(r['upper'])-11)**2<61
                for r in group),
            above_pair=sum(r['pair_comparison']=='strictly_above_pair' for r in group),
            below_pair=sum(r['pair_comparison']=='strictly_below_pair' for r in group),
            overlap=sum(r['pair_comparison']=='overlap' for r in group)))
    (ROOT/'certificates.json').write_text(json.dumps(results,indent=2))
    (ROOT/'certificate_summary.json').write_text(json.dumps(summary,indent=2))
    print('PASS:',len(results),'halo comparisons and 200-bit rational power enclosures.')
    print('Exact pair inequalities:',exact1,exact2)
    print('Certified A1 roots and operator upper bounds:',len(summary['operator_bounds']))
    for row in summary['groups']:
        print(row['tau'],row['best_certified_lower']['approximate_midpoint'],
              'above_pair',row['above_pair'],'below_pair',row['below_pair'],
              'overlap',row['overlap'])


if __name__=='__main__':main()
