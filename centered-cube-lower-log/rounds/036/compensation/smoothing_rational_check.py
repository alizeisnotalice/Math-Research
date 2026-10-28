"""Exact small-source check for the smoothed pair-diagonal formula."""
from fractions import Fraction as F
from pathlib import Path
import json

def primitive(a,delta,x):
    x=x%1;k=len(a)
    clip=lambda z:min(F(1),max(F(0),z))
    return sum((clip((x-z-shift+delta/2)/delta)for z in a for shift in [-1,0,1]),F(0))-k-k*x
def derivative(a,delta,x):
    x=x%1
    return sum((F(1,1)/delta for z in a for shift in [-1,0,1]if abs(x-z-shift)<delta/2),F(0))-len(a)
def direct_hh(a,delta,L):
    cuts=sorted({F(0),F(1),*[(z+sgn*delta/2-shift*L/2)%1 for z in a for sgn in [-1,1]for shift in [-1,1]]})
    total=F(0)
    for lo,hi in zip(cuts,cuts[1:]):
        x=(lo+hi)/2;ell=hi-lo
        h=primitive(a,delta,x+L/2)-primitive(a,delta,x-L/2)
        slope=derivative(a,delta,x+L/2)-derivative(a,delta,x-L/2)
        total+=ell*h*h+slope*slope*ell**3/12
    return total
checks=[]
for u in [[F(99,100),F(1,100),F(1,3)],[F(1,10),F(9,10),F(1,4),F(3,4),F(1,2)]]:
    K=len(u);a=[(j+v)/K for j,v in enumerate(u)]
    d=[(x-y)%1 for x in a for y in a if x!=y]
    for v in [8,32]:
        delta=F(1,4*K*v)
        for L in [F(10,9),F(7,6),F(11,9)]:
            t=L%1
            D=K*t/2-K*K*t*t/2+sum((max(F(0),t-z)for z in d),F(0))
            correction=F(0)
            for z in d:
                x=t-z
                if -delta<x<delta:
                    hinge=((x+delta)**3-2*max(F(0),x)**3)/(6*delta*delta)
                    correction+=hinge-max(F(0),x)
            short=sum(((delta-z)**3/(3*delta*delta)for z in d if z<delta),F(0))
            formula=2*(D-K*delta/6+correction-short)
            direct=direct_hh(a,delta,L)
            assert direct==formula
            checks.append({'K':K,'v':v,'L':str(L),'exact_HH':str(formula)})
out={'scope':'12 exact rational smoothed-HH checks; no random simulation.','all_equal':True,'checks':checks}
Path(__file__).with_name('smoothing_rational_results.json').write_text(json.dumps(out,indent=2)+'\n')
print('12 exact rational smoothed pair/direct piecewise-integral checks passed')
