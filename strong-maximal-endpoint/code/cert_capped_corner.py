"""Exact actual-halo witness for a 4x4 two-level input; no solver needed."""
from fractions import Fraction as Q
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
import json, math, sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
F=[[Q(16) for _ in range(4)] for _ in range(4)];F[0][0]=Q(4)
TAU=Q(10001,10000);M=64
assert sum(map(sum,F),Q(0))/16>TAU

def canonical(a,b,c,d):
    a,b,c,d=map(Q,(a,b,c,d))
    if c==0:
        assert a==0 and d!=0
        return (Q(0),b/d,Q(0),Q(1))
    if b*c-a*d==0:return(Q(0),a/c,Q(0),Q(1))
    return(a/c,b/c,Q(1),d/c)
def value(g,x):
    a,b,c,d=g;assert c*x+d!=0
    return(a*x+b)/(c*x+d)
def lower_integral(gs,lo,hi):
    gs=list(dict.fromkeys(gs))
    total=Q(0);cells=(hi-lo)*M;assert cells.denominator==1
    for k in range(int(cells)):
        left=lo+Q(k,M);right=left+Q(1,M);mid=(left+right)/2
        idx=max(range(len(gs)),key=lambda j:value(gs[j],mid))
        g=gs[idx];assert (g[2]*left+g[3])*(g[2]*right+g[3])>0
        height=max(Q(0),min(value(g,left),value(g,right)))
        total+=height/M
    return total

def strip(f):
    total=Q(0)
    for j in range(4):
        gs=[canonical(0,0,0,1)]
        for a in range(4):
            v=[sum(f[i][jj] for i in range(a,4)) for jj in range(4)]
            p=[Q(0)];
            for vv in v:p.append(p[-1]+vv)
            vc=v[j]
            for l in range(j+1):
                for r in range(j+1,5):gs.append(canonical(0,a-4+(p[r]-p[l])/((r-l)*TAU),0,1))
                gs.append(canonical(vc/TAU+a-4,(p[j]-vc*j-p[l])/TAU-(a-4)*l,1,-l))
            for r in range(j+1,5):gs.append(canonical(-vc/TAU-a+4,(p[r]-p[j]+vc*j)/TAU+(a-4)*r,-1,r))
        total+=lower_integral(gs,Q(j),Q(j+1))
    return total

def corner(f):
    gs=[canonical(0,0,0,1)];hi=Q(4)
    for a in range(4):
        for b in range(4):
            mass=sum(f[i][j] for i in range(a,4) for j in range(b,4))
            gs.append(canonical(b-4,mass/TAU-(b-4)*a,1,-a))
            hi=max(hi,a+mass/(TAU*(4-b)))
    return lower_integral(gs,Q(4),Q(math.ceil(hi)))

def transpose(f):return [list(x) for x in zip(*f)]
def flipx(f):return f[::-1]
def flipy(f):return [r[::-1] for r in f]

def log_bounds(a,m=64):
    u=(a-1)/(a+1);low=2*sum((u**(2*k+1)/Q(2*k+1) for k in range(m)),Q(0))
    rem=2*u**(2*m+1)/(Q(2*m+1)*(1-u*u));return low,low+rem

def decimal(q,upper=False):
    with localcontext() as c:
        c.prec=90;c.rounding=ROUND_CEILING if upper else ROUND_FLOOR
        d=Decimal(q.numerator)/Decimal(q.denominator)
        return str(d.quantize(Decimal('1e-30'),rounding=c.rounding))

strips=[strip(f) for f in [F,flipx(F),transpose(F),flipx(transpose(F))]]
corners=[corner(f) for f in [F,flipx(F),flipy(F),flipx(flipy(F))]]
area=Q(16)+sum(strips)+sum(corners)
lo,up=log_bounds(Q(2));cost_up=244+968*up
candidate=area/cost_up
baseline=(Q(1,16)+16*up)/(1+4*lo)
assert candidate>baseline
eps=Q(1,100);amp=1+eps
la,ua=log_bounds(amp);lb,ub=log_bounds(Q(4,3))
near_volume=16+4*eps*(1-eps*eps)*(15+3*lb)
near_cost=1+15*amp*(1+ua)
near_ratio=near_volume/near_cost
near_base=(1/amp+4*ua)/(1+la)
assert near_ratio>near_base
near={'epsilon':str(eps),'uniform_average_margin':str(3*eps**3/(4+4*eps)),
 'witness_volume_lower':decimal(near_volume),'cost_upper':decimal(near_cost,True),
 'ratio_lower':decimal(near_ratio),'constant_rectangle_upper':decimal(near_base,True),
 'strict_gap_lower':decimal(near_ratio-near_base),
 'first_variation_coefficient_lower':decimal((6*lb-1)/8),
 'first_variation_coefficient_upper':decimal((6*ub-1)/8,True)}
result={'near_unit_case':near,'input_amplitudes':[[int(x) for x in r] for r in F], 'average_threshold':str(TAU),
 'subdivisions_per_unit':M,'log2_series_terms':64,'interior_area':16,'interior_average':str(sum(map(sum,F),Q(0))/16),
 'strip_area_lower':[decimal(x) for x in strips],'corner_area_lower':[decimal(x) for x in corners],
 'area_lower_numerator':str(area.numerator),'area_lower_denominator':str(area.denominator),
 'area_lower_decimal':decimal(area),'actual_Orlicz_cost_upper':decimal(cost_up,True),
 'candidate_actual_halo_ratio_lower':decimal(candidate),
 'all_constant_rectangle_amplitudes_at_most_16_ratio_upper':decimal(baseline,True),
 'strict_gap_lower':decimal(candidate-baseline),
 'status':'exact rational witness volume and explicit log-series bounds; no global C2 optimality claim'}
(Path(__file__).resolve().parents[1]/'rounds'/'r023'/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if not k.endswith(('numerator','denominator'))},indent=2))
