#!/usr/bin/env python3
"""Exact finite-grid verification of the saved correlated 4x4 candidate.

No MATLAB result is trusted for geometry: this rebuilds every grid rectangle,
overlap mass, and receiving-cell union from the stored rational candidate.
"""
from fractions import Fraction
from decimal import Decimal
from itertools import combinations
from pathlib import Path
import hashlib, json

ROOT=Path('/Users/zhengzhihao/Desktop/T/work')
SOURCE=ROOT/'sm_correlated_grid_probe.json'; OUT=ROOT/'cert_correlated_grid.json'
E=(-16,-4,-1,0,1,2,3,4,5,8,20); DELTA=Fraction(1001,1000); M=256
AMPS=(Fraction(0),Fraction(1),Fraction(3,2),Fraction(2),Fraction(4),Fraction(8),Fraction(16))
def q(x): return Fraction(Decimal(str(x)))
def rects():
    return [(a,b,c,d) for a,b in combinations(E,2) for c,d in combinations(E,2)]
def overlap(r,i,j):
    a,b,c,d=r
    return max(0,min(i+1,b)-max(i,a))*max(0,min(j+1,d)-max(j,c))
def mass(r,g): return sum(Fraction(overlap(r,i,j))*g[i][j] for i in range(4) for j in range(4))
def area(r): return Fraction((r[1]-r[0])*(r[3]-r[2]))
def covered_cells(rs):
    cells=[]
    for i in range(len(E)-1):
        for j in range(len(E)-1):
            if any(r[0]<=E[i] and E[i+1]<=r[1] and r[2]<=E[j] and E[j+1]<=r[3] for r in rs): cells.append((i,j))
    return cells
def cellarea(cell):
    i,j=cell; return Fraction((E[i+1]-E[i])*(E[j+1]-E[j]))
def halo(rs): return sum(map(cellarea,covered_cells(rs)))
def good_rects(g): return [r for r in rects() if mass(r,g)>DELTA*area(r)]
def lnbounds(a):
    assert a>=1
    if a==1:return Fraction(0),Fraction(0)
    u=(a-1)/(a+1); s=sum(u**(2*k+1)/Fraction(2*k+1) for k in range(M+1))
    lo=2*s; rem=2*u**(2*M+3)/(Fraction(2*M+3)*(1-u*u))
    return lo,lo+rem
def costbounds(g):
    lo=hi=Fraction(0)
    for row in g:
        for a in row:
            if a==0: continue
            x,y=lnbounds(a) if a>1 else (Fraction(0),Fraction(0))
            lo+=a*(1+x); hi+=a*(1+y)
    return lo,hi
def dec(x,digits=45):
    from decimal import getcontext, ROUND_FLOOR, ROUND_CEILING, localcontext
    with localcontext() as c:
        c.prec=digits+12; c.rounding=ROUND_FLOOR; lo=Decimal(x.numerator)/Decimal(x.denominator)
    with localcontext() as c:
        c.prec=digits+12; c.rounding=ROUND_CEILING; hi=Decimal(x.numerator)/Decimal(x.denominator)
    return [str(lo),str(hi)]
def rstr(x): return f'{x.numerator}/{x.denominator}'
def baseline():
    best=None; allrows=[]
    for x1 in range(4):
      for x2 in range(x1,4):
       for y1 in range(4):
        for y2 in range(y1,4):
         for a in AMPS[1:]:
            g=[[a if x1<=i<=x2 and y1<=j<=y2 else Fraction(0) for j in range(4)] for i in range(4)]
            h=halo(good_rects(g)); clo,chi=costbounds(g); upper=h/clo
            row={'block':[x1+1,x2+1,y1+1,y2+1],'amplitude':rstr(a),'halo_area':rstr(h),
                 'cost_lower_decimal':dec(clo)[0],'ratio_upper_decimal':dec(upper)[1]}
            allrows.append(row)
            if best is None or upper>best[0]:best=(upper,row)
    return best,allrows
def main():
    raw=SOURCE.read_bytes(); saved=json.loads(raw,parse_float=Decimal)
    b=saved['best_candidate']; g=[[q(a) for a in row] for row in b['amplitude_grid']]
    selected=[tuple(r) for r in b['selected_rectangles']]
    margins=[mass(r,g)-DELTA*area(r) for r in selected]
    assert all(x>0 for x in margins)
    selected_halo=halo(selected); assert selected_halo==Fraction(384)
    full=good_rects(g); full_halo=halo(full)
    clo,chi=costbounds(g); candidate_lower=full_halo/chi
    (bupper,brow), allbase=baseline(); assert candidate_lower>bupper
    exact_full_margins=[mass(r,g)-DELTA*area(r) for r in full]
    result={
      'scope':'Exact certificate only for this finite 4x4 input grid and the listed finite output-edge rectangle family; no claim about C_2 or continuous halos.',
      'source_file':str(SOURCE),'source_sha256':hashlib.sha256(raw).hexdigest(),
      'out_edges':list(E),'threshold':'1001/1000','log_series':{'terms_M':M,'formula':'log A=2 sum_{k=0}^M u^(2k+1)/(2k+1)+R, u=(A-1)/(A+1), 0<=R<=2u^(2M+3)/((2M+3)(1-u^2))'},
      'amplitude_grid_rational':[[rstr(a) for a in row] for row in g],
      'selected_rectangles':[list(r) for r in selected],
      'selected_rectangle_margins_mass_minus_threshold_area':[rstr(x) for x in margins],
      'selected_union_area_exact':rstr(selected_halo),'full_good_rectangle_count':len(full),
      'full_good_union_area_exact':rstr(full_halo),'candidate_cost_decimal_outward':[dec(clo)[0],dec(chi)[1]],
      'candidate_ratio_lower_decimal_outward':dec(candidate_lower),
      'axisblock_constant_family_size':len(allbase),'best_axisblock_constant':brow,
      'best_baseline_ratio_upper_decimal_outward':dec(bupper),
      'strict_ratio_gap_lower_decimal_outward':dec(candidate_lower-bupper),
      'all_selected_strict':True,'comparison_certified':True}
    OUT.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
if __name__=='__main__':main()
