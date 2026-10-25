"""Physical algebraic L-infinity clearance; never use tangent-cone ranks as metric."""
from pathlib import Path
from itertools import product
import json
from round005_pure_satellite_check import AA,ZERO,interval,compare,add,scale,baseline_boxes,rvec
P=Path(__file__).resolve().parent

def center(ns):
    out=[]
    for j,n in enumerate(ns):
        raw=scale(tuple(map(int,rvec(n))),2);ll,hh=interval(raw);al,ah=interval(AA[j]);rat=[x/y for x,y in product([ll,hh],[al,ah])];m=min(rat).numerator//min(rat).denominator;assert m==max(rat).numerator//max(rat).denominator
        z=add(raw,scale(AA[j],-m))
        if compare(scale(z,2),AA[j])>0:z=add(AA[j],scale(z,-1))
        out.append(z)
    return out

def main():
    selected=json.loads((P/'round006_mixed_breakpoint_results.json').read_text())['rows'];selected=[max(selected,key=lambda r:r['prediction'])]+json.loads((P/'round006_screen_results.json').read_text())['runs'];boxes=baseline_boxes();rows=[]
    for row in selected:
        c=center(row['source_masses']);distance=None
        for shifts in product([-1,0,1],repeat=4):
            for lo,hi in boxes:
                d=ZERO
                for j in range(4):
                    disp=scale(AA[j],shifts[j]);l=add(lo[j],disp);h=add(hi[j],disp)
                    g=add(l,scale(c[j],-1)) if compare(c[j],l)<0 else (add(c[j],scale(h,-1)) if compare(c[j],h)>0 else ZERO)
                    if compare(g,d)>0:d=g
                if distance is None or compare(d,distance)<0:distance=d
        ll,hh=interval(distance);assert ll>0
        threshold=(2*ll)**4
        rows.append(dict(source_masses=row['source_masses'],distance_coefficient_vector=list(distance),physical_clearance_lower=str(ll),physical_clearance_upper=str(hh),clearance_float=float((ll+hh)/2),static_singleton_disjoint_epsilon_upper_lower_bound=str(threshold),threshold_float=float(threshold),scope='Static baseline U only; moving-box remainder must still be bounded.'))
    out=dict(status='rigorous baseline physical clearance, not a full moving-family certificate',rows=rows);(P/'round006_clearance_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({**out,'rows':[{k:v for k,v in r.items() if k not in ['physical_clearance_lower','physical_clearance_upper','static_singleton_disjoint_epsilon_upper_lower_bound','distance_coefficient_vector']} for r in rows]},indent=2))
if __name__=='__main__':main()
