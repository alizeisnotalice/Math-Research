#!/usr/bin/env python3
"""Round-003 floating discovery screen; it does not certify a lower bound."""
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import argparse
import hashlib, json, time

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / 'results' / 'bounds_ledger.json'
TOL = 1e-13

def boxes_for(a):
    raw=[]
    for ell in product(range(1,5), repeat=6):
        r=prod(ell)**(1/6)/2.0; b=[]
        for aj,k in zip(a,ell):
            if k%2: lo,hi=0.0,min(aj/2,r-(k-1)*aj/2)
            else: lo,hi=max(0.0,k*aj/2-r),aj/2
            if hi <= lo+TOL: break
            b.extend((lo,hi))
        else: raw.append((ell,tuple(b)))
    unique={b:ell for ell,b in raw}
    keep=[]
    for b,ell in unique.items():
        if any(c != b and all(c[2*j] <= b[2*j]+TOL and b[2*j+1] <= c[2*j+1]+TOL for j in range(6)) for c in unique):
            continue
        keep.append((ell,b))
    return raw,keep

def volume(boxes):
    @lru_cache(None)
    def rec(B):
        if not B:return 0.0
        if len(B[0])==2:
            order=sorted(B); a,b=order[0]; total=0.0
            for c,d in order[1:]:
                if c > b+TOL: total,a,b=total+b-a,c,d
                elif d>b: b=d
            return total+b-a
        cuts=sorted({z for b in B for z in b[:2]}); ans=0.0
        for lo,hi in zip(cuts,cuts[1:]):
            if hi <= lo+TOL: continue
            active=tuple(sorted({b[2:] for b in B if b[0] <= lo+TOL and hi <= b[1]+TOL}))
            ans += (hi-lo)*rec(active)
        return ans
    return rec(tuple(sorted(boxes)))

def candidate_list():
    h=1/500
    for p0,p1 in product((-h,h), repeat=2):
        for tail in ((h,h,h,h),(-h,-h,-h,-h),(h,h,-h,-h),(-h,-h,h,h)):
            yield (p0,p1,*tail)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--ledger',type=Path,default=LEDGER)
    args=parser.parse_args()
    started=time.monotonic()
    u=2**(1/6); v=3**(1/6)
    base=(1+u**4/2,(2+u**5*v)/4,*((1+u)/2 for _ in range(4)))
    ledger=json.loads(args.ledger.read_text())
    frozen=ledger['lower']['exact_interval'][1]
    cap=20*60; out=[]
    for i,delta in enumerate(candidate_list()):
        if time.monotonic()-started > cap: break
        spacing=tuple(x+d for x,d in zip(base,delta))
        raw,keep=boxes_for(spacing)
        ratio=64*volume(tuple(b for _,b in keep))
        out.append({'index':i,'delta_rational':[f'{round(x*500):d}/500' for x in delta],
                    'raw_boxes':len(raw),'kept_boxes':len(keep),'ratio_float':ratio,
                    'margin_over_display_upper':ratio-2.06513594276569})
    winner=max(out,key=lambda x:x['ratio_float']) if out else None
    result={'status':'completed' if len(out)==16 else 'wall_time_stop','scope':'numeric discovery only; no C6 claim',
            'base_formula':['1+2^(4/6)/2','(2+2^(5/6)*3^(1/6))/4','(1+2^(1/6))/2 x4'],
            'frozen_comparator_full_exact_upper':frozen,'candidates':out,'winner':winner,
            'promotion_threshold_margin':1e-7,'promotion_required':bool(winner and winner['margin_over_display_upper']>1e-7),
            'elapsed_seconds':time.monotonic()-started,'ledger_sha256':hashlib.sha256(args.ledger.read_bytes()).hexdigest(),
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('screen_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
