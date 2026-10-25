#!/usr/bin/env python3
"""Summarize the frozen Sobol H-epsilon values without quadrature claims."""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import hashlib
import json
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent


def dist2_s2(row, den):
    x=tuple(Fraction(v,den) for v in row); best=None
    for size in range(1,7):
        for tail in combinations(range(1,6),size-1):
            ones={0,*tail}; rest=[i for i in range(6) if i not in ones]
            a=sum((x[i] for i in rest),Fraction())/len(rest) if rest else Fraction()
            candidate=sum(((x[i]-1)**2 for i in ones),Fraction())+sum(((x[i]-a)**2 for i in rest),Fraction())
            best=candidate if best is None or candidate<best else best
    return best


def stats(rows):
    q=mp.mpf(6)/5; hs=[mp.mpf(r['H_epsilon_binary64']) for r in rows]; powers=[h**q for h in hs]
    total=sum(powers); top=sorted(powers,reverse=True)[:max(1,int(mp.ceil(len(rows)/100)))]
    return {"mean_H_to_q":mp.nstr(total/len(rows),40),"formal_norm_H_q":mp.nstr((total/len(rows))**(1/q),40),
      "maximum_H":mp.nstr(max(hs),40),"top_one_percent_count":len(top),"top_one_percent_share":mp.nstr(sum(top)/total,40)}


def main():
    result=json.loads((HERE/'full512_results.json').read_text()); meta=json.loads((HERE/'input_metadata.json').read_text())
    rows=result['points']; assert len(rows)==512
    positive=all(mp.mpf(r['H_epsilon_binary64'])>0 for r in rows)
    mp.mp.dps=80; q=mp.mpf(6)/5; target=mp.mpf(1)/16000
    by_rep=[]
    for rep in range(4): by_rep.append(stats(rows[128*rep:128*(rep+1)]))
    powers=[mp.mpf(r['H_epsilon_binary64'])**q for r in rows]; index=max(range(len(rows)),key=lambda i:powers[i])
    distances=[dist2_s2(r['numerators'],r['denominator']) for r in rows]
    output={
      "status":"numerical_route_selection_only", "points":len(rows),"replicates":4,"points_per_replicate":128,
      "q":"6/5","all_H_positive":positive,"target_norm_H_q":"1/16000","target_I_H":mp.nstr(target**q,40),
      "aggregate":stats(rows),"replicates_summary":by_rep,
      "max_contribution_index":index,"max_contribution_label":meta['labels'][index],"max_contribution_point":rows[index],
      "max_cancellation_ratio":mp.nstr(max(mp.mpf(r['cancellation_ratio']) for r in rows),30),
      "nearest_S2_distance_min":mp.nstr(mp.sqrt(min(distances)),40),
      "nearest_S2_distance_at_max_contribution":mp.nstr(mp.sqrt(distances[index]),40),
      "warnings":["Fixed scrambled-Sobol samples supply no integration error certificate.","H_epsilon is not G; no tail or global g estimate is included.","C++ values are binary64 exploration pending one Arb point calibration."],
      "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"input_sha256":hashlib.sha256((HERE/'full512_results.json').read_bytes()).hexdigest(),
    }
    (HERE/'probe_summary.json').write_text(json.dumps(output,indent=2)+'\n'); print(json.dumps(output,indent=2))


if __name__=='__main__': main()
