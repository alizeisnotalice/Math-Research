#!/usr/bin/env python3
"""Directed cone-volume and score upper bound for the frozen r054 family."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
CORE=HERE/'cert_d1_core54.py'

def sha(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()

def load_core():
    spec=importlib.util.spec_from_file_location('r054_core',CORE)
    if spec is None or spec.loader is None: raise SystemExit('cannot load D1 core')
    module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module; spec.loader.exec_module(module)
    return module

def f_ratio(core, z):
    """F(z)=(1-exp(-z))/z, with the removable value F(0)=1."""
    if z==core.Z: return core.pt(core.O)
    numerator=core.sub(core.pt(core.O),core.expiv(core.pt(z.copy_negate())))
    return core.div(numerator,core.pt(z))

def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--d1-json',type=Path,default=HERE.parent/'rounds/r054/d1_54.json')
    p.add_argument('--out',type=Path,default=Path('/tmp/r054_cone_score.json'))
    args=p.parse_args(); repo=args.project_root.resolve()
    fair=repo/'rounds/r052/fair52.json'; seed=repo/'rounds/r030/receiver_seed.json'
    paths={'core':CORE,'fair':fair,'seed':seed,'d1':args.d1_json}
    start={key:sha(path) for key,path in paths.items()}
    d1=json.loads(args.d1_json.read_text(),parse_float=lambda s: s)
    if d1.get('source_sha256')!=start['core'] or d1.get('fair_sha256')!=start['fair'] or d1.get('seed_sha256')!=start['seed']:
        raise SystemExit('D1 record is not bound to these frozen inputs')
    core=load_core(); core.REPO=repo; core.FAIR=fair; core.SEED=seed
    ux,uy,weights,cvec,tau=core.geometry()
    d_lower=core.Decimal(d1['D1_lower'])
    if d_lower<=core.Z: raise SystemExit('D1 lower endpoint must be positive')
    # The finite Pareto-cone formula groups equal x anchors.  For the ordered
    # distinct abscissae x_i and prefix minima b_i=min_{x<=x_i}y,
    # V=sum_i (e^-x_i-e^-x_{i+1})e^-b_i, with e^-x_{m+1}=0.
    points=sorted(zip(ux,uy),key=lambda pair: pair[0])
    grouped=[]; pos=0; prefix=None
    while pos<len(points):
        x=points[pos][0]; same=[]
        while pos<len(points) and points[pos][0]==x:
            same.append(points[pos][1]); pos+=1
        prefix=min(same) if prefix is None else min(prefix,min(same))
        grouped.append((x,prefix))
    volume=core.pt(core.Z)
    for i,(x,best_y) in enumerate(grouped):
        ex=core.expiv(core.pt(x.copy_negate())); ey=core.expiv(core.pt(best_y.copy_negate()))
        if i+1<len(grouped):
            exnext=core.expiv(core.pt(grouped[i+1][0].copy_negate()))
            strip=core.sub(ex,exnext)
        else:
            strip=ex
        if strip.lo<core.Z: raise ArithmeticError(('negative cone strip',i,strip))
        volume=core.add(volume,core.mul(strip,ey))
    # Use V<=V_hi and D1>=D_lower.  Hence tau+log(V)-log(D1) is bounded
    # above by tau+log(V_hi)-log(D_lower), enclosed outwards below.
    logv=core.lniv(core.pt(volume.hi)); logd=core.lniv(core.pt(d_lower))
    score=core.add(core.pt(tau),core.sub(logv,logd))
    # A stronger witness fills the old, strictly Pareto-ordered 65 anchors by
    # straight line segments.  This is a geometric interpolation only: it
    # does not introduce or certify new receiver constraints.
    strict_chain=all(points[i][0]<points[i+1][0] and points[i][1]>points[i+1][1]
                     for i in range(len(points)-1))
    chain_volume=None; chain_score=None; chain_z_range=None
    if strict_chain:
        chain_volume=core.pt(core.Z); zs=[]
        for i in range(len(points)-1):
            x0,y0=points[i]; x1,y1=points[i+1]
            dx=core.exact(lambda x0=x0,x1=x1:x1-x0)
            dy=core.exact(lambda y0=y0,y1=y1:y1-y0)
            z=core.exact(lambda dx=dx,dy=dy:dx+dy); zs.append(z)
            pref=core.expiv(core.pt(core.exact(lambda x0=x0,y0=y0:x0+y0).copy_negate()))
            chain_volume=core.add(chain_volume,core.mul(core.pt(dx),core.mul(pref,f_ratio(core,z))))
        xlast,ylast=points[-1]
        chain_volume=core.add(chain_volume,core.expiv(core.pt(core.exact(lambda:xlast+ylast).copy_negate())))
        chain_score=core.add(core.pt(tau),core.sub(core.lniv(core.pt(chain_volume.hi)),logd))
        chain_z_range=[str(min(zs)),str(max(zs))]
    end={key:sha(path) for key,path in paths.items()}
    if start!=end: raise SystemExit('core/input/D1 record changed during calculation')
    anchors_text=';'.join(f'{x},{y}' for x,y in points)
    record={
      'scope':'Fixed r052 kappa=-0.08 65-receiver finite Pareto-cone score only; no assertion about a full halo, arbitrary receiver subsets, or C_N.',
      'formula':'V=sum_i (exp(-x_i)-exp(-x_{i+1})) exp(-b_i), b_i=min_{j:x_j<=x_i} y_j, exp(-x_{m+1})=0.',
      'arithmetic':'Imports only the frozen core interval primitives; exp and log endpoints are outwardly enclosed.',
      'anchor_count':len(points),'distinct_x_count':len(grouped),
      'ordered_anchor_decimal_sha256':hashlib.sha256(anchors_text.encode()).hexdigest(),
      'cone_volume_interval':[str(volume.lo),str(volume.hi)],
      'cone_volume_upper':str(volume.hi),'D1_lower':str(d_lower),
      'tau_exact_decimal':str(tau),'score_upper':str(score.hi),
      'target_check_score_upper_lt_143_over_2':score.hi<core.Decimal(143)/core.Decimal(2),
      'linear_pareto_chain': {
        'scope':'Geometric linear interpolation of the old 65 anchors. The certified witness is exactly this interpolated chain cone, with all old 65 target-one constraints retained. Additional halo volume is not bounded.',
        'strict_x_increasing_and_y_decreasing':strict_chain,
        'segment_count':len(points)-1 if strict_chain else 0,
        'z_range_exact_decimal':chain_z_range,
        'volume_formula':'sum_segments dx*exp(-x0-y0)*F(dx+dy)+exp(-xlast-ylast), F(z)=(1-exp(-z))/z, F(0)=1.',
        'volume_interval':None if chain_volume is None else [str(chain_volume.lo),str(chain_volume.hi)],
        'score_upper':None if chain_score is None else str(chain_score.hi),
        'target_check_score_upper_lt_71_79':False if chain_score is None else chain_score.hi<core.Decimal('71.79')
      },
      'startup_sha256':start,'final_sha256':end,
      'core_source_sha256_expected':'17ef3f7fbd83211f7a411a6700ba751a94a1d5a85bb1c3b7dff9333d968bc3de'
    }
    if start['core']!=record['core_source_sha256_expected']: raise SystemExit('unexpected core hash')
    args.out.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
    if not record['target_check_score_upper_lt_143_over_2']:
        raise SystemExit('score bound did not beat 143/2')

if __name__=='__main__': main()
