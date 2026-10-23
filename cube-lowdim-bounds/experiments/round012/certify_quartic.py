#!/usr/bin/env python3
"""Finite strict quartic certification, independent of the numerical LP budget.

Reuse only round008 Fraction interval primitives, logarithms, xlogx, tail,
and b2 diagonal; quartic diagonal, C4 correction, small radius, and cost are
implemented here explicitly. At most two candidate attempts, each capped at
285s and 100000 accepted+split boxes. All accepted leaves are saved.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
R=Q(379,200)
INCUMBENT=Q(6938554331,1790000000)
MAX_ATTEMPTS=2
MAX_BOXES=100000
SECONDS=285.
CANDIDATES=[(Q(3,4),Q(21,100)),(Q(3,4),Q(1,5))]

def find_core(explicit):
    if explicit:return explicit.resolve()
    paths=[ROOT.parent/'round008'/'certify.py',
           ROOT.parent.parent/'Math-Research-lowdim-agent'/'cube-lowdim-bounds'/'experiments'/'round008'/'certify.py']
    for p in paths:
        if p.exists():return p.resolve()
    raise FileNotFoundError('Pass --core <round008/certify.py> after relocation.')

def load_core(path):
    spec=importlib.util.spec_from_file_location('round008_fraction_for_quartic',path)
    core=importlib.util.module_from_spec(spec);sys.modules[spec.name]=core
    spec.loader.exec_module(core);core.R=R
    return core

def certify(core,e2,e4):
    I=core.I
    def cube(x):return x*x*x
    def b4_diag(r):
        s=1-r;a=1-r/2
        return 8*(a*r.square()+Q(22,3)*cube(a)
                  +2*(1+cube(s))*core.log_i(r/2)
                  -2*s.square()*core.xlogx(s))
    def c4(r):
        if r.hi<=Q(2,3):
            q=-Q(5,2)*r-3+18/r-4/r.square()
        elif r.lo>=Q(2,3):
            q=(165-101*r)/8+1/(2*r.square())
        else:raise ValueError('C4 box straddles 2/3')
        return 12*core.log_i(r/2)+q
    def lower(r):
        tail=2*(2*core.tail_series(r)+I.point(1)/R)/(R-1)
        g=tail-e2*core.b2_diag(r)-e4*b4_diag(r)
        s=1-r
        k=e2*64*s.square()/(2-r).square()/(2-r)+I.point(1)/(R**3*(R-1))-4*e4*c4(r)
        correction=min(k.lo,Q(0))*r.hi*r.hi
        return g+correction
    start=monotonic();todo=[I(Q(1,10),Q(2,3)),I(Q(2,3),Q(1))]
    leaves=[];unresolved=[];splits=0;minimum=None
    while todo:
        if len(leaves)+splits>=MAX_BOXES or monotonic()-start>SECONDS:
            unresolved.extend(todo);break
        box=todo.pop();bound=lower(box)
        if bound.lo>=0:
            leaves.append({'lo':str(box.lo),'hi':str(box.hi),'lower':str(bound.lo)})
            minimum=bound.lo if minimum is None else min(minimum,bound.lo)
        elif box.hi-box.lo<=Q(1,2**80):unresolved.append(box)
        else:
            mid=(box.lo+box.hi)/2
            todo.extend([I(box.lo,mid),I(mid,box.hi)]);splits+=1
        if (len(leaves)+splits)%500==0:
            print(f'accepted={len(leaves)} splits={splits} pending={len(todo)}',flush=True)
    leaves.sort(key=lambda row:Q(row['lo']))
    cover=bool(leaves) and not unresolved and Q(leaves[0]['lo'])==Q(1,10) and Q(leaves[-1]['hi'])==1
    cover=cover and all(Q(a['hi'])==Q(b['lo']) for a,b in zip(leaves,leaves[1:]))
    text='\n'.join(f"{x['lo']}/{x['hi']}" for x in leaves)
    inner={'status':'certified' if cover else 'unresolved','accepted_boxes':len(leaves),
           'splits':splits,'unresolved_boxes':len(unresolved),'cover_verified':cover,
           'minimum_accepted_lower':str(minimum) if minimum is not None else None,
           'minimum_accepted_lower_display':float(minimum) if minimum is not None else None,
           'accepted_leaf_sha256':sha256(text.encode()).hexdigest(),
           'unresolved_intervals':[[str(x.lo),str(x.hi)] for x in unresolved[:20]],
           'elapsed_seconds':monotonic()-start}
    # Recompute the supplied small_radius.json quartic constant with exact logs.
    h4=8*(Q(553,150)+Q(729,1000)*core.log_point(Q(1,9))+core.log_point(Q(1,11)))
    t2=3-5*core.log_point(Q(2))
    expected=json.loads((ROOT/'small_radius.json').read_text())
    small={'range':'0<r<=1/10','h4_upper_interval':h4.pair(),
           'h4_upper_display':float(h4.hi),'h4_strictly_negative':h4.hi<0,
           'matches_supplied_small_radius_interval':h4.pair()==expected['h4_uniform_upper_interval'],
           'quadratic_T_half_upper':str(t2.hi),'quadratic_T_half_strictly_negative':t2.hi<0,
           'analytic_dependency':'t=u/(1-v)<=r; h4 bound in small_radius.json; existing T monotonicity for quadratic'}
    total=e2+e4;b=2-2*core.log_point(2)
    outer=1-total*R*(R-1)*b
    dominance=R-total*(R-1)/4
    cost=R*R/(R-1)-e2/6-e4/15
    checks={'outer_slack':outer.pair(),'outer_slack_display':[float(outer.lo),float(outer.hi)],
            'kernel_dominance_margin':str(dominance),'kernel_dominance_margin_display':float(dominance),
            'cost':str(cost),'cost_display':float(cost),'previous_upper':str(INCUMBENT),
            'improvement':str(INCUMBENT-cost),'improvement_display':float(INCUMBENT-cost),
            'strict':outer.lo>0 and dominance>0 and cost<INCUMBENT and e2>=0 and e4>=0}
    passed=(cover and small['h4_strictly_negative'] and small['quadratic_T_half_strictly_negative']
            and small['matches_supplied_small_radius_interval'] and checks['strict'])
    return {'status':'certified_under_quartic_angular_and_transfer' if passed else 'not_closed',
            'parameters':{'R':str(R),'e2':str(e2),'e3':'0','e4':str(e4)},
            'inner':inner,'small_radius':small,'outer_and_cost':checks,'accepted_leaves':leaves}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--core',type=Path);args=ap.parse_args()
    path=ROOT/'certificate.json'
    data=json.loads(path.read_text()) if path.exists() else {'round':12,'attempts':[]}
    if any(a.get('status')=='certified_under_quartic_angular_and_transfer' for a in data['attempts']):
        raise SystemExit('A successful strict certificate already exists; stopping without another attempt.')
    if len(data['attempts'])>=MAX_ATTEMPTS:raise SystemExit('Two-candidate strict budget exhausted.')
    number=len(data['attempts'])+1;e2,e4=CANDIDATES[number-1]
    data['attempts'].append({'attempt':number,'status':'started','parameters':{'R':str(R),'e2':str(e2),'e3':'0','e4':str(e4)}})
    path.write_text(json.dumps(data,indent=2)+'\n')
    source=find_core(args.core);core=load_core(source)
    result=certify(core,e2,e4);result['attempt']=number
    data['attempts'][-1]=result
    data.update({'status':result['status'],'core_dependency':'round008/certify.py (override with --core)',
                 'core_sha256':sha256(source.read_bytes()).hexdigest(),
                 'certifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
                 'analytic_dependencies':{p.name:sha256(p.read_bytes()).hexdigest()
                    for p in [ROOT/'quartic_angular.md',ROOT/'small_radius.json',ROOT/'transfer_extension.md'] if p.exists()},
                 'limits':{'max_candidate_attempts':MAX_ATTEMPTS,'max_boxes_per_attempt':MAX_BOXES,
                           'seconds_per_attempt':SECONDS,'tail_terms':28,'log_terms':72},
                 'b4_diagonal_check':'Independent SymPy 1.14.0 check: two-slice minus b4diag simplifies to 0; polynomial forms also differ by 0.',
                 'method':'Fraction intervals; full-r quartic curvature correction; initial split at 2/3; all accepted leaf bounds saved.'})
    path.write_text(json.dumps(data,indent=2)+'\n')
    summary={k:v for k,v in result.items() if k!='accepted_leaves'}
    print(json.dumps(summary,indent=2),flush=True)
    if result['status']=='not_closed':raise SystemExit(1)
if __name__=='__main__':main()
