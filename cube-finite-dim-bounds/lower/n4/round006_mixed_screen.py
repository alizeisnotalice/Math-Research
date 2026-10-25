"""Exact reflection-quotiented mixed phases, theta, and numerical face screening."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
import json,time,hashlib
import numpy as np
from scipy.optimize import minimize
from round005_exposure import A,geometry,phase_model,prediction,rvec
from round005_pure_satellite_check import AA,ZERO,interval,compare,add,scale,baseline_boxes,normal
P=Path(__file__).resolve().parent
SOURCES=[1,2,4,16,24,36]

def mixed_coordinates():
    choices=[]
    for j in range(4):
        dd={}
        for n in SOURCES:
            raw=scale(tuple(map(int,rvec(n))),2);ll,hh=interval(raw);al,ah=interval(AA[j]);ratios=[x/y for x,y in product([ll,hh],[al,ah])]
            ml=min(ratios).numerator//min(ratios).denominator;mh=max(ratios).numerator//max(ratios).denominator;assert ml==mh
            z=add(raw,scale(AA[j],-ml))
            if compare(scale(z,2),AA[j])>0:z=add(AA[j],scale(z,-1))
            dd.setdefault(z,[]).append(n)
        choices.append(sorted(dd.items(),key=lambda x:sum(interval(x[0]))))
    return choices

def all_theta(choices):
    boxes=baseline_boxes();periodic=[]
    for shift in product([-1,0,1],repeat=4):
        disp=[scale(AA[j],shift[j]) for j in range(4)]
        for lo,hi in boxes:periodic.append((tuple(add(lo[j],disp[j]) for j in range(4)),tuple(add(hi[j],disp[j]) for j in range(4))))
    ranks=[]
    for j in range(4):
        vals={x[0] for x in choices[j]}
        for lo,hi in periodic:vals.update([lo[j],hi[j]])
        vals=sorted(vals,key=lambda x:sum(interval(x)))
        for l,h in zip(vals,vals[1:]):assert interval(l)[1]<interval(h)[0]
        ranks.append({x:i for i,x in enumerate(vals)})
    lows=np.array([[ranks[j][lo[j]] for j in range(4)] for lo,hi in periodic]);highs=np.array([[ranks[j][hi[j]] for j in range(4)] for lo,hi in periodic]);rows=[]
    for combo in product(*choices):
        vectors=[x[0] for x in combo];z=np.array([ranks[j][vectors[j]] for j in range(4)]);active=np.all((lows<=z)&(highs>=z),axis=1)
        theta=normal.pure_satellite_term(tuple(map(int,z)),1,[(tuple(map(int,l)),tuple(map(int,h))) for l,h in zip(lows[active],highs[active])])
        phi=[float(sum(interval(v))/2)/A[j] for j,v in enumerate(vectors)]
        rows.append(dict(source_masses=[x[1][0] for x in combo],equivalent_sources=[x[1] for x in combo],phi=phi,theta=str(theta),active_old_images=int(sum(active))))
    return rows

def main():
    start=time.monotonic();choices=mixed_coordinates();rows=all_theta(choices)
    print('THETA',len(rows),[len(c) for c in choices],{x:sum(r['theta']==x for r in rows) for x in sorted(set(r['theta'] for r in rows))},flush=True)
    patches,_,geom=geometry();old=geom['baseline_joint_ratio'];directions=[np.zeros(4),np.ones(4)/4]+[np.eye(4)[j]*s/4 for j in range(4) for s in [-1,1]]
    # Positive-theta and mixed-source points first; this changes order, never duplicates a phase.
    rows.sort(key=lambda r:(-float(F(r['theta'])),len(set(r['source_masses']))*-1,tuple(r['source_masses'])))
    scored=[];runs=[];stop='finite screening complete';candidate=None
    for i,row in enumerate(rows):
        st=time.monotonic();model=phase_model(np.array(row['phi']),patches);scores=[prediction(model,d,old)+float(F(row['theta'])) for d in directions];k=int(np.argmax(scores))
        item=dict(**row,predictions=scores,best_prediction=scores[k],direction=directions[k].tolist(),seconds=time.monotonic()-st);scored.append(item)
        if i%16==0:print('SCORE',i,item,flush=True)
        if scores[k]>1e-6:candidate=item;stop='positive numerical direction; stop for finite epsilon certificate';break
    if candidate is None:
        for row in sorted(scored,key=lambda x:-x['best_prediction'])[:2]:
            model=phase_model(np.array(row['phi']),patches);theta=float(F(row['theta']))
            res=minimize(lambda d:-(prediction(model,d,old)+theta),row['direction'],method='Nelder-Mead',bounds=[(-.5,.5)]*4,options={'maxfev':600,'xatol':1e-8,'fatol':1e-9})
            rr=dict(source_masses=row['source_masses'],phi=row['phi'],theta=row['theta'],direction=res.x.tolist(),prediction=float(-res.fun),nfev=res.nfev,stop=str(res.message));runs.append(rr);print('RUN',rr,flush=True)
            if -res.fun>1e-6:candidate=rr;stop='positive optimized numerical direction; stop for finite certificate';break
    out=dict(status='numerical mixed-phase screening plus exact theta; not family bound',stop=stop,unique_coordinate_counts=[len(c) for c in choices],total_mixed_phases=len(rows),theta_histogram={x:sum(r['theta']==x for r in rows) for x in sorted(set(r['theta'] for r in rows))},scored=scored,runs=runs,candidate=candidate,geometry=geom,elapsed=time.monotonic()-start)
    (P/'round006_screen_results.json').write_text(json.dumps(out,indent=2)+'\n')
    if candidate is not None:(P/'round006_direction_packet.json').write_text(json.dumps(candidate,indent=2)+'\n')
    print('FINAL',stop,'elapsed',out['elapsed'],flush=True)
if __name__=='__main__':main()
