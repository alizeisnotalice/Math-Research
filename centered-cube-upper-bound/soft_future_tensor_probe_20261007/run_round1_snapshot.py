from pathlib import Path
import argparse,importlib.util,json,math,time,hashlib
import numpy as np
D=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('model',D/'model.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def solve(point,rd):
    n=point['n'];N,LN,SN,BN={1:(128,17,65,17),2:(256,33,129,33),3:(512,65,257,65)}[rd]
    Ls=np.unique(np.concatenate((point['L_breakpoints'],np.geomspace(1,2,LN))))
    sigmas=np.unique(np.concatenate((np.linspace(0,min(1,32/n),SN),np.linspace(0,1,BN),np.minimum(1,2.**np.arange(5,math.ceil(math.log2(n))+1)/n))))
    Z=np.array(point['Z']);dx=np.array(point['dx']);t0=time.time()
    h,b,e,cur=m.coordinates(Z,dx,Ls,N)
    blo=np.maximum(0,b-e);bhi=b+e
    mid=[];low=[];high=[]
    for s in sigmas:
        mid.append(m.response(h,b,s));low.append(m.response(h,blo,s));high.append(m.response(h,bhi,s))
    mid=np.array(mid);low=np.array(low);high=np.array(high)
    gm=mid.max(axis=1);gl=low.max(axis=1);gh=high.max(axis=1)
    viol=np.where(gl>1)[0];eligible=np.where(gm>=1)[0]
    if len(eligible) and eligible[-1]<len(sigmas)-1:
        j=int(eligible[-1]);a=float(sigmas[j]);c=float(sigmas[j+1])
        for _ in range(32):
            ss=(a+c)/2
            if m.response(h,b,ss).max()>=1:a=ss
            else:c=ss
        raw=(a+c)/2;winner=int(np.argmax(m.response(h,b,raw)))
        root=dict(status='FINITE_L_LAST_VISIBLE_ROOT_ONLY',sigma_lo=a,sigma_hi=c,sigma_raw=raw,n_sigma=n*raw,Ls=float(Ls[winner]),index=winner)
    else:
        root=dict(status='NO_BRACKET_ENDPOINT_OR_NO_VISIBLE_ROOT',sigma_raw=None)
    # Numerically enclosed joint full-future cap, with analytic scalar quadrature
    # error but UNVERIFIED machine rounding. No formal continuous certificate.
    # L cells cover entire[1,2] and every hard arrival is a cut; upper at u is closed.
    LF=np.exp(n*np.log(Ls[1:]/Ls[:-1]))
    def cap_trial(s0):
        cuts=np.unique(np.concatenate(([s0],sigmas[sigmas>s0],[1.])))
        stack=[(float(a),float(c),0) for a,c in zip(cuts[:-1],cuts[1:])];tested=0;accepted=[];unknown=[]
        while stack and tested<384:
            a,c,depth=stack.pop();tested+=1
            ub=m.sigma_cell_upper(h,blo,bhi,a,c)
            ubfull=ub[1:]*LF
            largest=float(max(ub[0],ubfull.max()))
            if largest<=1:accepted.append(largest)
            elif depth<12 and n*(c-a)>1e-3:
                mm=(a+c)/2;stack.extend([(a,mm,depth+1),(mm,c,depth+1)])
            else:unknown.append(dict(sigma_cell=[a,c],upper=largest,L_cell_index=int(np.argmax(ubfull))))
        if stack:unknown.append(dict(status='CELL_BUDGET_REMAINING',cells=len(stack)))
        return dict(status='NUMERICALLY_ENCLOSED_FULL_FUTURE_NOT_FLOAT_CERTIFIED' if not unknown else 'UNKNOWN',sigma_start=s0,accepted_cells=len(accepted),tested_cells=tested,remaining_or_failed=unknown,max_accepted_upper=max(accepted) if accepted else None)
    trials=[];Krec=None
    if root['sigma_raw'] is not None:
        for dk in [.5,1.,2.,4.,8.,16.]:
            s0=min(1,root['sigma_hi']+dk/n)
            trial=cap_trial(s0);trial['n_sigma_gap']=dk;trials.append(trial)
            if trial['status'].startswith('NUMERICALLY_ENCLOSED'):break
        ix=root['index'];s=root['sigma_raw']
        Klo=m.weighted_future_average(h[:,ix,:],blo[:,ix,:],s);Kmid=m.weighted_future_average(h[:,ix,:],b[:,ix,:],s);Khi=m.weighted_future_average(h[:,ix,:],bhi[:,ix,:],s)
        Krec=dict(alpha=math.ceil(math.sqrt(n)),normalized_K_lower=Klo,normalized_K_mid=Kmid,normalized_K_upper=Khi,
            method='Exact positive normalized Bernstein coefficient integration in exact arithmetic; includes scalar quadrature bound, machine rounding not certified.',sigma=s,L=float(Ls[ix]),fullsoft_at_same_L=float(m.response(h,b,1.)[ix]))
    # Actual scalar hard jumps at each arrival retained. Isolated left trace has
    # full coordinate hard counts excluding all equal arrival atoms.
    traces=[]
    for L in point['L_breakpoints'][1:-1]:
        hc=[];lc=[]
        for q in m.Q:
            dd=2*np.abs(dx[:,None]+(Z[:,None]-(4//q)*np.arange(64*q)[None,:])/4)
            hc.append(np.sum(dd<=L,axis=1)/(q*L));lc.append(np.sum(dd<L,axis=1)/(q*L))
        ix=int(np.searchsorted(Ls,L));sr=root['sigma_raw'] if root['sigma_raw'] is not None else 0.
        rv=m.response(np.array(hc)[:,None,:],b[:,ix:ix+1,:],sr)[0];lv=m.response(np.array(lc)[:,None,:],b[:,ix:ix+1,:],sr)[0]
        traces.append([float(L),float(lv),float(rv)])
    np.savez_compressed(D/f'round{rd}_{point["id"]}.npz',L=Ls,sigma=sigmas,response_mid=mid,response_lower=low,response_upper=high,hard_jump_traces=np.array(traces),future_grid_mid=np.maximum.accumulate(gm[::-1])[::-1],future_grid_lower=np.maximum.accumulate(gl[::-1])[::-1],future_grid_upper=np.maximum.accumulate(gh[::-1])[::-1])
    return dict(point_id=point['id'],n=n,round=rd,spectral_N=N,L_nodes=len(Ls),hard_arrival_L_nodes=len(point['L_breakpoints']),all_hard_arrival_coverage=True,sigma_nodes=len(sigmas),seconds=time.time()-t0,max_coordinate_analytic_error_plus_unverified_guard=float(e.max()),max_coordinate_curvature_bound=float(cur.max()),
       finite_grid_max_F1=dict(lower=float(gl[-1]),mid=float(gm[-1]),upper=float(gh[-1])),finite_grid_endpoint_G1_half_condition=bool(gh[-1]<=.5),
       finite_visible_FIRST=root,future_cap_trials=trials,future_average=Krec,last_clear_finite_sigma_violation=float(sigmas[viol[-1]]) if len(viol) else None,
       candidate_middle_class=dict(threshold_1_over_n=1/n,S0=1/(16*math.sqrt(n)),raw_above_1_over_n=bool(root['sigma_raw']>1/n) if root['sigma_raw'] is not None else None,raw_above_S0=bool(root['sigma_raw']>1/(16*math.sqrt(n))) if root['sigma_raw'] is not None else None),
       scope='Full future numerical enclosure uses all L cells and adaptive sigma cells, but float function evaluation is not outward certified. Finite visible root not necessarily continuous actual FIRST. K/q belongs to that candidate, not a certified actual FIRST pair. No actualgeom gates or new weak-type fee.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--round',type=int,required=True);args=ap.parse_args();rd=args.round
    path=D/f'round{rd}.json';assert not path.exists(),'Immutable round refuses overwrite'
    points=json.loads((D/'frozen_points.json').read_text())['points'];rows=[];tic=time.time()
    for point in points:
        row=solve(point,rd);rows.append(row);(D/f'round{rd}_partial.json').write_text(json.dumps(rows,indent=2))
        print(json.dumps(dict(round=rd,id=row['point_id'],seconds=row['seconds'],root=row['finite_visible_FIRST'],K=row['future_average'],cap_status=[q['status'] for q in row['future_cap_trials']])),flush=True)
    path.write_text(json.dumps(dict(round=rd,rows=rows,elapsed_seconds=time.time()-tic),indent=2))
