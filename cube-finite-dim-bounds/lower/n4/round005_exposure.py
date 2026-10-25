"""Critical-phase exposed-face model; numerical directional predictions, not a proof."""
from pathlib import Path
from itertools import product
import json,time,hashlib
import numpy as np
from scipy.optimize import minimize
from staggered_search import A
P=Path(__file__).resolve().parent
U=2**.25;V=3**.25
AV=np.zeros((4,16),int);AV[:,0]=[4,2,2,2];AV[0,8]=2;AV[1,13]=1;AV[2,4]=AV[3,4]=2
BASIS=np.array([U**i*V**j/4 for i in range(4) for j in range(4)])
def rvec(n):
    p=q=0
    while n%2==0:p+=1;n//=2
    while n%3==0:q+=1;n//=3
    assert n==1
    z=np.zeros(16,int);z[4*(p%4)+q%4]=2*2**(p//4)*3**(q//4);return z

def geometry():
    starts=np.array(list(product(range(-2,2),repeat=4)));raw=[];coords=[{tuple(np.zeros(16,int)),tuple(AV[j])} for j in range(4)]
    for ell in product(range(1,5),repeat=4):
        n=int(np.prod(ell));rv=rvec(n);r=float(rv@BASIS)
        for k in starts:
            lo=k+np.array(ell)-1-r/A;hi=k+r/A
            if np.any(np.minimum(hi,1)<=np.maximum(lo,0)):continue
            lk=[];hk=[]
            for j in range(4):
                l=(k[j]+ell[j]-1)*AV[j]-rv;h=k[j]*AV[j]+rv
                lk.append(tuple(l) if lo[j]>0 else tuple(np.zeros(16,int)))
                hk.append(tuple(h) if hi[j]<1 else tuple(AV[j]));coords[j].update([lk[-1],hk[-1]])
            raw.append((tuple(lk),tuple(hk),n,r))
    raw=list(dict.fromkeys(raw));ordered=[sorted(c,key=lambda k:np.array(k)@BASIS) for c in coords]
    cuts=[np.array([np.array(k)@BASIS/A[j] for k in ordered[j]]) for j in range(4)]
    maps=[{k:i for i,k in enumerate(c)} for c in ordered]
    li=np.array([[maps[j][b[0][j]] for j in range(4)] for b in raw]);hi=np.array([[maps[j][b[1][j]] for j in range(4)] for b in raw]);ns=np.array([b[2] for b in raw]);rs=np.array([b[3] for b in raw])
    patches=[];area={}
    for axis in range(4):
        others=[j for j in range(4) if j!=axis];ii=np.array(list(product(*(range(len(cuts[j])-1) for j in others))))
        lower=np.stack([cuts[j][ii[:,k]] for k,j in enumerate(others)],axis=1);upper=np.stack([cuts[j][ii[:,k]+1] for k,j in enumerate(others)],axis=1)
        cover=np.all((li[:,None,others]<=ii)&(hi[:,None,others]>ii),axis=2)
        for c in range(1,len(cuts[axis])-1):
            left=np.any(cover&((li[:,axis]<c)&(hi[:,axis]>=c))[:,None],axis=0)
            right=np.any(cover&((li[:,axis]<=c)&(hi[:,axis]>c))[:,None],axis=0)
            crossing=np.any(cover&((li[:,axis]<c)&(hi[:,axis]>c))[:,None],axis=0)
            relevant=(left^right)|(left&right&~crossing)
            for q in np.flatnonzero(relevant):
                sign=0 if left[q] and right[q] else (1 if left[q] else -1)
                plus=sorted(set(map(int,ns[(hi[:,axis]==c)&cover[:,q]]))) if left[q] else []
                minus=sorted(set(map(int,ns[(li[:,axis]==c)&cover[:,q]]))) if right[q] else []
                unique=sorted(set(plus+minus))
                patch=dict(axis=axis,sign=sign,c=float(cuts[axis][c]),lo=lower[q],hi=upper[q],counts=unique,plus=plus,minus=minus,cell=tuple(map(int,ii[q])))
                patches.append(patch)
                for n in unique:area[n]=area.get(n,0)+float(np.prod(upper[q]-lower[q]))
    # Baseline volume directly from the exact arrangement cells, independent of union recursion.
    ii=np.array(list(product(*(range(len(c)-1) for c in cuts))));covered=np.any(np.all((li[:,None,:]<=ii)&(hi[:,None,:]>ii),axis=2),axis=0)
    cellvol=np.prod(np.stack([np.diff(cuts[j])[ii[:,j]] for j in range(4)],axis=1),axis=1)
    old=float(np.prod(A)*sum(cellvol[covered]));assert abs(old-1.9017426111009554)<1e-12
    raw_patch_count=len(patches);groups={}
    for p in patches:groups.setdefault((p['axis'],p['sign'],p['c'],tuple(p['plus']),tuple(p['minus'])),set()).add(p['cell'])
    compressed=[]
    for (axis,sign,c,plus,minus),cells in groups.items():
        counts=sorted(set(plus+minus))
        others=[j for j in range(4) if j!=axis]
        while cells:
            lo=min(cells);hi=list(lo)
            for j in range(3):
                while True:
                    trial=hi[j]+1
                    layer=list(product(*(range(lo[k],hi[k]+1) if k!=j else [trial] for k in range(3))))
                    if all(x in cells for x in layer):hi[j]=trial
                    else:break
            for x in product(*(range(lo[k],hi[k]+1) for k in range(3))):cells.remove(x)
            compressed.append(dict(axis=axis,sign=sign,c=c,lo=np.array([cuts[j][lo[k]] for k,j in enumerate(others)]),hi=np.array([cuts[j][hi[k]+1] for k,j in enumerate(others)]),counts=list(counts),plus=list(plus),minus=list(minus)))
    return compressed,area,dict(raw_witness_boxes=len(raw),cut_counts=list(map(len,cuts)),raw_exposed_patches=raw_patch_count,exposed_patches=sum(p['sign']!=0 for p in compressed),internal_contact_patches=sum(p['sign']==0 for p in compressed),total_patches=len(compressed),baseline_joint_ratio=old,min_distinct_normalized_endpoint_gap=min(float(np.diff(c).min()) for c in cuts))

def phase_model(phi,patches):
    out=[[] for _ in range(4)]
    for p in patches:
        j=p['axis'];other=[k for k in range(4) if k!=j];ns=np.array(p['counts']);rr=ns**.25/2;RR=rr[:,None]/A[other]
        cuts=[]
        for k,axis in enumerate(other):
            cc=[p['lo'][k],p['hi'][k]]
            for R in RR[:,k]:
                for shift in range(-3,4):
                    for x in [phi[axis]+shift-R,phi[axis]+shift+R]:
                        if p['lo'][k]+1e-12<x<p['hi'][k]-1e-12:cc.append(x)
            cuts.append(np.unique(cc))
        ii=np.array(list(product(*(range(len(c)-1) for c in cuts))))
        mid=np.stack([(cuts[k][ii[:,k]]+cuts[k][ii[:,k]+1])/2 for k in range(3)],axis=1)
        weights=np.prod(np.stack([np.diff(cuts[k])[ii[:,k]] for k in range(3)],axis=1),axis=1)
        transverse=np.prod(np.maximum(0,np.floor(mid[:,None,:]+RR-phi[other]+1e-10)-np.ceil(mid[:,None,:]-RR-phi[other]-1e-10)+1),axis=2)
        Rj=rr/A[j];lower=p['c']-Rj-phi[j];upper=p['c']+Rj-phi[j]
        closed=np.maximum(0,np.floor(upper+1e-10)-np.ceil(lower-1e-10)+1)
        nearplus=np.maximum(0,np.floor(upper-1e-10)-np.ceil(lower-1e-10)+1)
        nearminus=np.maximum(0,np.floor(upper+1e-10)-np.ceil(lower+1e-10)+1)
        factor=1/(8*ns**.75*A[j]);cgain=transverse*closed*factor
        gp=transverse*nearplus*factor;gm=transverse*nearminus*factor
        mp=np.broadcast_to(np.array([n in p['plus'] for n in ns]),cgain.shape)
        mm=np.broadcast_to(np.array([n in p['minus'] for n in ns]),cgain.shape)
        fp=(closed>nearplus)[None,:]&(transverse>0)&mp
        fm=(closed>nearminus)[None,:]&(transverse>0)&mm
        out[j].append((weights,cgain,gp,gm,np.broadcast_to(Rj,cgain.shape),mp,mm,fp,fm,np.full(len(weights),p['sign']==0)))
    packed=[]
    for parts in out:
        width=max(x[1].shape[1] for x in parts);cols=[[] for _ in range(10)]
        for part in parts:
            pad=width-part[1].shape[1]
            for k,x in enumerate(part):
                if k in [0,9]:cols[k].append(x)
                else:cols[k].append(np.pad(x,((0,0),(0,pad)),constant_values=(-np.inf if k in [1,2,3] else 0)))
        packed.append(tuple(np.concatenate(x) for x in cols))
    return packed

def prediction(model,d,old):
    flux=0
    for j in range(4):
        weights,cg,gp,gm,R,mp,mm,fp,fm,seam=model[j]
        plus=np.max(np.where(mp,gp-R*d[j],-np.inf),axis=1)
        minus=np.max(np.where(mm,gm-R*d[j],-np.inf),axis=1)
        born=np.maximum(0,np.max(np.where(fp|fm,cg-R*d[j],-np.inf),axis=1))
        inward=np.maximum(plus,minus)
        gain=inward+np.maximum(0,born-np.maximum(inward,-born))
        # Old internal contacts can open gaps; only pre-existing crossing boxes shield them.
        L=plus;right=-minus
        gap=np.maximum(0,right-L)-np.maximum(0,np.minimum(right,born)-np.maximum(L,-born))
        gain=np.where(seam,-gap,gain)
        flux+=float(weights@gain)
    return float(np.prod(A)*flux+old*sum(d)-old)

def main():
    start=time.monotonic();patches,areas,geom=geometry();print('GEOMETRY',geom,flush=True)
    old=geom['baseline_joint_ratio'];phases=[np.zeros(4)];phase_specs=[dict(kind='coincident')]
    for n in sorted(areas,key=lambda n:-areas[n]):
        rr=n**.25/2
        for signs in product([-1,1],repeat=4):
            phi=np.mod(np.array(signs)*2*rr/A,1)
            if not any(np.max(abs(phi-q))<1e-10 for q in phases):
                phases.append(phi);phase_specs.append(dict(kind='critical',primary_mass=n,signs=list(signs),formula='phi_j=frac(sign_j*N0^(1/4)/a_j)'))
            if len(phases)>=97:break
        if len(phases)>=97:break
    directions=[np.zeros(4),np.ones(4)/4]+[np.eye(4)[j]*s/4 for j in range(4) for s in [-1,1]]
    rows=[];models={};stop='finite batch completed'
    for k,phi in enumerate(phases):
        st=time.monotonic();model=phase_model(phi,patches);models[k]=model
        scores=[prediction(model,d,old) for d in directions];best=int(np.argmax(scores))
        row=dict(phase_index=k,phase_formula=phase_specs[k],phi=phi.tolist(),predictions=scores,best_direction=directions[best].tolist(),best_prediction=scores[best],model_seconds=time.monotonic()-st);rows.append(row)
        if k%8==0:print(k,row,flush=True)
        if row['model_seconds']>5:stop='phase model >5seconds finite stop';break
    calibration=prediction(models[0],np.ones(4)/4,old);assert abs(calibration)<1e-9,calibration
    runs=[]
    for row in sorted(rows[1:],key=lambda r:-r['best_prediction'])[:2]:
        model=models[row['phase_index']]
        res=minimize(lambda d:-prediction(model,d,old),row['best_direction'],method='Nelder-Mead',bounds=[(-.5,.5)]*4,options={'maxfev':600,'xatol':1e-8,'fatol':1e-9})
        runs.append(dict(phase_formula=row['phase_formula'],phi=row['phi'],direction=res.x.tolist(),prediction=float(-res.fun),nfev=res.nfev,stop=str(res.message)))
        print('RUN',runs[-1],flush=True)
    out=dict(status='numerical exposed-face directional predictions only; not rigorous derivatives',geometry=geom,stop=stop,phase_count=len(rows),coincident_scaled_calibration=calibration,phases=rows,runs=runs,elapsed=time.monotonic()-start)
    (P/'round005_exposure_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('FINAL',stop,'best',max(r['prediction'] for r in runs),flush=True)
if __name__=='__main__':main()
