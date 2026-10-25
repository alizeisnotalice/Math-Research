"""Dual-satellite jointly partitioned exposed/contact model, finite structure pairs."""
from pathlib import Path
from itertools import product
import json,time
import numpy as np
from round005_exposure import A,geometry,prediction,phase_model
from round006_direction_breakpoints import candidates,axis_value
P=Path(__file__).resolve().parent

def base_phases():
    data=json.loads((P/'round006_screen_results.json').read_text())['scored'];sources=[[1,36,16,24],[36,36,16,24],[1,4,16,24]]
    return [np.array(next(r['phi'] for r in data if r['source_masses']==s)) for s in sources]
def reflect(p,indices):
    q=p.copy();q[indices]=1-q[indices];return q

def pairs():
    a,b,c=base_phases();swap=a[[0,1,3,2]]
    return [('A_complement',a,1-a),('A_R34',a,reflect(a,[2,3])),('A_R1',a,reflect(a,[0])),('A_B',a,b),('A_swap34',a,swap),('A_complementB',a,1-b),('B_complementC',b,1-c),('C_R2',c,reflect(c,[1])),('A_R2',a,reflect(a,[1])),('B_R34',b,reflect(b,[2,3]))]

def joint_model(phases,shares,patches):
    out=[[] for _ in range(4)]
    for p in patches:
        j=p['axis'];other=[k for k in range(4) if k!=j];ns=np.array(p['counts']);rr=ns**.25/2;RR=rr[:,None]/A[other];cuts=[]
        for k,axis in enumerate(other):
            cc=[p['lo'][k],p['hi'][k]]
            for phi in phases:
                for R in RR[:,k]:
                    for shift in range(-3,4):
                        for x in [phi[axis]+shift-R,phi[axis]+shift+R]:
                            if p['lo'][k]+1e-12<x<p['hi'][k]-1e-12:cc.append(x)
            cuts.append(np.unique(cc))
        ii=np.array(list(product(*(range(len(c)-1) for c in cuts))))
        mid=np.stack([(cuts[k][ii[:,k]]+cuts[k][ii[:,k]+1])/2 for k in range(3)],axis=1)
        weights=np.prod(np.stack([np.diff(cuts[k])[ii[:,k]] for k in range(3)],axis=1),axis=1)
        cg=np.zeros((len(mid),len(ns)));gp=cg.copy();gm=cg.copy();Rj=rr/A[j]
        for phi,share in zip(phases,shares):
            transverse=np.prod(np.maximum(0,np.floor(mid[:,None,:]+RR-phi[other]+1e-10)-np.ceil(mid[:,None,:]-RR-phi[other]-1e-10)+1),axis=2)
            lower=p['c']-Rj-phi[j];upper=p['c']+Rj-phi[j]
            closed=np.maximum(0,np.floor(upper+1e-10)-np.ceil(lower-1e-10)+1)
            nearplus=np.maximum(0,np.floor(upper-1e-10)-np.ceil(lower-1e-10)+1);nearminus=np.maximum(0,np.floor(upper+1e-10)-np.ceil(lower+1e-10)+1)
            factor=share/(8*ns**.75*A[j]);cg+=transverse*closed*factor;gp+=transverse*nearplus*factor;gm+=transverse*nearminus*factor
        mp=np.broadcast_to(np.array([n in p['plus'] for n in ns]),cg.shape);mm=np.broadcast_to(np.array([n in p['minus'] for n in ns]),cg.shape)
        fp=(cg>gp+1e-12)&mp;fm=(cg>gm+1e-12)&mm
        out[j].append((weights,cg,gp,gm,np.broadcast_to(Rj,cg.shape),mp,mm,fp,fm,np.full(len(weights),p['sign']==0)))
    packed=[]
    for parts in out:
        width=max(x[1].shape[1] for x in parts);cols=[[] for _ in range(10)]
        for part in parts:
            pad=width-part[1].shape[1]
            for k,x in enumerate(part):cols[k].append(x if k in [0,9] else np.pad(x,((0,0),(0,pad)),constant_values=(-np.inf if k in [1,2,3] else 0)))
        packed.append(tuple(np.concatenate(x) for x in cols))
    return packed

def optimize_direction(model,old):
    J=float(np.prod(A));direction=[];counts=[]
    for j in range(4):
        cuts,_=candidates(model[j]);scores=[axis_value(model[j],d,J,old) for d in cuts];i=int(np.argmax(scores));direction.append(cuts[i]);counts.append(len(cuts))
    return direction,prediction(model,np.array(direction),old)+1,counts

def main():
    st=time.monotonic();patches,_,geom=geometry();old=geom['baseline_joint_ratio'];a,_,_=base_phases();one=phase_model(a,patches);coincident=joint_model([a,a],[.2,.8],patches)
    checks=[np.zeros(4),np.ones(4)/4,np.array([.125,.1031116711436,0,0])]
    calibration=max(abs(prediction(one,d,old)-prediction(coincident,d,old)) for d in checks);assert calibration<1e-11
    rows=[];positive=None
    for name,phi1,phi2 in pairs():
        for share in [.8,.5,.2]:
            model=joint_model([phi1,phi2],[share,1-share],patches);direction,val,ncuts=optimize_direction(model,old)
            row=dict(name=name,phi1=phi1.tolist(),phi2=phi2.tolist(),share1=share,theta_physical=1,direction=direction,prediction=val,axis_breakpoints=ncuts);rows.append(row);print(row,flush=True)
            if val>1e-6:positive=row;(P/'round007_direction_packet.json').write_text(json.dumps(row,indent=2)+'\n');break
        if positive:break
    out=dict(status='joint three-orbit finite model, not sum of single-satellite predictions',rows=rows,positive=positive,coincident_phase_regression_max_error=calibration,elapsed=time.monotonic()-st)
    (P/'round007_model_results.json').write_text(json.dumps(out,indent=2)+'\n');print('FINAL',out['elapsed'],'positive',positive,flush=True)
if __name__=='__main__':main()
