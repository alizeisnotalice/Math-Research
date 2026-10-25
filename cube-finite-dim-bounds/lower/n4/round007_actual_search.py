"""Finite dual-satellite actual structure batch and14-variable refinement."""
from pathlib import Path
import time,json
import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc
from round005_exposure import A,geometry
from round007_joint_model import pairs
from round007_three_orbit import evaluate,face_keys,capture
P=Path(__file__).resolve().parent;OLD=1.9017426111009554

def save_candidate(row,selected,stage):
    (P/'round007_actual_candidate.json').write_text(json.dumps(dict(status='numerical candidate; exact certificate needed',stage=stage,row=row,selected_keys=[list(k) for k in selected]),separators=(',',':'))+'\n')

def main():
    start=time.monotonic();data=json.loads((P/'round007_model_results.json').read_text());patches,_,_=geometry();pools={};rows=[];runs=[];stop='finite batch complete';candidate=None
    try:
        if data['positive']:
            structures=[data['positive']];sigmas=[1e-5,1e-6,1e-7]
        else:structures=data['rows'];sigmas=[.0001,.001,.01]
        for item in structures:
            name=item['name'];phi1=np.array(item['phi1']);phi2=np.array(item['phi2']);d=np.array(item['direction']);share=item['share1']
            if name not in pools:
                pools[name],meta=face_keys([phi1,phi2],patches);print('POOL',name,meta,flush=True)
            for sigma in sigmas:
                row,keys,selected=evaluate(A*(1+sigma*d),phi1,phi2,sigma,share,pool=pools[name]);row['name']=name;rows.append(row)
                print('STRUCT',len(rows),row,flush=True)
                if row['ratio']>OLD+1e-11 and row['increment_quotient']>1e-5:save_candidate(row,selected,'structure');candidate=row;raise RuntimeError('positive finite union: certificate priority')
        # Best two independent structures, finite refinement; retain discovered structure keys.
        chosen=[]
        for row in sorted(rows,key=lambda r:-r['increment_quotient']):
            if all(row['name']!=q['name'] or abs(row['share1']-q['share1'])>.01 for q in chosen):chosen.append(row)
            if len(chosen)==2:break
        for initial in chosen:
            sigma=initial['sigma'];d=(np.array(initial['a'])/A-1)/sigma;phi1=np.array(initial['phi1']);phi2=np.array(initial['phi2']);share=initial['share1'];pool=set(pools[initial['name']])
            ph=np.array([np.zeros(4),phi1,phi2]);weights=np.array([1,sigma*share,sigma*(1-share)])
            pool|=capture(np.array(initial['a']),ph,weights,qmc.Sobol(4,scramble=True,seed=20261025707).random_base2(8))
            x=np.r_[d,phi1,phi2,sigma,share];bounds=[(-.5,.5)]*4+[(0,1)]*8+[(1e-7,.02),(.05,.95)]
            simplex=np.tile(x,(15,1))
            for j in range(14):
                step=.03 if j<4 else(.01 if j<12 else(.1*sigma if j==12 else.05))
                simplex[j+1,j]=min(bounds[j][1],x[j]+step)
                if simplex[j+1,j]==x[j]:simplex[j+1,j]=max(bounds[j][0],x[j]-step)
            def obj(y):
                nonlocal candidate
                a=A*(1+y[12]*y[:4]);row,keys,selected=evaluate(a,y[4:8],y[8:12],y[12],y[13],pool=pool);row['name']='opt_'+initial['name'];rows.append(row)
                if len(rows)%20==0:print('OPT',len(rows),row,flush=True)
                if row['ratio']>OLD+1e-11 and row['increment_quotient']>1e-5:save_candidate(row,selected,'optimization');candidate=row;raise RuntimeError('positive finite union: certificate priority')
                return -row['increment_quotient']
            res=minimize(obj,x,method='Nelder-Mead',bounds=bounds,options={'maxfev':240,'xatol':1e-8,'fatol':1e-8,'initial_simplex':simplex});r=dict(name=initial['name'],x=res.x.tolist(),nfev=res.nfev,quotient=float(-res.fun),stop=str(res.message));runs.append(r);print('RUN',r,flush=True)
            y=res.x;pool|=capture(A*(1+y[12]*y[:4]),np.array([np.zeros(4),y[4:8],y[8:12]]),np.array([1,y[12]*y[13],y[12]*(1-y[13])]),qmc.Sobol(4,scramble=True,seed=20261025707).random_base2(8))
            row,keys,selected=evaluate(A*(1+y[12]*y[:4]),y[4:8],y[8:12],y[12],y[13],pool=pool,power=12,seed=20261025708);row['name']='final_'+initial['name'];rows.append(row);print('FINAL REFRESH',row,flush=True)
            if row['ratio']>OLD+1e-11 and row['increment_quotient']>1e-5:save_candidate(row,selected,'final');candidate=row;raise RuntimeError('positive finite union: certificate priority')
    except RuntimeError as e:stop=str(e)
    result=dict(status='actual three-orbit finite joint unions, not input/family bounds',stop=stop,rows=rows,runs=runs,candidate=candidate,elapsed=time.monotonic()-start)
    (P/'round007_actual_results.json').write_text(json.dumps(result,indent=2)+'\n');print('FINAL',stop,'best',max(r['ratio'] for r in rows),flush=True)
if __name__=='__main__':main()
