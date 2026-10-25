"""Registered16 rational perturbations and2x240 finite superperiod refinements."""
from pathlib import Path
from fractions import Fraction as F
import time,json,numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
from round008_superperiod import evaluate,decode,capture
P=Path(__file__).resolve().parent;OLD=1.9017426111009554

def save_candidate(row,selected,stage):
    (P/'round008_candidate.json').write_text(json.dumps(dict(status='numerical candidate, needs exact R4 upper comparison',stage=stage,row=row,selected_keys=[list(k) for k in selected]),separators=(',',':'))+'\n')

def main():
    start=time.monotonic();embedding=json.loads((P/'round008_embedding_results.json').read_text());assert embedding['old_certificate_coefficients_exact_match'];reg=json.loads((P/'round008_preregister.json').read_text());rows=[];runs=[];stop='finite batch complete';candidate=None
    baseline,_,_=evaluate(np.zeros(9),capturing=False);assert abs(baseline['ratio']-OLD)<1e-12
    print('BASELINE',baseline,flush=True)
    try:
        cases=[]
        for k,spec in enumerate(reg['directions_and_amplitudes']):
            x=np.r_[np.zeros(4),float(F(spec['epsilon'])),list(map(lambda q:float(F(q)),spec['h']))];row,keys,selected=evaluate(x);row['case']=k;rows.append(row);cases.append((row,keys));print('CASE',k,row,flush=True)
            if row['ratio']>OLD+1e-11:save_candidate(row,selected,'structure');candidate=row;raise RuntimeError('positive numerical candidate: exact certification priority')
        # These16 seeds have distinct abs-mass/displacement types; omit the zero baseline.
        starts=sorted(cases,key=lambda q:-q[0]['ratio'])[:2]
        for initial,pool in starts:
            x=np.array(initial['x']);bounds=[(-.03,.03)]*4+[(-.2,.2)]+[(-.15,.15)]*4;simplex=np.tile(x,(10,1))
            for j in range(9):simplex[j+1,j]+= .002 if j<4 else(.01 if j==4 else.005)
            def obj(y):
                nonlocal candidate
                row,keys,selected=evaluate(y,pool=pool);row['case']='opt_'+str(initial['case']);rows.append(row)
                if len(rows)%25==0:print('OPT',len(rows),row,flush=True)
                if row['ratio']>OLD+1e-11:save_candidate(row,selected,'optimization');candidate=row;raise RuntimeError('positive numerical candidate: exact certification priority')
                return -row['ratio']
            res=minimize(obj,x,method='Nelder-Mead',bounds=bounds,options={'maxfev':240,'xatol':1e-8,'fatol':1e-10,'initial_simplex':simplex});run=dict(case=initial['case'],x=res.x.tolist(),ratio=float(-res.fun),nfev=res.nfev,stop=str(res.message));runs.append(run);print('RUN',run,flush=True)
            pp,qq,ph,ww=decode(res.x);pool=set(pool)|capture(qq,ph,ww,qmc.Sobol(4,scramble=True,seed=20261025808).random_base2(8));final,keys,selected=evaluate(res.x,pool=pool,power=12,seed=20261025809);final['case']='final_'+str(initial['case']);rows.append(final);print('FINAL POINT',final,flush=True)
            if final['ratio']>OLD+1e-11:save_candidate(final,selected,'final');candidate=final;raise RuntimeError('positive numerical candidate: exact certification priority')
    except RuntimeError as e:stop=str(e)
    result=dict(status='finite superperiod witnesses only; no input/family upper bound',stop=stop,baseline=baseline,rows=rows,runs=runs,candidate=candidate,elapsed=time.monotonic()-start)
    (P/'round008_search_results.json').write_text(json.dumps(result,indent=2)+'\n');print('END',stop,'best',max(r['ratio'] for r in rows),flush=True)
if __name__=='__main__':main()
