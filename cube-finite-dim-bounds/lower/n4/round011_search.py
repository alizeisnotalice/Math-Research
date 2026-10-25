"""Eight registered coupled chess directions; one gated finite refinement."""
from pathlib import Path
from fractions import Fraction as F
import json,time,numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc
from round011_chess_engine import evaluate,decode,capture
P=Path(__file__).resolve().parent;OLD=1.9017426111009554

def save_candidate(row,selected,stage):
    (P/'round011_candidate.json').write_text(json.dumps(dict(status='numerical only; exact R4 upper comparison/audit required',stage=stage,row=row,selected_keys=[list(k) for k in selected]),separators=(',',':'))+'\n')

def main():
    st=time.monotonic();embedding=json.loads((P/'round011_embedding_results.json').read_text());assert embedding['old16_coefficients_exact_match'];reg=json.loads((P/'round011_preregister.json').read_text());rows=[];run=None;candidate=None;stop='eight directions complete; no capture gain, refinement skipped'
    baseline,_,_=evaluate(np.zeros(9),capturing=False);assert abs(baseline['ratio']-OLD)<1e-12;print('BASELINE',baseline,flush=True);structures=[]
    try:
        for k,spec in enumerate(reg['directions']):
            x=np.r_[np.zeros(4),float(F(spec['epsilon'])),list(map(lambda q:float(F(q)),spec['h']))];old,_,_=evaluate(x,capturing=False);row,keys,selected=evaluate(x);row.update(case=k,embedded_only_ratio=old['ratio'],capture_gain=row['ratio']-old['ratio']);rows.append(row);structures.append((row,keys));print('STRUCT',k,row,flush=True)
            if row['ratio']>OLD+1e-11:save_candidate(row,selected,'structure');candidate=row;raise RuntimeError('new numerical candidate: exact certificate priority')
        benefited=[q for q in structures if q[0]['capture_gain']>1e-10]
        if benefited:
            initial,pool=max(benefited,key=lambda q:q[0]['ratio']);x=np.array(initial['x']);bounds=[(-.03,.03)]*4+[(-.2,.2)]+[(-.1,.1)]*4;simplex=np.tile(x,(10,1))
            for j in range(9):simplex[j+1,j]+=.001 if j<4 else(.01 if j==4 else.002)
            def obj(y):
                nonlocal candidate
                row,keys,selected=evaluate(y,pool=pool);row['case']='opt';rows.append(row)
                if len(rows)%25==0:print('OPT',len(rows),row,flush=True)
                if row['ratio']>OLD+1e-11:save_candidate(row,selected,'optimization');candidate=row;raise RuntimeError('new numerical candidate: exact certificate priority')
                return -row['ratio']
            result=minimize(obj,x,method='Nelder-Mead',bounds=bounds,options={'maxfev':320,'initial_simplex':simplex,'xatol':1e-8,'fatol':1e-10});run=dict(start_case=initial['case'],x=result.x.tolist(),ratio=float(-result.fun),nfev=result.nfev,stop=str(result.message));print('RUN',run,flush=True)
            pp,qq,ph,ww=decode(result.x);pool=set(pool)|capture(qq,ph,ww,qmc.Sobol(4,scramble=True,seed=20261025111).random_base2(8));final,keys,selected=evaluate(result.x,pool=pool,power=12,seed=20261025112);final['case']='final';rows.append(final);print('FINAL POINT',final,flush=True);stop='eight directions, capture-gated refinement and final capture complete'
            if final['ratio']>OLD+1e-11:save_candidate(final,selected,'final');candidate=final;raise RuntimeError('new numerical candidate: exact certificate priority')
    except RuntimeError as error:stop=str(error)
    output=dict(status='finite chess witnesses only; no input/family upper bound',stop=stop,baseline=baseline,rows=rows,run=run,candidate=candidate,elapsed=time.monotonic()-st)
    (P/'round011_search_raw.json').write_text(json.dumps(output,indent=2)+'\n');print('END',stop,'best',max(r['ratio'] for r in rows),flush=True)
if __name__=='__main__':main()
