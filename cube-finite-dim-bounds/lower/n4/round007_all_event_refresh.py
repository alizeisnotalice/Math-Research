"""Two fixed final points: all passing cube entry events at the same4096 centers."""
from pathlib import Path
import json,time,numpy as np
from scipy.stats import qmc
from round005_exposure import A,geometry
from round007_three_orbit import PRIMARY,PURE,PAIRS,face_keys,capture,boxes,volume
from round007_actual_search import OLD,save_candidate
P=Path(__file__).resolve().parent

def all_events(a,phases,weights,t):
    events=[]
    for orb in range(3):
        for j in range(4):events.extend(np.abs((t[:,j,None]-phases[orb,j]-np.arange(-2,4))*a[j]).T)
    rr=np.array(events).T;ll=np.ceil(t[:,None,None,:]-phases[None,None,:,:]-rr[:,:,None,None]/a-1e-11).astype(np.int16);hh=np.floor(t[:,None,None,:]-phases[None,None,:,:]+rr[:,:,None,None]/a+1e-11).astype(np.int16)
    counts=np.prod(np.maximum(0,hh-ll+1),axis=3);mass=counts@weights;good=(mass>np.maximum((2*rr)**4,1e-30)*(1+1e-10))&(rr<=2.5);keys=set()
    for i,j in np.argwhere(good):
        l=ll[i,j].copy();h=hh[i,j].copy()
        for orb in range(3):
            if np.any(h[orb]<l[orb]):l[orb]=0;h[orb]=-1
        keys.add(tuple(map(int,np.r_[l.ravel(),h.ravel()])))
    return keys,int(sum(good.ravel()))

def main():
    start=time.monotonic();data=json.loads((P/'round007_actual_results.json').read_text());patches,_,_=geometry();struct=[r for r in data['rows'] if not r['name'].startswith(('opt_','final_'))];chosen=[]
    for row in sorted(struct,key=lambda r:-r['increment_quotient']):
        if all(row['name']!=q['name'] or abs(row['share1']-q['share1'])>.01 for q in chosen):chosen.append(row)
        if len(chosen)==2:break
    finals=[r for r in data['rows'] if r['name'].startswith('final_')];results=[]
    for initial,final in zip(chosen,finals):
        st=time.monotonic();pool,_=face_keys([np.array(initial['phi1']),np.array(initial['phi2'])],patches)
        for row in [initial,final]:
            ph=np.array([np.zeros(4),row['phi1'],row['phi2']]);ww=np.array([1,row['sigma']*row['share1'],row['sigma']*(1-row['share1'])]);pool|=capture(np.array(row['a']),ph,ww,qmc.Sobol(4,scramble=True,seed=20261025707).random_base2(8))
        a=np.array(final['a']);ph=np.array([np.zeros(4),final['phi1'],final['phi2']]);ww=np.array([1,final['sigma']*final['share1'],final['sigma']*(1-final['share1'])]);extra,eventcount=all_events(a,ph,ww,qmc.Sobol(4,scramble=True,seed=20261025708).random_base2(12));keys=PRIMARY|PURE|PAIRS|pool|extra;b,selected=boxes(a,ph,ww,keys);ratio=float(volume(b)[0]/sum(ww))
        row={**final,'ratio':ratio,'increment_quotient':(ratio-OLD)/final['sigma'],'old_two_events_ratio':final['ratio'],'improvement':ratio-final['ratio'],'passing_entry_events':eventcount,'all_event_key_count':len(extra),'key_count':len(keys),'selected_boxes':len(b),'seconds':time.monotonic()-st};results.append(row);print(row,flush=True)
        if ratio>OLD+1e-11 and row['increment_quotient']>1e-5:save_candidate(row,selected,'all-event-fixed-final');raise RuntimeError('positive finite union: certificate priority')
    (P/'round007_all_event_results.json').write_text(json.dumps(dict(status='two fixed final points, same4096 centers, all passing events; no optimization',rows=results,elapsed=time.monotonic()-start),indent=2)+'\n')
if __name__=='__main__':main()
