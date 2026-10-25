"""One24-point Sobol batch followed by at most32 deterministic rational proposals."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,os,sys,warnings
import numpy as np
import scipy
from scipy.stats import qmc
from flint import ctx
from kernel import GreenKernel,sha,TABLE
P=Path(__file__).resolve().parent;CACHE=Path('/tmp/astra_green018')
DEN=104729;SEED=202610260181
def canonical(obj):return json.dumps(obj,sort_keys=True,separators=(',',':')).encode()
def save(path,obj):
    tmp=Path(str(path)+'.partial');tmp.write_bytes(canonical(obj));os.replace(tmp,path)
def point(tail):return (F(1),)+tuple(F(n,DEN) for n in tail)
def validate(tail,engine):
    if not (0<tail[-1] and tail[0]<DEN and all(a>b for a,b in zip(tail,tail[1:]))):return 'outside_strict_sorted_positive_face'
    try:engine.projections(point(tail),distinct=True)
    except ValueError as e:return str(e)
    return None
def main():
    began=time.monotonic();ctx.prec=256;engine=GreenKernel()
    qualification=json.loads((P/'qualification.json').read_text());assert qualification['status']=='new3504_kernel_qualification_pass'
    inputs_path=P/'candidate_inputs.json'
    if inputs_path.exists():inputs=json.loads(inputs_path.read_text())
    else:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',UserWarning)
            uniforms=qmc.Sobol(d=4,scramble=True,seed=SEED).random(40)
        tails=[sorted(map(int,np.floor(DEN*row)),reverse=True) for row in uniforms]
        inputs={'seed':SEED,'denominator':DEN,'sobol_candidates':40,'tails':tails,'numpy':np.__version__,'scipy':scipy.__version__}
        save(inputs_path,inputs)
    assert inputs['seed']==SEED and inputs['denominator']==DEN and len(inputs['tails'])==40
    sources={str(p):sha(p) for p in [Path(__file__),P/'kernel.py',P/'registration.md',P/'qualification.json',inputs_path,TABLE]}
    fingerprint=hashlib.sha256(canonical(sources)).hexdigest();cp=CACHE/'batch_checkpoint.json'
    state=json.loads(cp.read_text()) if cp.exists() else {'fingerprint':fingerprint,'phase':'sobol','next_candidate':0,'records':[],
          'sobol_decisions':[],'local_decisions':[],'updates':[],'local_attempt':0,'best_record':None,'group_base':None,'group_records':[],
          'runtime_segments':[]}
    assert state['fingerprint']==fingerprint
    def checkpoint():save(cp,state)
    def stop_if_due():
        if time.monotonic()-began>120:
            state['runtime_segments'].append({'pid':os.getpid(),'seconds':time.monotonic()-began,'phase':state['phase']});checkpoint()
            print(json.dumps({'status':'checkpoint','pid':os.getpid(),'phase':state['phase'],'records':len(state['records'])}));return True
        return False
    def best(indices):return min(indices,key=lambda i:(F(state['records'][i]['G_upper_dyadic']),state['records'][i]['point']))
    if state['phase']=='sobol':
        for i in range(state['next_candidate'],40):
            tail=inputs['tails'][i];reason=validate(tail,engine);decision={'candidate':i,'tail_numerators':tail,'rejection':reason}
            if reason is None:
                record=engine.evaluate(point(tail),distinct=True);record.update(origin='sobol',candidate=i,tail_numerators=tail)
                decision['record_index']=len(state['records']);state['records'].append(record)
            state['sobol_decisions'].append(decision);state['next_candidate']=i+1;checkpoint()
            if len(state['records'])==24:
                state.update(phase='local',best_record=best(range(24)),sobol_best_record=best(range(24)));checkpoint();break
            if stop_if_due():return
        if state['phase']=='sobol':
            state['phase']='insufficient_24_valid_within_40';checkpoint()
    if state['phase']=='local':
        steps=[DEN//d for d in (8,16,32,64)]
        for attempt in range(state['local_attempt'],32):
            group=attempt//2;coordinate=group%4;step=steps[group//4];sign=1 if attempt%2==0 else -1
            if attempt%2==0:
                state['group_base']=list(state['records'][state['best_record']]['tail_numerators']);state['group_records']=[state['best_record']]
            tail=state['group_base'].copy();tail[coordinate]+=sign*step;tail.sort(reverse=True)
            reason=validate(tail,engine);decision={'attempt':attempt,'base_tail':state['group_base'],'coordinate':coordinate,'step':sign*step,
                                                  'tail_numerators':tail,'rejection':reason}
            if reason is None:
                found=next((i for i,r in enumerate(state['records']) if r['tail_numerators']==tail),None)
                if found is None:
                    record=engine.evaluate(point(tail),distinct=True);record.update(origin='local',attempt=attempt,tail_numerators=tail)
                    found=len(state['records']);state['records'].append(record);decision['already_evaluated']=False
                else:decision['already_evaluated']=True
                decision['record_index']=found;state['group_records'].append(found)
            state['local_decisions'].append(decision);state['local_attempt']=attempt+1
            if attempt%2:
                old=state['best_record'];new=best(state['group_records']);state['best_record']=new
                if new!=old:
                    state['updates'].append({'after_attempt':attempt,'from_record':old,'to_record':new,
                                             'strict_interval_improvement':F(state['records'][new]['G_upper_dyadic'])<F(state['records'][old]['G_lower_dyadic'])})
            checkpoint()
            if stop_if_due():return
        state['phase']='complete';checkpoint()
    state['runtime_segments'].append({'pid':os.getpid(),'seconds':time.monotonic()-began,'phase':state['phase']});checkpoint()
    result={'status':state['phase'],'precision_bits':256,'denominator':DEN,'seed':SEED,
            'sobol_candidates_examined':state['next_candidate'],'accepted_sobol_points':sum(r['origin']=='sobol' for r in state['records']),
            'local_proposals_attempted':state['local_attempt'],'new_local_evaluations':sum(r['origin']=='local' for r in state['records']),
            'sobol_decisions':state['sobol_decisions'],'local_decisions':state['local_decisions'],'records':state['records'],
            'updates':state['updates'],'sobol_best_record':state.get('sobol_best_record'),'best_record':state['best_record'],
            'selected_value':state['records'][state['best_record']] if state['best_record'] is not None else None,
            'new_global_g':False,'new_C5_bound':False,'source_hashes':sources,'checkpoint_fingerprint':fingerprint,
            'runtime_segments':state['runtime_segments'],'python':sys.executable}
    save(P/'batch_result.json',result)
    print(json.dumps({k:result[k] for k in ['status','sobol_candidates_examined','accepted_sobol_points','local_proposals_attempted','new_local_evaluations','best_record','selected_value','runtime_segments']},indent=2))
if __name__=='__main__':main()
