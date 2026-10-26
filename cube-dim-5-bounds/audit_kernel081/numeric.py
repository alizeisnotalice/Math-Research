#!/usr/bin/env python3
"""Finite, resumable independent integer partial-fraction value audit."""
import argparse,gzip,json,os,time,resource
from pathlib import Path
from fractions import Fraction as F
from common import require,file_sha,digest
from arith import direct_values,dyadic_outward,P,LOG_TERMS,ATAN_TERMS

HERE=Path(__file__).resolve().parent


def write(path,obj):
    temp=Path(str(path)+'.tmp');temp.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n');temp.replace(path)


def load_inputs():
    reg=json.loads((HERE/'registration.json').read_text())
    require(file_sha(reg['manifest'])==reg['manifest_sha256'],'manifest hash mismatch')
    manifest=json.loads(Path(reg['manifest']).read_text())
    require(file_sha(reg['fan'])==reg['fan_sha256'],'true fan hash mismatch')
    fan=json.loads(gzip.decompress(Path(reg['fan']).read_bytes()))
    jobs=[]
    for batch in manifest['batches']:
        root=Path(batch['directory']);p=root/'point_manifest.json'
        require(file_sha(p)==batch['files']['point_manifest.json'],'point manifest hash mismatch')
        rows=json.loads(p.read_text())['point_checkpoints']
        require(len(rows)==batch['selected'],'point manifest count mismatch')
        for ordinal,row in enumerate(rows):jobs.append({'batch':batch['batch'],'ordinal':ordinal,**row})
    require(len(jobs)==reg['finite_new_points']==1318,'finite point count mismatch')
    require(len({tuple(j['point']) for j in jobs})==1318,'new points repeat')
    last=manifest['batches'][-1];overlay=Path(last['directory'])/'state.json'
    require(file_sha(overlay)==last['files']['state.json'],'last overlay mismatch')
    new_cache=json.loads(overlay.read_text())['new_cache']
    require(len(new_cache)==1318,'new cache count mismatch')
    return reg,manifest,fan,jobs,new_cache


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--count',type=int,default=64);args=parser.parse_args()
    require(1<=args.count<=128,'finite invocation cap is128 points')
    started=time.monotonic();pid=os.getpid()
    reg,manifest,fan,jobs,new_cache=load_inputs()
    codes={f:file_sha(HERE/f) for f in ['registration.json','common.py','arith.py','numeric.py']}
    path=HERE/'numeric_result.json'
    if path.exists():
        result=json.loads(path.read_text());require(result['code_locks']==codes,'successful-prefix code/input lock changed')
    else:
        result={'task':'audit_kernel081_numeric','status':'pending_values','code_locks':codes,'manifest_sha256':reg['manifest_sha256'],
                'fan_sha256':reg['fan_sha256'],'arithmetic':'integer outward only; direct rational poles; fixed-point atanh logs; Machin pi',
                'fixed_bits':P,'log_terms':LOG_TERMS,'atan_terms':ATAN_TERMS,'components':['H','physical_dx2','physical_dx3','physical_dx4','physical_dx5'],
                'shared_production_code':False,'shared_Arb_or_log_library':False,'shared_inputs':'true017fan/normalization/positive-eps resolvent analytic identity',
                'points':[],'invocations':[],'new_old_cache_evaluations':0,'MATLAB_calls':0}
    start=len(result['points']);end=min(start+args.count,len(jobs))
    print(json.dumps({'pid':pid,'numeric_start':start,'numeric_stop_exclusive':end,'total':len(jobs),'all5_components':True}),flush=True)
    for index in range(start,end):
        job=jobs[index]
        require(file_sha(job['path'])==job['sha256'],'raw point checkpoint hash mismatch')
        raw=json.loads(Path(job['path']).read_text());key=','.join(job['point']);record=new_cache[key]
        require(raw['stage']=='complete' and raw['point']==job['point'] and raw['record']==record,'raw point/record/new cache mismatch')
        point_started=time.monotonic();values,geometry=direct_values(job['point'],fan)
        saved=[]
        for independent,(lo,hi) in zip(values,record['values']):
            a,b=F(lo),F(hi);require(a<=independent[0]<=independent[1]<=b,'independent enclosure not contained in original ball')
            outward=dyadic_outward(independent);require(a<=F(outward[0])<=F(outward[1])<=b,'saved independent dyadic enclosure escaped original ball')
            saved.append(list(outward))
        require(str(raw['record']['attempts'][0]['minimum_integer_projection_gap'])==geometry['minimum_integer_projection_gap'],'independent projection gap differs')
        item={'index':index,'batch':job['batch'],'ordinal':job['ordinal'],'point':job['point'],'point_checkpoint_sha256':job['sha256'],
              'independent_values_512_dyadic':saved,'all5_strict_containment':True,'geometry':geometry,'seconds':time.monotonic()-point_started}
        result['points'].append(item);result['verified_count']=len(result['points'])
        result['status']='numeric_complete_all1318_H_plus4grad' if len(result['points'])==len(jobs) else 'pending_values'
        # Every completed point becomes an immutable successful prefix before next point.
        write(path,result)
        if (index+1)%16==0:print(json.dumps({'pid':pid,'verified_prefix':index+1}),flush=True)
    for p,h in reg['protected_old_locks'].items():require(file_sha(p)==h,'old protected file changed')
    invocation={'pid':pid,'start':start,'stop_exclusive':end,'count':end-start,'wall_seconds':time.monotonic()-started,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    result['invocations'].append(invocation);result['old_locks_unchanged']=True
    write(path,result)
    print(json.dumps({**invocation,'status':result['status'],'verified_count':len(result['points']),'numeric_result_sha256':file_sha(path)}),flush=True)


if __name__=='__main__':main()
