#!/usr/bin/env python3
"""175 one exclusive actual dispatch; never retries mathematical computation."""
from pathlib import Path
import argparse,datetime,os,subprocess,sys,time
from common import sha,read,put
HERE=Path(__file__).resolve().parent

def snapshot(plan):
    src=read(HERE/'sources.json');pre=read(plan['preflight']['path'])
    refs=dict(src['source_locks'])
    for i,r in enumerate(pre['shards']+[pre['topology']]):refs[f'prepared175_{i}']=r
    p166=read(src['source_locks']['metadata166']['path'])
    for i,r in enumerate(p166['shards']+[p166['topology']]):refs[f'prepared166_{i}']=r
    root=Path(src['source_locks']['completion166']['path']).parent
    for name,r in read(src['source_locks']['outputs166']['path'])['files'].items():refs['166actual/'+name]={'path':str(root/name),**r}
    for name,h in plan['code_and_source_hashes'].items():refs['175code/'+name]={'path':str(HERE/name),'sha256':h}
    refs['175plan']={'path':str(HERE/'plan.json'),'sha256':sha(HERE/'plan.json')}
    result={}
    for k,r in refs.items():
        p=Path(r['path']);value={'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}
        assert value['sha256']==r['sha256']
        if 'bytes' in r:assert value['bytes']==r['bytes']
        result[k]=value
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--plan-sha',required=True);ap.add_argument('--base',type=Path,required=True);a=ap.parse_args()
    assert sha(HERE/'plan.json')==a.plan_sha;plan=read(HERE/'plan.json');assert sha(__file__)==plan['code_and_source_hashes']['dispatch.py']
    dest=a.base/'dispatch001';dest.mkdir(exist_ok=False);out=a.base/'run001';assert not out.exists()
    put(dest/'registration.json',{'task':'175_single_dispatch','pid':os.getpid(),'plan_sha256':a.plan_sha,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'output':str(out),'attempts':1,'wall_cap':None,'authorization_sha256':plan['authorization_sha256']})
    before=snapshot(plan);put(dest/'sources_before.json',before);start=time.monotonic()
    log=(dest/'stdout.log').open('xb');err=(dest/'stderr.log').open('xb')
    p=subprocess.Popen([sys.executable,'-B',str(HERE/'consumer.py'),'run','--output',str(out),'--plan-sha',a.plan_sha],stdout=log,stderr=err)
    put(dest/'master_pid.json',{'pid':p.pid,'controller_pid':os.getpid()});rc=p.wait();log.close();err.close()
    put(dest/'master_exit.json',{'pid':p.pid,'exitcode':rc,'seconds':time.monotonic()-start})
    if rc!=0:put(dest/'stopped.json',{'reason':'nonzero master','exitcode':rc,'prefix_preserved':True,'automatic_retry':False});raise SystemExit(rc)
    after=snapshot(plan);put(dest/'sources_after.json',after);assert before==after
    outputs={}
    for q in sorted(out.rglob('*.json')):
        assert q.name!='fault.json';outputs[str(q.relative_to(out))]={'bytes':q.stat().st_size,'sha256':sha(q)}
    put(dest/'output_manifest.json',{'files':outputs,'bytes':sum(x['bytes'] for x in outputs.values())})
    put(dest/'completion.json',{'task':'upper_secant175_actual001','controller_pid':os.getpid(),'master_pid':p.pid,'master_exitcode':rc,
        'seconds':time.monotonic()-start,'sources_before_sha256':sha(dest/'sources_before.json'),'sources_after_sha256':sha(dest/'sources_after.json'),
        'all_inputs_before_after_unchanged':True,'completion_sha256':sha(out/'completion.json'),'output_manifest_sha256':sha(dest/'output_manifest.json'),
        'worker_exitcodes':read(out/'worker_exitcodes.json'),'global_Cn_claim':False})
    print(read(dest/'completion.json'))
