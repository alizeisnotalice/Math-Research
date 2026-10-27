#!/usr/bin/env python3
"""One authorized179A byte encode and independent cold decode;no retries."""
from pathlib import Path
from datetime import datetime,timezone
import ast,hashlib,json,os,subprocess,sys,time

HERE=Path(__file__).resolve().parent
PLAN='c037e917cf87506e0a7e3e0b0969d9854704ced4e47e36c2de201368ee5e5504'
PREP='ff3caec1a3f90f39954bc573b891c40790061daa5bd5ce3d5bdedad59a37b4af'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())

if __name__=='__main__':
    assert sha(HERE/'plan.json')==PLAN and sha(HERE/'manifest.json')==PREP
    for n,x in read(HERE/'manifest.json')['files'].items():assert sha(HERE/n)==x['sha256']
    for n in ['encode179.py','decode179.py','run179.py']:ast.parse((HERE/n).read_text())
    base=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_combined_archive179');base.mkdir(exist_ok=True)
    out=base/'run001';decoded=base/'decode001';dispatch=base/'dispatch001';prereg=base/'preregistration001.json'
    assert not out.exists() and not decoded.exists() and not dispatch.exists() and not prereg.exists();dispatch.mkdir()
    repo=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5');plan=read(HERE/'plan.json')
    put(prereg,{'task':'179A_one_byte_pack_and_independent_decode','authorization':'Primary explicitly approved profileA plan and fixed32projections/302literals/1173targets;no numerical algorithms',
        'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),'plan_sha256':PLAN,'preparation_manifest_sha256':PREP,'profile':'A',
        'encode_sha256':sha(HERE/'encode179.py'),'decode_sha256':sha(HERE/'decode179.py'),'controller_sha256':sha(__file__),
        'selected_sources':{'path':str(HERE/'selection_final.json'),'sha256':plan['selection_final_sha256'],'logical_paths':1525,'literal_blobs':302,'original_math_targets':1173},
        'repository_refs':{'path':str(HERE/'repository_refs.json'),'sha256':plan['repository_refs_sha256'],'count':95},
        'field_recipes_sha256':plan['field_recipes_sha256'],'projection_plan_sha256':plan['projection_plan_sha256'],'projection_payloads':32,
        'pack_group_cap_bytes':33554432,'physical_pack_hard_gate_bytes':50000000,'exclusive_encode_output':str(out),'exclusive_decode_output':str(decoded),
        'no_math_no_H_binary_decode':True,'no_automatic_retry_fallback_or_extra_sources':True})
    psha=sha(prereg);start=time.monotonic()
    commands=[('encode',[sys.executable,'-B',str(HERE/'encode179.py'),'--output',str(out),'--repository',str(repo),'--preregistration',str(prereg),'--preregistration-sha',psha])]
    results=[]
    try:
        for label,cmd in commands:
            with (dispatch/(label+'.stdout')).open('xb') as sout,(dispatch/(label+'.stderr')).open('xb') as serr:
                p=subprocess.Popen(cmd,stdout=sout,stderr=serr);put(dispatch/(label+'_pid.json'),{'pid':p.pid,'utc':datetime.now(timezone.utc).isoformat(),'command':cmd});print(label.upper()+'_PID',p.pid,flush=True);code=p.wait()
            r={'stage':label,'pid':p.pid,'exitcode':code,'stdout_sha256':sha(dispatch/(label+'.stdout')),'stderr_sha256':sha(dispatch/(label+'.stderr'))};results.append(r);put(dispatch/(label+'_exit.json'),r)
            assert code==0,(label,code)
        manifest=out/'archive_manifest.json';msha=sha(manifest)
        cmd=[sys.executable,'-B',str(HERE/'decode179.py'),'--archive',str(out),'--manifest-sha',msha,'--repository',str(repo),'--output',str(decoded)]
        with (dispatch/'decode.stdout').open('xb') as sout,(dispatch/'decode.stderr').open('xb') as serr:
            p=subprocess.Popen(cmd,stdout=sout,stderr=serr);put(dispatch/'decode_pid.json',{'pid':p.pid,'utc':datetime.now(timezone.utc).isoformat(),'command':cmd});print('DECODE_PID',p.pid,flush=True);code=p.wait()
        r={'stage':'independent_decode','pid':p.pid,'exitcode':code,'stdout_sha256':sha(dispatch/'decode.stdout'),'stderr_sha256':sha(dispatch/'decode.stderr')};results.append(r);put(dispatch/'decode_exit.json',r);assert code==0,code
        assert sha(prereg)==psha and sha(HERE/'plan.json')==PLAN and sha(HERE/'manifest.json')==PREP
        reg=read(prereg);assert sha(HERE/'encode179.py')==reg['encode_sha256'] and sha(HERE/'decode179.py')==reg['decode_sha256'] and sha(__file__)==reg['controller_sha256']
        put(dispatch/'completion.json',{'status':'one_encode_independent_decode_both_exit0','controller_pid':os.getpid(),'seconds':time.monotonic()-start,'stages':results,
            'preregistration_sha256':psha,'archive_manifest_sha256':msha,'encode_receipt_sha256':sha(out/'encode_receipt.json'),'decode_receipt_sha256':sha(decoded/'decode_receipt.json'),
            'restoration_checks_sha256':sha(decoded/'restoration_checks.json'),'mathematical_algorithm_calls':0,'H_binary_decodes':0,'attempts':1})
        print(json.dumps(read(dispatch/'completion.json'),sort_keys=True),flush=True)
    except BaseException as ex:
        put(dispatch/'fault.json',{'controller_pid':os.getpid(),'exception':repr(ex),'completed_stages':results,'prefix_preserved':True,'action':'stop,no retry,profileB or new sources'})
        raise
