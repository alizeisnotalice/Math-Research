"""Run exactly one released remaining-case reference, never retry or resume."""
import signal
import time
START=time.monotonic()
def timeout(signum,frame):raise TimeoutError('180-second full reference runner limit')
signal.signal(signal.SIGALRM,timeout)
signal.setitimer(signal.ITIMER_REAL,180)

from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import os
import subprocess
import sys

P=Path(__file__).resolve().parent
child=None
receipt=None
folder=None
def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save_new(path,obj):
    with Path(path).open('x') as f:json.dump(obj,f,indent=2);f.write('\n')
def write_receipt():
    if receipt is None:return
    receipt['runner_wall_seconds']=time.monotonic()-START
    temp=folder/'receipt.json.tmp';temp.write_text(json.dumps(receipt,indent=2)+'\n');temp.replace(folder/'receipt.json')
def dyadic(raw):
    assert set(raw)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
    scale=F(2)**int(raw['common_binary_exponent'])
    lo,hi=int(raw['lower_mantissa'])*scale,int(raw['upper_mantissa'])*scale
    assert lo<=hi
    return lo,hi

try:
    if len(sys.argv)!=3 or sys.argv[1] not in ('low-primary','S3-primary'):
        raise ValueError('Usage: run_reference_once.py CASE MAIN_FINAL_RELEASE.json')
    case=sys.argv[1];release_path=Path(sys.argv[2]).resolve();release=read(release_path)
    assert release['schema']=='F061-main-full-gradient-release-v1'
    assert release['authorized'] is True and release['released_by']=='main'
    assert release['cases']==['low-primary','S3-primary'] and release['blocks_per_case']==1
    assert release['program_wall_timeout_seconds']==180 and release['max_remaining_native_processes']==2
    assert release['preparation_manifest_sha256']==sha(P/'preparation_manifest.json')
    seal=read(P/'preparation_manifest.json')
    assert seal['status']=='SEALED_PREPARATION_PASS_AWAITING_FINAL_RELEASE'
    for name,digest in seal['artifact_sha256'].items():assert sha(P/name)==digest,name
    for path,digest in read(P/'input_manifest.json')['sha256'].items():assert sha(path)==digest,path
    author=Path(seal['author_directory'])
    assert sha(author/'pilot_manifest.json')==seal['author_pilot_manifest_sha256']
    for name,digest in read(author/'pilot_manifest.json')['artifact_sha256'].items():assert sha(author/name)==digest,name
    for path,digest in read(author/'input_manifest.json')['source_sha256'].items():assert sha(path)==digest,path
    assert read(P/'preparation_results.json')['status']=='PASS_NINE_PREPARATION_CONFIGURATIONS_NO_FULL_REFERENCE_RELEASE'
    binary=P/'scratch'/'physical_pf_reference'
    assert sha(binary)==seal['binary_sha256']
    fan=read(author/'registration.json')['fan_path']
    benchmark_path=P/f'benchmark_{case}.json';benchmark=read(benchmark_path)
    assert sha(benchmark_path)==seal['benchmark_sha256'][case]
    assert benchmark['case_id']==case and benchmark['mode']=='benchmark' and benchmark['computed_tuples']==128
    root=P/'full_reference';root.mkdir(exist_ok=True)
    expected=[] if case=='low-primary' else ['low-primary']
    assert sorted(x.name for x in root.iterdir())==expected,'Unique ordered batch required'
    if expected:
        prior=read(root/'low-primary'/'receipt.json')
        assert prior['status']=='COMPLETE_SIX_GRADIENT_COMPONENTS_MATCH_NATIVE'
        assert prior['release_sha256']==sha(release_path)
        assert prior['comparison_sha256']==sha(root/'low-primary'/'comparison.json')
    folder=root/case;folder.mkdir()
    command=[str(binary),'remaining',case,fan,'0','1']
    receipt={'status':'REGISTERED_NOT_LAUNCHED','case_id':case,'command':command,
             'release_sha256':sha(release_path),'release_path':str(release_path),
             'preparation_manifest_sha256':sha(P/'preparation_manifest.json'),
             'source_sha256':sha(P/'physical_pf_reference.cpp'),'binary_sha256':sha(binary),
             'fan_sha256':sha(fan),'benchmark_sha256':sha(benchmark_path),
             'program_wall_timeout_seconds':180,'new_remaining_native_calls':0,
             'complete_reference_fields':0,'automatic_retries':0,'author_native_calls':0,
             'global_g_certificate':False,'new_C6_bound':False}
    write_receipt()
    stdout_path=folder/'native_stdout.json';stderr_path=folder/'native_stderr.txt'
    with stdout_path.open('xb') as out,stderr_path.open('xb') as err:
        child=subprocess.Popen(command,cwd=P,stdout=out,stderr=err,start_new_session=True)
        receipt.update(status='RUNNING_REMAINING_TUPLES',owned_child_pid=child.pid,
                       owned_process_group_id=child.pid,new_remaining_native_calls=1);write_receipt()
        code=child.wait(timeout=max(0.001,180-(time.monotonic()-START)))
    receipt.update(native_exit_code=code,native_stdout_sha256=sha(stdout_path),native_stderr_sha256=sha(stderr_path));write_receipt()
    assert code==0 and stderr_path.stat().st_size==0,'Native failure, no retry'
    remaining=read(stdout_path)
    assert remaining['case_id']==case and remaining['mode']=='remaining'
    assert remaining['precision_bits']==768 and remaining['block']==0 and remaining['blocks']==1
    assert remaining['computed_tuples']==967552 and remaining['physical_fields']==6
    assert remaining['tuple_field_contributions']==5805312 and remaining['skip_benchmark_in_remaining_mode'] is True
    native=read(author/'production'/case/'native.json')
    assert benchmark['zero_base_tuples']+remaining['zero_base_tuples']==native['zero_node_terms']
    assert benchmark['repeated_base_tuples']+remaining['repeated_base_tuples']==native['confluent_terms']
    rows=[]
    for k,(b,r,n) in enumerate(zip(benchmark['gradient_partial_normalized'],remaining['gradient_partial_normalized'],native['gradient'])):
        bl,bu=dyadic(b);rl,ru=dyadic(r);nl,nu=dyadic(n);lo,hi=bl+rl,bu+ru
        rows.append({'physical_field':k+1,'reference_lower_exact':str(lo),'reference_upper_exact':str(hi),
                     'native_lower_exact':str(nl),'native_upper_exact':str(nu),
                     'native_contains_complete_reference':nl<=lo<=hi<=nu})
    comparison={'status':'COMPLETE_REFERENCE_COMPARISON','case_id':case,'window':{'kind':'B','a':'1/250','T':'1/2'},
                'benchmark_tuples_reused':128,'remaining_tuples_new':967552,'complete_tuples_each_field':967680,
                'physical_fields':6,'complete_tuple_field_contributions':5806080,
                'zero_base_tuples':native['zero_node_terms'],'repeated_base_tuples':native['confluent_terms'],
                'source_native_sha256':sha(author/'production'/case/'native.json'),'fields':rows,
                'all_six_native_intervals_contain_reference':all(r['native_contains_complete_reference'] for r in rows),
                'Euler_used_as_gradient_proof':False,'global_g_certificate':False,'new_C6_bound':False}
    save_new(folder/'comparison.json',comparison)
    assert len(rows)==6 and comparison['all_six_native_intervals_contain_reference'],'Complete component comparison failed; stop'
    for name,digest in seal['artifact_sha256'].items():assert sha(P/name)==digest,name
    receipt.update(status='COMPLETE_SIX_GRADIENT_COMPONENTS_MATCH_NATIVE',complete_reference_fields=6,
                   new_remaining_tuples=967552,benchmark_tuples_reused=128,
                   complete_tuples_each_field=967680,native_kernel_seconds=remaining['kernel_seconds'],
                   comparison_sha256=sha(folder/'comparison.json'))
    print(json.dumps({'case_id':case,'status':receipt['status'],'complete_reference_fields':6,
                      'new_remaining_tuples':967552,'benchmark_tuples_reused':128,
                      'native_kernel_seconds':remaining['kernel_seconds'],'global_g_certificate':False}))
except BaseException as error:
    signal.setitimer(signal.ITIMER_REAL,0)
    if child is not None and child.poll() is None:
        try:os.killpg(child.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        child.wait()
        if receipt is not None:receipt['killed_and_waited_only_owned_process_group']=True
    if receipt is not None:
        receipt.update(status='FAILED_OR_TIMEOUT_STOP_NO_RETRY',error=repr(error),complete_reference_fields=0)
        for name in ('native_stdout.json','native_stderr.txt','comparison.json'):
            path=folder/name
            if path.exists():receipt[name+'_retained_sha256']=sha(path)
    raise
finally:
    signal.setitimer(signal.ITIMER_REAL,0)
    write_receipt()
