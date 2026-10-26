#!/usr/bin/env python3
"""118 one finite saved-evidence audit: no binary, compiler, kernel or Arb call."""
import sys
sys.dont_write_bytecode=True
import os,json,hashlib,difflib,tarfile,re,time
from pathlib import Path
from fractions import Fraction as Q
from datetime import datetime,timezone
OUT=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T')
AUTHOR=ROOT/'Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round118/band_engine_performance'
OLD=ROOT/'Math-Research-cube6-audit/cube-dim-6-bounds/upper/round016/band_engine'
START=time.monotonic();CPU=time.process_time();HASHES={}
def need(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):
    p=Path(p);v=hashlib.sha256(p.read_bytes()).hexdigest();old=HASHES.setdefault(str(p),v);need(old==v,'Changing evidence '+str(p));return v
def unique(pairs):
    d={}
    for k,v in pairs:need(k not in d,'Duplicate JSON key');d[k]=v
    return d
def read(p):sha(p);return json.loads(Path(p).read_text(),object_pairs_hook=unique,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON number '+x)))
def save(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
def resolve(p):return Path(p) if Path(p).is_absolute() else ROOT/p
def interval(d):
    need(set(d)=={'lower_mantissa','upper_mantissa','common_binary_exponent'},'Complete dyadic fields')
    for v in d.values():need(isinstance(v,str) and re.fullmatch(r'0|-?[1-9][0-9]*',v) and len(v)<=10000,'Canonical finite integer token')
    e=int(d['common_binary_exponent']);need(abs(e)<=10000,'Registered saved exponent range')
    t=Q(2)**e;l=int(d['lower_mantissa'])*t;u=int(d['upper_mantissa'])*t;need(l<=u,'Dyadic endpoint order');return l,u
def fields(d):
    need(len(d['gradient'])==6,'Six physical gradients')
    return {'band':d['band'],**{f'gradient[{i}]':d['gradient'][i] for i in range(6)},'p_a':d['p_a'],'p_T':d['p_T'],'Euler':d['Euler_xgrad_plus5B_plusTpT_minusapa']}
def scaled(iv,c):return (iv[0]*c,iv[1]*c) if c>=0 else (iv[1]*c,iv[0]*c)
def euler(d,command):
    x=[Q(int(z),int(command[7])) for z in command[8:14]];a=Q(int(command[3]),int(command[4]));T=Q(int(command[5]),int(command[6]))
    vals=[scaled(interval(d['band']),5),scaled(interval(d['p_T']),T),scaled(interval(d['p_a']),-a)]
    vals.extend(scaled(interval(g),q) for g,q in zip(d['gradient'],x));l=sum(v[0] for v in vals);u=sum(v[1] for v in vals)
    need(l<=0<=u,'Saved exact componentwise Euler consistency')
def main():
    reg=read(OUT/'registration.json');need(sha(Path(__file__))==reg['checker_sha256'],'Exact registered checker')
    save('started.json',dict(pid=os.getpid(),utc=datetime.now(timezone.utc).isoformat(),exclusive_create=True,registration_sha256=sha(OUT/'registration.json'),checker_sha256=sha(Path(__file__)),native_calls=0,compile_calls=0,finite_cases=4,dyadic_comparisons=40))
    need(sha(AUTHOR/'manifest_v2.json')==reg['author_manifest_sha256'] and sha(AUTHOR/'final_receipt_v2.json')==reg['author_final_receipt_sha256'],'Original final118 seal')
    manifest=read(AUTHOR/'manifest_v2.json')
    for row in manifest['artifacts']:
        p=AUTHOR/row['path'];need(p.resolve().is_relative_to(AUTHOR.resolve()),'Whitelisted artifact path')
        need(sha(p)==row['sha256'] and p.stat().st_size==row['bytes'],'Every original118 artifact')
    for row in manifest['external_pins_postrun']:need(sha(row['path'])==row['expected_sha256']==row['actual_sha256'] and row['matches'],'External runtime/fan/source pin')
    for name in ('source_manifest.json','candidate_v2_source_manifest.json'):
        for row in read(AUTHOR/name)['sources']:need(sha(resolve(row['path']))==row['sha256'],'Saved source provenance '+name)
    refs=read(AUTHOR/'saved_512_baselines_v2.json')['cases'];testreg=read(AUTHOR/'test_registration_v2.json');run=read(AUTHOR/'run_receipt_v2.json')
    final=read(AUTHOR/'final_receipt_v2.json');build=read(AUTHOR/'build_receipt_128.json');started=read(AUTHOR/'started_v2.json')
    need(final['manifest_sha256']==sha(AUTHOR/'manifest_v2.json') and final['run_receipt_sha256']==sha(AUTHOR/'run_receipt_v2.json') and final['report_sha256']==sha(AUTHOR/'REPORT_V2.md'),'Final receipt links')
    need(started['registration_sha256']==sha(AUTHOR/'test_registration_v2.json') and started['reference_sha256']==sha(AUTHOR/'saved_512_baselines_v2.json') and started['source_manifest_sha256']==sha(AUTHOR/'candidate_v2_source_manifest.json'),'Original exclusive-start bindings')
    old=(OLD/'band_engine.cpp').read_text();sha(OLD/'band_engine.cpp');new=(AUTHOR/'band_engine_precision.cpp').read_text()
    needle='static constexpr slong PREC=512;';replacement='#ifndef BAND_ENGINE_PRECISION_BITS\n#define BAND_ENGINE_PRECISION_BITS 512\n#endif\nstatic constexpr slong PREC=BAND_ENGINE_PRECISION_BITS;'
    need(old.count(needle)==1 and new==old.replace(needle,replacement),'Exact precision-only source patch')
    patch=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(OLD/'band_engine.cpp'),tofile='band_engine_precision.cpp'))
    need(patch==(AUTHOR/'precision.patch').read_text(),'Frozen patch exact match')
    need(testreg['patched_source_sha256']==sha(AUTHOR/'band_engine_precision.cpp') and testreg['patch_sha256']==sha(AUTHOR/'precision.patch'),'Registered source and patch')
    # Check the actual include bytes against the pinned archive, without extracting
    # or compiling: a header-archive hash alone does not identify the include tree.
    archive=OLD.parents[1]/'round010/arb_engine/headers_min.tar.gz';sha(archive);headers=[]
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar.getmembers():
            if not member.isfile():continue
            need(not Path(member.name).is_absolute() and '..' not in Path(member.name).parts,'Header archive path')
            p=OLD/'scratch'/member.name;need(p.resolve().is_relative_to((OLD/'scratch').resolve()),'Actual include path')
            content=tar.extractfile(member).read();need(p.read_bytes()==content,'Actual compiler include bytes equal sealed archive');headers.append({'path':str(p),'sha256':sha(p)})
    need(bool(headers),'Complete nonempty header archive')
    need(build['successful'] and build['bits']==128 and build['binary_sha256']==sha(AUTHOR/'band_engine_128'),'Original compiled binary hash')
    need(len(build['phases'])==2 and [p['phase'] for p in build['phases']]==['compile','post_link'] and all(p['exit_code']==0 and p['stderr']=='' for p in build['phases']),'Original successful compile/link receipts')
    compilecmd=build['phases'][0]['command'];need(compilecmd==testreg['builds'][0]['compile_command'] and '-DBAND_ENGINE_PRECISION_BITS=128' in compilecmd and not any('fast-math' in x for x in compilecmd),'Registered128 compile flags')
    expected=reg['case_ids'];need([r['case_id'] for r in refs]==[r['case_id'] for r in run['runs']]==[r['case_id'] for r in testreg['cases']]==expected,'Four exact ordered point-window inputs')
    old_sources=read(AUTHOR/'candidate_source_manifest.json')['sources'];comparisons=[]
    for baseline,record,planned in zip(refs,run['runs'],testreg['cases']):
        case=baseline['case_id'];command=record['command'];oldcommand=baseline['command']
        need(command==planned['native_command_128'] and command[1:]==oldcommand[1:] and len(command)==14,'Same exact rational input/window/fan')
        need([str(Q(int(x),int(command[7]))) for x in command[8:14]]==[str(Q(x)) for x in baseline['center']],'Physical point identity')
        need(Q(int(command[3]),int(command[4]))>0 and Q(int(command[3]),int(command[4]))<Q(int(command[5]),int(command[6])),'Positive-time band domain')
        need(record['exit_code']==0 and record['stderr']=='' and record['parse_error'] is None,'Recorded success, no partial/error output')
        output=Path(record['output_path']);need(output.resolve().is_relative_to((AUTHOR/'runs_v2').resolve()) and sha(output)==record['stdout_sha256'] and output.stat().st_size==record['stdout_bytes'],'Complete new stdout bytes')
        native=read(output);prior=baseline['native'];need(native['status']==prior['status']=='all_real_nodes_supported' and native['precision_bits']==128 and prior['precision_bits']==512,'Correct native status/precision metadata')
        for key in ('full_group_terms','computed_sign_pairs','confluent_terms','zero_node_terms','unique_node_values'):need(native[key]==prior[key],'Precision-independent exact topology '+key)
        need(native['full_group_terms']==967680 and native['computed_sign_pairs']==483840,'Complete logical/physical fan group')
        if case.startswith('C114.'):
            original=next(Path(row['path']) for row in old_sources if Path(row['path']).name=='native.json' and Path(row['path']).parent.name==case)
            for row in old_sources:
                if Path(row['path']).parent==original.parent:need(sha(row['path'])==row['sha256'],'Original high-window source bytes')
        else:original=ROOT/'Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round115/two_band_pilot/jobs/pilot1/native.json'
        original_native=read(original);original_receipt=read(original.parent/'native_receipt.json')
        need(sha(original)==baseline['reference_native_sha256']==planned['reference_native_sha256'] and fields(original_native)==fields(prior['interval_fields']),'Full saved512 intervals match raw original')
        need(original_receipt['native_exit_code']==0 and original_receipt['stderr_empty'] and original_receipt['native_sha256']==sha(original) and original_receipt['command'][1:]==command[1:],'Saved512 raw receipt/input identity')
        for field,value in fields(native).items():
            l,u=interval(value);a,b=interval(fields(prior['interval_fields'])[field]);need(l<=a<=b<=u,'Exact128 contains saved512 '+case+' '+field)
            comparisons.append({'case':case,'field':field,'contains_saved512':True,'outer_lower':str(l),'outer_upper':str(u),'saved512_lower':str(a),'saved512_upper':str(b)})
        need(interval(native['Euler_xgrad_plus5B_plusTpT_minusapa'])[0]<=0<=interval(native['Euler_xgrad_plus5B_plusTpT_minusapa'])[1],'Native Euler includes0')
        euler(native,command);euler(prior['interval_fields'],command)
    need(len(comparisons)==40 and not run['conditional_256_needed_cases'] and final['new_native_calls']==4 and final['fallback_256_calls']==0 and final['old_512_reruns']==0,'Finite original calibration counts')
    for p,v in list(HASHES.items()):need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v,'Postcheck changed source/output')
    result=dict(status='PASS118_EXACT_SAVED_PRECISION_PATCH_AND40_CONTAINMENTS',source_hashes=HASHES,headers_equal_archive=headers,
        exact_source_change='PREC compile-time parameterization only',comparisons=comparisons,all4_exact_topology_counts_match=True,
        static_Arb_enclosure_reliability_independent_of_precision=True,C016_theory_and_runtime_contract_still_required=True,
        measured_widths_and_throughput_not_universal=True,global_speedup_guarantee=False,full_gradient_independent_reference_gap_unchanged=True,
        ambiguous_finite_Euler_failure_not_silently_retryable=True,production117_changed=False,proposed_next_round_fallback_only=True,
        new_kernel_calls=0,new_compile_calls=0,new_native_binary_launches=0,children_created=0,source_writes=0,C6_improved=False)
    save('checked_saved.json',result)
    save('run_receipt.json',dict(status='PASS',pid=os.getpid(),utc=datetime.now(timezone.utc).isoformat(),attempts=1,
        wall_after_result_before_receipt=time.monotonic()-START,cpu_before_receipt=time.process_time()-CPU,
        runtime_limit=None,native_calls=0,compile_calls=0,result_sha256=sha(OUT/'checked_saved.json')))
    print(json.dumps({'status':result['status'],'containments':40,'headers_checked':len(headers),'native_calls':0,'compile_calls':0}))
if __name__=='__main__':
    try:main()
    except BaseException as error:
        save('failure_receipt.json',dict(status='STOP_NO_RETRY',pid=os.getpid(),error=repr(error),wall=time.monotonic()-START,cpu=time.process_time()-CPU,native_calls=0,compile_calls=0))
        raise
