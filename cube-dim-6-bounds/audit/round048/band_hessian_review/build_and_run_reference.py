from pathlib import Path
import json,hashlib,time,subprocess,tarfile,os
P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
clock=json.loads((P/'input_manifest.json').read_text());S=Path(clock['source_directory'])
assert not(P/'reference_execution.json').exists()
receipt={'status':'PREPARING','pid':os.getpid(),'started_unix':time.time(),'reference_components':[[2,2],[5,6]],'full_21_independent_reference':False,'author_binary_executions':0}
try:
 for f,h in clock['input_sha256'].items():assert sha(f)==h,f
 for f,h in clock['source_files_sha256'].items():assert sha(f)==h,f
 m=json.loads((S/'manifest.json').read_text())
 for f,h in m['whitelist_sha256'].items():assert sha(S/f)==h,f
 for f,h in m['external_readonly_sha256'].items():assert sha(f)==h,f
 record=json.loads((S/'build_record.json').read_text());archive=Path(record['headers_archive'])
 target=P/'scratch';target.mkdir(exist_ok=True)
 with tarfile.open(archive) as t:t.extractall(target,filter='data')
 lib=Path(next(p for p in record['runtime_dylib_sha256'] if 'libflint.' in p));binary=target/'two_components'
 command=['/usr/bin/clang++','-O3','-std=c++20','-I',str(target/'include'),str(P/'reference.cpp'),str(lib),'-Wl,-rpath,'+str(lib.parent),'-o',str(binary)]
 b=subprocess.run(command,capture_output=True,text=True);receipt['compile']={'command':command,'exit_code':b.returncode,'stdout':b.stdout,'stderr':b.stderr};assert b.returncode==0
 patch=['/usr/bin/install_name_tool','-change','/DLC/flint/.dylibs/libflint.24.0.dylib','@rpath/libflint.24.0.dylib',str(binary)]
 q=subprocess.run(patch,capture_output=True,text=True);receipt['patch']={'command':patch,'exit_code':q.returncode,'stderr':q.stderr};assert q.returncode==0
 receipt.update(status='RUNNING_ONE_BATCH',source_sha256=sha(P/'reference.cpp'),binary_sha256=sha(binary))
 (P/'reference_execution.json').write_text(json.dumps(receipt,indent=2)+'\n')
 fan=next(p for p in m['external_readonly_sha256'] if p.endswith('fan_integer_input.txt'))
 remaining=clock['started_unix']+1200-time.time();assert remaining>60
 start=time.monotonic();run=subprocess.run([str(binary),fan],capture_output=True,text=True,timeout=min(600,remaining-30))
 receipt.update(exit_code=run.returncode,process_seconds=time.monotonic()-start,stderr=run.stderr)
 (P/'reference_stdout.txt').write_text(run.stdout)
 assert run.returncode==0,run.stderr
 result=json.loads(run.stdout);assert result['logical_terms_each']==967680 and result['precision_bits']>=768 and not result['full_21_independent_reference']
 (P/'reference_results.json').write_text(json.dumps(result,indent=2)+'\n')
 receipt.update(status='COMPLETE_TWO_REFERENCES',finished_unix=time.time(),audit_elapsed_seconds=time.time()-clock['started_unix'])
 print(json.dumps({k:receipt[k] for k in ['status','pid','process_seconds','audit_elapsed_seconds']}))
except BaseException as e:
 receipt.update(status='FAILED_STOP_NO_RETRY',error=repr(e),finished_unix=time.time());raise
finally:(P/'reference_execution.json').write_text(json.dumps(receipt,indent=2)+'\n')
