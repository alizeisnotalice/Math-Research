#!/usr/bin/env python3
"""One first/second full check, <=180s process-group deadline, six processes.

Static wrapper only until the parent authorizes execution. No retry, repair,
LP or pre-checker mathematical evaluation; mass gates use saved Fractions.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import monotonic,time

ROOT=Path(__file__).resolve().parent
SECONDS=180.
INPUTS={
    '../round25/candidate.json':'7081cfef6fb47e5d47c95d3740c7a6ba48d48bb6197928a84d8cb6340c4166d5',
    '../round25/structured_kernel.json':'00c50260612ec3e8b15e810af2fe62269c10a09db0a7cd136c624df8a90bc02f',
    '../round25/exact_cost.json':'0120e7a55961ae7858541fb95319a1e4454c3c4e953c603478e51e488b7f8b33',
    '../round25/candidate_saved_verification.json':'b3be829ebf4b20c7cacca750d67a297a23914a1627087c82e2b1c974962ea790',
    '../round25/verify_candidate_saved.py':'0ef08169ce9d25517c8d7d00e2567a89863c71790538cbc89419756b5d64f97b',
    '../round25/solve_result.json':'50bed7523f577382aee952be871ad248f3efb38a46441d4a50d7596c3592b239',
    '../round25/solve.run.json':'e5e76d45566ccad6555fab1899a58ac0eadf87ba9b7829293b55effa2d2f9541',
    '../round25/solve.py':'e9074068c755b7635f7fcdfb2acf262e6a1d48e49143b4490e5d581978c753cd'}
FIRST={
    '../round17/fixed_checker/cert_run.py':'079906db354febdd57d537f0d747428bc90bd37fa88297449cd1d46c7446310d',
    '../round17/fixed_checker/hlcert.py':'4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05'}
SECOND={
    '../upstream-v4/checker/independent/hl_check.py':'0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8',
    '../upstream-v4/checker/independent/hl_core.py':'c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'}


def digest(path):return sha256(path.read_bytes()).hexdigest()


def read(relative):return json.loads((ROOT/relative).read_text())


def require(condition,message):
    if not condition:raise ValueError(message)


def dump(path,value):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n');tmp.replace(path)


def main():
    require(len(sys.argv)==2 and sys.argv[1] in ('first','second'),'Choose exactly first or second.')
    stage=sys.argv[1];meta=ROOT/(stage+'.run.json');log=ROOT/(stage+'.log');report=ROOT/(stage+'_report.json')
    if any(p.exists() for p in (meta,log,report)):
        raise SystemExit('Existing stage evidence: no repeated checker.')
    start=monotonic();deadline=start+SECONDS;wrapper_before=digest(Path(__file__))
    record={'status':'started','stage':stage,'budget_seconds':180,'started_unix':time(),
        'processes':6,'claim':'359087/100000','claim_decimal':'3.59087','mathematical_runs_started':0,
        'wrapper_sha256_before':wrapper_before,'input_hashes_before':{},'source_hashes_before':{}}
    with meta.open('x') as stream:stream.write(json.dumps(record,indent=2)+'\n')
    proc=None
    def guard():
        if monotonic()>=deadline:raise TimeoutError('180s total process deadline reached.')
    def check_pins(pins):
        result={}
        for relative,wanted in pins.items():
            guard();actual=digest(ROOT/relative)
            require(actual==wanted,'Pinned input/source mismatch: '+relative);result[relative]=actual
        return result
    source_pins=FIRST if stage=='first' else SECOND
    try:
        record['input_hashes_before']=check_pins(INPUTS)
        record['source_hashes_before']=check_pins(source_pins)
        solve=read('../round25/solve_result.json');solve_run=read('../round25/solve.run.json')
        candidate=read('../round25/candidate.json');structured=read('../round25/structured_kernel.json')
        cost=read('../round25/exact_cost.json');audit=read('../round25/candidate_saved_verification.json')
        require(solve['status']=='one_finite_lp_and_exact_local_check_closed' and solve['candidate_exported'] is True,
            'Finite discovery did not export a candidate.')
        require(solve_run['status']=='terminal' and solve_run['exit_code']==0 and
            solve_run['source_sha256_before']==solve_run['source_sha256_after']==INPUTS['../round25/solve.py'],
            'Discovery run/source binding mismatch.')
        require(solve['candidate_sha256']==INPUTS['../round25/candidate.json'] and
            solve['exact_cost_sha256']==INPUTS['../round25/exact_cost.json'],'Discovery candidate/cost hash mismatch.')
        require(audit['status']=='passed_saved_candidate_arithmetic_audit' and audit['worst_rows_reconstructed']==1 and
            audit['candidate_sha256']==INPUTS['../round25/candidate.json'] and audit['structured_sha256']==INPUTS['../round25/structured_kernel.json'] and
            audit['exact_cost_sha256']==INPUTS['../round25/exact_cost.json'] and audit['verifier_sha256']==INPUTS['../round25/verify_candidate_saved.py'],
            'Independent saved-value candidate audit missing.')
        require(set(candidate)=={'alpha','R','h','N','c','a'} and all(candidate[k]==structured[k] for k in candidate),
            'Plain/structured schema or data mismatch.')
        require((Q(candidate['alpha']),Q(candidate['R']),Q(candidate['h']),candidate['N'])==(Q(47,40),Q(33,20),Q(1,32),80) and
            Q(structured['beta'])==Q(candidate['alpha']),'Wrong fixed kernel family.')
        require(cost['certificate_sha256']==INPUTS['../round25/structured_kernel.json'] and cost['improves_baseline'] is True and
            Q(cost['upper'])<Q(cost['baseline_lower']) and Q(cost['upper'])<Q(359087,100000),
            'Exact cost/claim gate failed; do not change claim.')
        relative_script='../round17/fixed_checker/cert_run.py' if stage=='first' else '../upstream-v4/checker/independent/hl_check.py'
        args=['--Q','8','--maxdepth','7','--Kdepth','6','--procs','6','--only','MKA'] if stage=='first' else \
            ['--claim','3.59087','--procs','6','--maxdepth','8','--kdepth','10','--only','c1,c2,c3,c4']
        script=(ROOT/relative_script).resolve();candidate_path=(ROOT/'../round25/candidate.json').resolve()
        command=[sys.executable,'-I',str(script),str(candidate_path),*args,'--out',str(report)]
        record['command_relative']=['python','-I',relative_script,'../round25/candidate.json',*args,'--out',report.name]
        record['working_directory_relative']=str(Path(relative_script).parent)
        record['input_gate_passed']=True;record['elapsed_before_launch_seconds']=monotonic()-start;dump(meta,record);guard()
        with log.open('x') as stream:
            proc=subprocess.Popen(command,cwd=script.parent,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            record.update(pid=proc.pid,owned_process_group=proc.pid,mathematical_runs_started=1)
            try:require(os.getpgid(proc.pid)==proc.pid,'Checker did not own its new process group.')
            except ProcessLookupError:pass  # A very early exited child is still harvested below.
            dump(meta,record)
            try:
                code=proc.wait(timeout=max(.001,deadline-monotonic()))
                record.update(status='terminal',exit_code=code)
            except subprocess.TimeoutExpired:
                try:os.killpg(proc.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                proc.wait()
                record.update(status='timed_out',exit_code=proc.returncode)
    except Exception as exc:
        record.update(status='failed_gate_or_wrapper',reason=f'{type(exc).__name__}: {exc}')
    finally:
        if proc is not None:
            try:os.killpg(proc.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            if proc.poll() is None:proc.wait()
            record.setdefault('exit_code',proc.returncode)
        # Final hash/metadata work is not another mathematical call. Preserve it after a timeout.
        record['input_hashes_after']={relative:digest(ROOT/relative) for relative in INPUTS if (ROOT/relative).is_file()}
        record['source_hashes_after']={relative:digest(ROOT/relative) for relative in source_pins if (ROOT/relative).is_file()}
        record['inputs_unchanged']=record['input_hashes_before']==record['input_hashes_after']
        record['sources_unchanged']=record['source_hashes_before']==record['source_hashes_after']
        record['wrapper_sha256_after']=digest(Path(__file__))
        record['wrapper_unchanged']=record['wrapper_sha256_after']==wrapper_before
        if log.is_file():record.update(log_sha256=digest(log),log_bytes=log.stat().st_size)
        record['report_available']=report.is_file()
        if report.is_file():record.update(report_sha256=digest(report),report_bytes=report.stat().st_size)
        record['elapsed_seconds']=monotonic()-start
        if record['status']=='terminal' and not (record['inputs_unchanged'] and record['sources_unchanged'] and record['wrapper_unchanged']):
            record['status']='provenance_changed_invalid'
        dump(meta,record)
    print(json.dumps(record,ensure_ascii=False))
    return 0 if record['status']=='terminal' and record['exit_code']==0 else 1


if __name__=='__main__':sys.exit(main())
