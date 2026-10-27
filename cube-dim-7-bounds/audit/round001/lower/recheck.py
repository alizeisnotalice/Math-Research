#!/usr/bin/env python3
"""Portable, read-only exact recheck of the two sealed seven-dimensional seeds.

Default candidate: this project's lower/round001. All output goes to stdout;
no result, lock, receipt or bytecode file is written. --reference optionally
authenticates the frozen comparison against the complete original R075 file.
The historical check.py is imported only for its independent arithmetic and
full labelled-box union functions; its path-bound main() is never called.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, importlib.util, json, platform, sys, time

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
SEALED_CHECK_SHA256='211db5dbf67dabd619d753a0e8e0065be216e01480cff0f81b06567964c22e70'

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(condition,message):
    if not condition: raise ValueError(message)

class Parser(argparse.ArgumentParser):
    def error(self,message): raise ValueError(message)

def recheck(candidate_root,reference=None):
    require(not sys.flags.optimize,'Run without -O: sealed arithmetic uses assertions.')
    candidate=Path(candidate_root).expanduser().resolve()
    lock=json.loads((HERE/'input_lock.json').read_text())
    manifest=json.loads((HERE/'manifest.json').read_text())
    audit_before={name:sha(HERE/name) for name in manifest['files']}
    require(audit_before==manifest['files'],'Historical sealed audit file hash mismatch.')
    require(sha(HERE/'check.py')==SEALED_CHECK_SHA256,'Independent arithmetic hash mismatch.')
    candidate_before={name:sha(candidate/name) for name in lock['candidate_files']}
    require(candidate_before==lock['candidate_files'],'Candidate differs from the originally audited final lock.')
    # No author module is imported. Suppress bytecode caching to remain read-only.
    spec=importlib.util.spec_from_file_location('sealed_independent_cube7_audit',HERE/'check.py')
    audit=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    frozen=json.loads((candidate/'frozen_internal_reference.json').read_text())
    interval=tuple(map(Q,frozen['interval']))
    require(len(interval)==2 and interval[0]<=interval[1],'Invalid frozen reference interval.')
    reference_path=None; reference_hash=None
    if reference is None:
        comparison={'mode':'frozen_interval_only',
                    'original_reference_file_verified':False,
                    'statement':'Comparison uses the complete frozen rational interval; the original R075 file was not accessed in this run.'}
    else:
        reference_path=Path(reference).expanduser().resolve()
        reference_hash=sha(reference_path)
        require(reference_hash==lock['baseline_sha256']==frozen['sha256'],'Original R075 reference hash mismatch.')
        original=json.loads(reference_path.read_text())['candidate']['exact_interval']
        require(tuple(map(Q,original))==interval,'Original and frozen complete intervals differ.')
        comparison={'mode':'original_reference_hash_and_interval_verified',
                    'original_reference_file_verified':True,'reference_path':str(reference_path),
                    'reference_sha256':reference_hash}
    output={'status':'PASS exact portable read-only recheck','candidate_root':str(candidate),
            'comparison':comparison,'reference_interval':list(map(str,interval)),
            'algorithm':'independent full labelled raw boxes, fixed-coordinate event sweep; no author arithmetic imports',
            'historical_check_sha256':SEALED_CHECK_SHA256,'seeds':[]}
    for power in (4,5):
        t=time.time(); cpu=time.process_time()
        claim=json.loads((candidate/f'seed{power}_results.json').read_text())
        inner,outer,checks,minwidth,maxspacing,raw_count=audit.construct(power,claim)
        budget=[0]; volumes=[]; states=[]
        for direction,boxes in (('inner',inner),('outer',outer)):
            engine=audit.RawUnion(budget)
            volume=Q(128*engine.volume(boxes),audit.S**7)
            field='union_lower' if direction=='inner' else 'union_upper'
            require(volume==Q(claim[field]),f'Seed {power} {direction} union mismatch.')
            volumes.append(volume); states.append(len(engine.memo))
        require(volumes[0]<=volumes[1],'Inward/outward direction mismatch.')
        eps=Q(1,10**25); theta=1-Q(1,10**18); K=10**19
        require(minwidth>2*eps and maxspacing/2<2,'Smearing geometric budget failed.')
        loss=128*raw_count*14*eps*2**6
        finite=theta*(volumes[0]-loss)*Q(2*K+1,2*K+5)**7
        require(loss==Q(claim['finite_L1']['union_shrink_loss_upper']),'Smearing loss mismatch.')
        require(finite==Q(claim['finite_L1']['finite_ratio_lower']),'Finite L1 fraction mismatch.')
        require(finite>interval[1],'Finite L1 ratio does not clear the complete frozen reference upper.')
        if power==5: require(finite>Q(2135677,1000000),'Claimed simple strict bound failed.')
        output['seeds'].append({'power':power,'raw_positive_boxes':raw_count,
            'distinct_positive_boxes':len(inner),'all_source_coordinate_checks':checks,
            'memo_states_each_direction':states,'memo_states_both_directions':budget[0],
            'union_lower':str(volumes[0]),'union_upper':str(volumes[1]),
            'finite_L1_ratio_lower':str(finite),'finite_ratio_above_reference_upper':True,
            'wall_seconds':time.time()-t,'cpu_seconds':time.process_time()-cpu})
    require(Q(output['seeds'][1]['union_lower'])>Q(output['seeds'][0]['union_upper']),'Seed order separation failed.')
    require({name:sha(candidate/name) for name in lock['candidate_files']}==candidate_before,'Candidate changed during recheck.')
    require({name:sha(HERE/name) for name in manifest['files']}==audit_before,'Historical audit changed during recheck.')
    if reference_path is not None: require(sha(reference_path)==reference_hash,'Original R075 reference changed during recheck.')
    output.update({'all_locked_inputs_unchanged_after_run':True,
                   'seed5_finite_ratio_strictly_above_2135677_over_1000000':True,
                   'seed5_ideal_lower_strictly_above_seed4_ideal_upper':True,
                   'runtime':{'executable':sys.executable,'python':sys.version,'platform':platform.platform()}})
    return output

def main():
    parser=Parser(description=__doc__)
    parser.add_argument('--candidate-root',type=Path,default=HERE.parents[2]/'lower'/'round001',
                        help='candidate lower/round001 directory (default: relative to this project)')
    parser.add_argument('--reference',type=Path,
                        help='optional complete original R075 results.json; absent means frozen-interval comparison only')
    try:
        args=parser.parse_args()
        print(json.dumps(recheck(args.candidate_root,args.reference),indent=2),flush=True)
        return 0
    except Exception as exc:
        print(json.dumps({'status':'FAIL portable read-only recheck','error':str(exc)}),flush=True)
        return 1

if __name__=='__main__': sys.exit(main())
