#!/usr/bin/env python3
"""Compile/run only the three registered new graph moments; never old space batches.

All binaries, journals and raw output stay in --work-dir (outside Git).
If a 120-second run checkpoints, repeat with the same --work-dir --resume.
"""
from fractions import Fraction as F
from itertools import product, combinations_with_replacement
from math import comb, factorial
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import time


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rising(a,n):
    answer=F(1)
    for j in range(n): answer*=a+j
    return answer
def binomial(a,n):
    answer=F(1)
    for j in range(n): answer*=(a-j)/(j+1)
    return answer


def independent_single(n):
    beta=F(3,4)
    den=[binomial(beta,k)*(-1)**k*F(1+(-1)**k,2) for k in range(n+1)]
    inv=[F(1)]
    for k in range(1,n+1): inv.append(-sum((den[j]*inv[k-j] for j in range(1,k+1)),F()))
    return inv[n]*factorial(n)/rising(beta,n)


def independent_triangle(counts=(6,6,6)):
    # Average the four distinct sign groups (A,B,AB), not graph cumulants.
    beta=F(3,4)
    indices=sorted(product(*(range(k+1) for k in counts)),key=sum)
    den={(0,0,0):F(1)};inverse={(0,0,0):F(1)}
    for k in indices[1:]:
        n=sum(k)
        average=F(sum(A**k[0]*B**k[1]*(A*B)**k[2]
                      for A,B in product((-1,1),repeat=2)),4)
        multinomial=factorial(n)
        for j in k: multinomial//=factorial(j)
        den[k]=binomial(beta,n)*(-1)**n*multinomial*average
        inverse[k]=-sum((den[j]*inverse[tuple(a-b for a,b in zip(k,j))]
                         for j in product(*(range(a+1) for a in k))
                         if any(j) and den[j]),F())
    factor=1
    for k in counts: factor*=factorial(k)
    return inverse[counts]*factor/rising(beta,sum(counts))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--gmp-prefix',type=Path,default=Path('/opt/homebrew/Cellar/gmp/6.3.0'))
    parser.add_argument('--reference-kernel',type=Path)
    parser.add_argument('--results-path',type=Path)
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    registration=json.loads((root/'round012_registration.json').read_text())
    work=(args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix='round012-pilot-')))
    work.mkdir(parents=True,exist_ok=True)
    if any((p/'.git').exists() for p in (work,*work.parents)):
        raise ValueError('work directory must be outside every Git repository')
    header=root/'round012_graph_moments_gmp.hpp';cpp=root/'round012_graph_moments_cli.cpp'
    source_id=sha(header);journal=work/'moments.journal';binary=work/'graph_moments'
    if journal.exists() and not args.resume: raise ValueError('existing journal requires --resume')
    command=[shutil.which('clang++') or 'c++','-O2','-std=c++17','-Wall','-Wextra','-Werror',
             '-isystem',str(args.gmp_prefix/'include'),str(cpp),'-L'+str(args.gmp_prefix/'lib'),
             '-lgmpxx','-lgmp','-o',str(binary)]
    subprocess.run(command,check=True)
    requests=work/'requests.txt'
    requests.write_text(''.join(r['id']+' '+' '.join(map(str,r['counts']))+'\n'
                                for r in registration['pilot_requests']))
    runs=[]
    def run(max_new,label):
        output=work/(label+'.json')
        started=time.monotonic()
        completed=subprocess.run([str(binary),str(requests),str(journal),source_id,'120',str(max_new),str(output)],
                                 text=True,capture_output=True,timeout=125)
        if completed.returncode not in (0,75): raise RuntimeError(completed.stderr)
        raw=json.loads(output.read_text())
        runs.append({'label':label,'returncode':completed.returncode,'stdout':completed.stdout,
                     'wall_seconds':time.monotonic()-started,'result':raw})
        return raw
    if not args.resume:
        first=run(10,'forced_checkpoint')
        assert first['status']=='checkpoint' and first['cache']['computed']==10
        with journal.open('ab') as f: f.write(b'partial')
    final=run(0,'resume' if args.resume else 'after_recovery')
    if final['status']!='completed':
        checkpoint={'status':'checkpoint_continue_with_resume','work_dir':str(work),'runs':runs,
                    'source_id':source_id,'journal_sha256':sha(journal)}
        (work/'pilot_checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
        print(json.dumps(checkpoint,indent=2));raise SystemExit(75)
    if not args.resume:
        assert final['cache']['loaded']==10 and final['cache']['trailing_bytes_removed']==7
    values={r['id']:F(r['moment']) for r in final['results']}
    started=time.monotonic()
    single=independent_single(26);triangle=independent_triangle()
    assert values['single_edge_26']==single
    assert values['triangle_18']==triangle
    cross_seconds=time.monotonic()-started
    # New high degree only: a checked failure of the old signed-128-bit range.
    overflow_rows=[r['id'] for r in final['results'] if int(r['B']).bit_length()>127 or int(r['D']).bit_length()>127]
    assert 'single_edge_26' in overflow_rows
    partitions=list(combinations_with_replacement(range(13),5))
    partition_counts={str(n):sum(sum(p)==n for p in partitions) for n in range(24)}
    source_hashes={p.name:sha(p) for p in [header,cpp,root/'round012_registration.json',Path(__file__)]}
    if args.reference_kernel:
        assert sha(args.reference_kernel)=='2ccb0cae445f330c5155c70b90aa8c39a4e6ab994a45bc8494aed6918d0af9b3'
        source_hashes['frozen_astra004_graph_moments.cpp']=sha(args.reference_kernel)
    payload={
      'status':'three_new_high_degree_moments_and_recovery_certified',
      'scope':'No full high-degree spatial expansion, old 3389-moment baseline, Bernstein certificate, g or C5 improvement.',
      'source_sha256':source_hashes,'GMP_version':final['GMP_version'],
      'compile_command':command,'binary_sha256':sha(binary),'work_dir':str(work),
      'runs':runs,'independent_crosschecks':{
        'single_edge_26_two_group_fraction':str(single),'triangle_18_four_group_fraction':str(triangle),
        'both_exactly_equal':True,'elapsed_seconds':cross_seconds},
      'signed_128_bit_range_exceeded_by':overflow_rows,
      'quadratic_Z_scale':{
        'W22_raw_edge_compositions':comb(31,9),
        'W_powers_0_through_22_total_raw_compositions':comb(32,10),
        'degree26_raw_edge_compositions':comb(35,9),
        'spatial_partition_count_by_X_total_degree_max_coordinate12':partition_counts,
        'spatial_partitions_total_degree_up_to23_max_coordinate12':sum(partition_counts.values()),
        'full_face_Bernstein_grid_at_coordinate_degree12':13**4,
        'warning':'These are enumeration bounds/counts, not completed space coefficients or a runtime estimate.'},
      'checkpoint':{'journal':str(journal),'sha256':sha(journal),'bytes':journal.stat().st_size,
                    'protocol':'source-SHA256-locked append-only canonical B integer journal; per-record FNV checksum; partial trailing row removed; completed entries fsync checkpoints',
                    'resume_command':f'python -B {Path(__file__).name} --work-dir {work} --resume'},
      'frozen_reference_modified':False,'MATLAB_used':False}
    destination=args.results_path or root/'round012_pilot_results.json'
    destination.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'output':str(destination),'sha256':sha(destination),
                     'new_moments':final['results'],'cache':final['cache'],
                     'crosscheck_seconds':cross_seconds,'quadratic_Z_scale':payload['quadratic_Z_scale']},indent=2))


if __name__=='__main__': main()
