#!/usr/bin/env python3
"""Portable stdout-only exact endpoint/coverage audit. Never imports kernel.py."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,time,platform
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent;CUBE=HERE.parents[2]
p=argparse.ArgumentParser();p.add_argument('--numeric-dir',type=Path,default=CUBE/'upper/round012/kernel_full');p.add_argument('--pilot-dir',type=Path,default=CUBE/'upper/round011/kernel_pilot');p.add_argument('--theory-dir',type=Path,default=CUBE/'upper/round011/kernel_theory');p.add_argument('--point-source',type=Path,default=CUBE/'upper/round003/truncated_numeric/results.json');a=p.parse_args()
NUM=a.numeric_dir.resolve();PILOT=a.pilot_dir.resolve();THEORY=a.theory_dir.resolve();POINT=a.point_source.resolve()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
files=[x for x in sorted(NUM.iterdir())if x.is_file()]+[PILOT/n for n in ('manifest.json','run.py','results.json','testcases.json')]+[THEORY/n for n in ('kernel_proof.md','finite_contract.md')]+[POINT]
def snap():return {str(p):sha(p)for p in files}
start=time.perf_counter();before=snap()
for n,h in json.loads((NUM/'manifest.json').read_text())['files'].items():assert sha(NUM/n)==h,n
for line in (NUM/'prerun-sha256.txt').read_text().splitlines():
    h,n=line.split(maxsplit=1);named=Path(n)
    if not named.is_absolute():f=NUM/n
    elif named.name in ('kernel_proof.md','finite_contract.md'):f=THEORY/named.name
    elif named.name=='manifest.json' and 'kernel_pilot' in str(named):f=PILOT/'manifest.json'
    else:f=named
    assert sha(f)==h,n
copied=json.loads((NUM/'r011-manifest.json').read_text());actual=json.loads((PILOT/'manifest.json').read_text())
assert (NUM/'r011-results.json').read_bytes()==(PILOT/'results.json').read_bytes()
assert (NUM/'kernel.py').read_bytes()==(PILOT/'run.py').read_bytes()
for n,on in {'kernel.py':'run.py','r011-results.json':'results.json','cells.jsonl':'cells.jsonl','r009-rays.json':'r009-rays.json','r009-inputs.json':'r009-inputs.json','r010-manifest.json':'r010-manifest.json','testcases.json':'testcases.json'}.items():assert sha(NUM/n)==copied['files'][on]==actual['files'][on]
old=json.loads((NUM/'r011-results.json').read_text());T=json.loads((NUM/'testcases.json').read_text());R=json.loads((NUM/'result.json').read_text());receipt=json.loads((NUM/'receipt.json').read_text())
assert old['status']=='completed' and old['generic_planning_width_pass'] and old['promotions']==[]
assert R['status']==receipt['status']=='completed' and R['precision']==256 and receipt['precision_bits']==256
physical=json.loads(POINT.read_text())['minima'][0]
assert physical['lowest_sobol_index']==2710 and physical['point_floathex']==T['generic_case']['x_floathex']==R['point_floathex']
assert T['generic_case']['epsilon']==R['epsilon']=='1/1024'
exact_point=[str(Q.from_float(float.fromhex(h)))for h in R['point_floathex']]
assert all(Q(x)>0 for x in exact_point) and len({abs(Q(x))for x in exact_point})==7
pilot={x['chunk']:x for x in old['generic_chunks']};assert set(pilot)==set(range(4))
rows=[json.loads(s)for s in (NUM/'terminal.jsonl').read_text().splitlines()]
assert len(rows)==5040 and {x['chunk']for x in rows}==set(range(5040)) and [x['chunk']for x in rows[:4]]==[0,1,2,3]
assert all(x['status']=='completed' for x in rows)
def dy(d):return Q(int(d[0]))*Q(2)**int(d[1])
def ep(x):return[str(x.lower().fmpq()),str(x.upper().fmpq())]
keys=list(R['counts']);total={k:0 for k in keys};reused={k:0 for k in keys};new={k:0 for k in keys};L=Q(0);U=Q(0);negative=0;crossing=0;newcpu=Q(0)
for row in rows:
    j=row['chunk'];assert row['start_phase']==64*j and row['end_phase']==64*(j+1) and row['precision']==256
    assert row['source']==('reused_r011' if j<4 else 'new')
    counts=row['counts'];assert set(counts)==set(keys) and counts['top_requests']==64*212
    assert all(isinstance(v,int) and v>=0 for v in counts.values())
    assert 0<counts['unique_top']<=counts['top_requests'] and counts['log_evaluations']<=64*46
    if j<4:
        # Original file bytes are reused; normalized interval payload is exactly unchanged.
        assert row['interval']==pilot[j]['interval']
        assert [v['phase_id']for v in pilot[j]['phase_records']]==list(range(64*j,64*(j+1)))
        assert counts=={k:sum(v['counts'][k]for v in pilot[j]['phase_records'])for k in keys}
        assert row['cpu_seconds'] is None and row['wall_seconds'] is None and row['peak_rss_bytes'] is None
    else:
        assert row['cpu_seconds']>=0 and row['wall_seconds']>=0
        newcpu+=Q(str(row['cpu_seconds']))
    lo,hi=dy(row['interval']['lower']),dy(row['interval']['upper']);assert lo<=hi
    L+=lo;U+=hi;negative+=int(hi<0);crossing+=int(lo<=0<=hi)
    for k in keys:
        total[k]+=counts[k];(reused if j<4 else new)[k]+=counts[k]
assert total==R['counts'] and total['top_requests']==R['top_requests_before_cache']==68382720
assert R['chunks']==receipt['terminal_chunks']==5040 and R['new_chunks']==receipt['new_completed']==5036 and R['reused_chunks']==receipt['reused_chunks']==4
assert R['halfphases']==322560 and R['no_extra_normalization'] and R['not_a_global_g_certificate']
assert total['hermite_evaluations']==total['zero_ray_knots']==total['repeated_top_nodes']==0
reportedL,reportedU=dy(R['full_H']['lower']),dy(R['full_H']['upper'])
assert 0<reportedL<=L<=U<=reportedU and reportedU-reportedL<=Q(1,10**12)
assert dy(R['width']['upper'])>=reportedU-reportedL and R['width_pass'] and R['positive']
# The only Arb arithmetic:768-bit conversion of exact signed endpoint sums.
ctx.prec=768
low=arb(fmpq(L.numerator,L.denominator)).lower();high=arb(fmpq(U.numerator,U.denominator)).upper()
lowQ,highQ=Q(str(low.fmpq())),Q(str(high.fmpq()))
assert lowQ==L and highQ==U # these dyadic sums are exactly representable at768bits
assert reportedL<=lowQ<=highQ<=reportedU and 0<lowQ and highQ-lowQ<=Q(1,10**12)
# Resource bookkeeping, not a scientific performance or memory bound.
snapshots=[]
for line in (NUM/'resources.jsonl').read_text().splitlines():
    obj=json.loads(line);workers=[]
    for s in obj.get('processes',[]):
        b=s.split(None,4)
        if len(b)==5 and 'multiprocessing.spawn import spawn_main' in b[4]:workers.append((b[0],Q(b[2]),int(b[3])))
    if workers:snapshots.append({'workers':len(workers),'summed_pcpu':sum((z[1]for z in workers),Q(0)),'summed_rss_KiB':sum(z[2]for z in workers)})
assert snapshots and all(s['workers']==15 for s in snapshots)
parentcpu=Q(str(receipt['parent_cpu_seconds']));wall=Q(str(receipt['wall_seconds']));coreavg=(newcpu+parentcpu)/wall
assert snap()==before
out={'status':'PASS full-generic signed endpoint/coverage audit; no new phase','environment':{'python':platform.python_version(),'python_flint':flint.__version__,'final_endpoint_conversion_bits':768},'physical_point_exact_dyadic':exact_point,'physical_point_source_index':2710,'epsilon':'1/1024','kernel_frozen_bytes_equal_R011':True,'first4_original_result_file_bytes_and_interval_payloads_reused':True,'complete_chunk_coverage':[0,5039],'chunks':5040,'reused_chunks':4,'new_chunks':5036,'halfphases':322560,'normalization':'each saved chunk already includes2/(2pi)^7; final sum only, no second factor','negative_chunk_intervals':negative,'chunk_intervals_containing_zero':crossing,'exact_Fraction_sum_interval':[str(L),str(U)],'independent_768_outward_H_interval':[str(low.fmpq()),str(high.fmpq())],'independent_exact_width':str(U-L),'reported_full_H_interval':[str(reportedL),str(reportedU)],'reported_interval_contains_exact_signed_sum':True,'strict_positive_and_width1e12':True,'total_instrumented_counts_including_reuse':total,'reused_R011_instrumented_counts':reused,'new_R012_instrumented_counts':new,'count_scope':'unique counters are per-phase cache contexts, summed; not globally unique keys across phases','new_CPU_decimal_record_sum_exact':str(newcpu),'reported_new_CPU_seconds':str(receipt['sum_new_chunk_cpu_seconds']),'parent_CPU_seconds':str(receipt['parent_cpu_seconds']),'wall_seconds':str(receipt['wall_seconds']),'approx_combined_core_equivalents':float(coreavg),'worker_phase_snapshots':len(snapshots),'sample_worker_count':15,'sample_summed_CPU_percent_range':[str(min(s['summed_pcpu']for s in snapshots)),str(max(s['summed_pcpu']for s in snapshots))],'sample_summed_worker_RSS_KiB_range':[min(s['summed_rss_KiB']for s in snapshots),max(s['summed_rss_KiB']for s in snapshots)],'resource_caveat':'shared machine; summed RSS may count shared pages; neither15coresaturated nor systemwide peak is claimed','source_hashes':before,'source_before_after_equal':True,'new_DD_or_phase_evaluations':0,'whole_kernel_replays':0,'author_checker_executed':False,'new_MC':0,'global_g_proved':False,'new_I_budget':False,'new_C7_bound':False,'elapsed_audit_seconds':time.perf_counter()-start}
print(json.dumps(out,indent=2))
