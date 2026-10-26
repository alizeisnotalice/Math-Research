"""Independent saved-only B055 postrun audit. Never imports producer code."""
import sys
sys.dont_write_bytecode = True
import datetime
import hashlib
import json
import os
import subprocess
import time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

P = Path(__file__).resolve().parent
B = P.parent / 'sparse4_probe_precheck_review'
def load(p):
    return json.loads(Path(p).read_text())
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc_interval(lo, hi):
    return {'lower_exact': str(lo), 'upper_exact': str(hi), 'width_exact': str(hi-lo)}
def decimals(lo, hi, places=50):
    scale = 10 ** places
    l = lo.numerator * scale // lo.denominator
    u = -((-hi.numerator * scale) // hi.denominator)
    def fmt(n):
        sign = '-' if n < 0 else ''
        n = abs(n)
        return f'{sign}{n//scale}.{n%scale:0{places}d}'
    return {'lower': fmt(l), 'upper': fmt(u)}

def check():
    cfg=load(P/'input_manifest.json')
    C=Path(cfg['producer_directory'])
    reg=load(P/'registration.json')
    checks=[]
    def pin(group,path,digest):
        assert sha(path)==digest,(group,str(path),'sha256 mismatch')
        checks.append((group,str(path),digest))
    for name,h in cfg['audit_source_sha256'].items():
        pin('auditor_source',P/name,h)
    for path,h in cfg['additional_readonly_sha256'].items():
        pin('additional_manual_interface_source',path,h)
    # 1. Hash chain: do not execute either producer or former independent checker.
    pin('full_manifest',C/'full_manifest.json',reg['producer_full_manifest_sha256'])
    pin('full_report',C/'full_report.md',reg['producer_full_report_sha256'])
    full=load(C/'full_manifest.json')
    for name,h in full['whitelist_sha256'].items():pin('full_whitelist',C/name,h)
    pin('precheck_manifest',C/'precheck_manifest.json',full['precheck_manifest_sha256'])
    pre=load(C/'precheck_manifest.json')
    for name,h in pre['whitelist_sha256'].items():pin('precheck_local',C/name,h)
    for path,h in pre['external_readonly_sha256'].items():pin('precheck_external',path,h)
    pin('B055_manifest',B/'manifest.json',reg['precheck_manifest_sha256'])
    assert full['B055_manifest_sha256']==reg['precheck_manifest_sha256']
    bm=load(B/'manifest.json')
    for name,h in bm['whitelist_sha256'].items():pin('B055_sealed_file',B/name,h)
    for path,h in load(B/'input_manifest.json')['sha256'].items():pin('B055_original_input',path,h)
    for path,h in load(B/'header_binding.json')['headers_sha256'].items():pin('compile_header',path,h)
    build=load(C/'build_record.json')
    pin('binary',C/'scratch/sparse4_probe',build['binary_sha256'])
    assert build['binary_sha256']==pre['binary_sha256']
    assert bm['status']=='FROZEN_INDEPENDENT_PASS_PRECHECK_STRUCTURE_AND_SOURCE_ONLY'
    # 2. One completed native, release and tool-session receipts.
    e=load(C/'full_execution.json');s=load(C/'tool_session_receipt.json');r=load(C/'full_release.json')
    assert e['attempt']==s['attempt']==full['native_attempts']==r['max_native_attempts']==1
    assert s['native_restarts']==full['native_restarts']==0
    assert e['phase']=='EXITED' and e['exit_code']==s['terminal_tool_exit_code']==0
    assert s['status']=='COMPLETE_SAME_SESSION_EXIT_0'
    assert e['native_pid']==s['native_pid']==full['native_pid']==59613
    assert s['tool_session_id']==full['tool_session_id']==70248
    assert s['launch_chunk_id']=='f43ca9' and s['terminal_chunk_id']=='ea43f8'
    assert s['launch_command']=='python3 launch_once.py'
    assert e['command']==[str(C/'scratch/sparse4_probe'),'full',str(C/'native_input.txt'),str(C/'full_checkpoint.jsonl'),str(C/'full_token.txt')]
    assert r['status']=='MAIN_RELEASED_ONE_NATIVE_AFTER_B055' and r['B055_review_accepted_by_main'] is True
    assert r['B055_manifest_sha256']==sha(B/'manifest.json')
    assert Path(r['B055_manifest_path'])==B/'manifest.json'
    assert r['precheck_manifest_sha256']==e['precheck_manifest_sha256']==sha(C/'precheck_manifest.json')
    assert r['native_timeout_seconds']==e['timeout_seconds']==1800
    assert 0<e['process_wall_seconds']==s['process_wall_seconds']==full['process_wall_seconds']<1800
    assert r['parameters']=={'beta':'1/80','point_numerators':[50,50,50,45,27,0],'point_denominator':50,'T':'32/825','Arb_bits':512}
    for key,name in [('main_release_sha256','full_release.json'),('source_sha256','sparse4_probe.cpp'),('input_sha256','native_input.txt'),('stdout_sha256','native_full.stdout.json'),('stderr_sha256','native_full.stderr.txt')]:assert e[key]==sha(C/name)
    assert e['binary_sha256']==build['binary_sha256']
    assert (C/'native_full.stderr.txt').read_bytes()==b''
    assert (C/'full_token.txt').read_text()=='C055_MAIN_RELEASED_FINAL_DEPENDENCIES_CHECKED\n'
    release_time=datetime.datetime.fromisoformat(r['time_utc']).timestamp()
    assert release_time<e['start_unix']
    assert e['start_unix']+e['process_wall_seconds']<datetime.datetime.fromisoformat(full['created_utc']).timestamp()
    ps=subprocess.run(['ps','-p',str(e['native_pid']),'-o','pid=,stat=,command='],capture_output=True,text=True)
    assert not ps.stdout.strip(),('native PID currently exists',ps.stdout)
    # 3. Every saved checkpoint, exact field counts and monotone progress.
    n=load(C/'native_full.stdout.json');cp=[json.loads(x) for x in (C/'full_checkpoint.jsonl').read_text().splitlines()]
    assert len(cp)==61 and n['status']=='COMPLETE_ONE_JOINT_SCALAR_PROBE' and n['precision_bits']==512
    assert 0<n['kernel_seconds']<=e['process_wall_seconds']
    last_seconds=0;last_nodes=0
    for i,row in enumerate(cp[:-1],1):
        images=32*i;terms=images*13229
        assert row['phase']=='running' and row['complete'] is False and row['images_complete']==images
        assert row['paired_image_cone_terms']==row['G_DD5_evaluations']==row['H_DD5_evaluations']==terms
        assert row['logical_terms_each_field']==24*terms and row['total_DD5_evaluations']==2*terms
        assert last_seconds<row['seconds']<=n['kernel_seconds']
        assert last_nodes<=row['unique_nodes']<=n['unique_nodes']
        last_seconds=row['seconds'];last_nodes=row['unique_nodes']
    assert cp[-2]['images_complete']==1920 and cp[-2]['unique_nodes']==n['unique_nodes']==303744
    final=cp[-1]
    assert final['phase']=='computed' and final['complete'] is True
    for row in [final,n]:
        assert row['paired_image_cone_terms']==row['G_DD5_evaluations']==row['H_DD5_evaluations']==25399680
        assert row['total_DD5_evaluations']==50799360
    assert n['G_logical_terms']==n['H_logical_terms']==609592320==25399680*24
    assert final['seconds']==cp[-2]['seconds']==n['kernel_seconds']
    assert 0<=n['zero_node_paired_terms']==259937<=25399680
    assert 0<=n['confluent_paired_terms']==84987<=25399680
    comp=load(C/'point_comparison.json')
    assert comp['zero_node_logical_terms_each_field']==24*n['zero_node_paired_terms']==6238488
    assert comp['confluent_logical_terms_each_field']==24*n['confluent_paired_terms']==2039688
    # 4. Integer common-unit reconstruction, not producer's Fraction comparison implementation.
    # Native exponents are -511, -508, -512, hence common unit 2^-512.
    u=1<<512;iv={};integer_iv={}
    for field,exponent in [('G',-511),('H_T',-508),('K_T',-512)]:
        v=n[field];assert int(v['common_binary_exponent'])==exponent
        a,b=int(v['lower_mantissa']),int(v['upper_mantissa'])
        assert str(a)==v['lower_mantissa'] and str(b)==v['upper_mantissa'] and a<=b
        lo=a<<(exponent+512);hi=b<<(exponent+512)
        integer_iv[field]=(lo,hi);iv[field]=(Q(lo,u),Q(hi,u))
    gl,gu=integer_iv['G'];hl,hu=integer_iv['H_T'];kl,ku=integer_iv['K_T']
    sublo=gl-hu;subhi=gu-hl
    assert kl<=sublo<=subhi<=ku
    assert 0<hl<=hu<gl<=gu and 0<kl<=ku
    # 5. Exact threshold and margin certificate; all decimal display is directed by integer division.
    assert Q(32,825)/4096==Q(1,105600)
    gap_numerator=u-105600*ku
    assert gap_numerator>0
    threshold=Q(1,105600);margin=(threshold-iv['K_T'][1],threshold-iv['K_T'][0])
    assert margin[0]==Q(gap_numerator,105600*u)>0
    assert comp['threshold_exact']=='1/105600'
    assert comp['classification']==full['classification']=='EXCLUDES_ONLY_THIS_FIXED_FRAME_TIME_NECESSARY_CONDITION'
    for key in ['global_mT_certified','global_g_or_I_certified','new_Cn_bound']:
        assert comp[key] is False and full[key] is False
    # 6. Bind the already proved true-S2 geometry and the fixed-cube interface; no group/pair reenumeration.
    old=load(B/'results.json');s2=load(C/'S2_certificate.json')
    assert old['true_S2']['minimum_squared_pair_distance']==s2['pair_min_squared']=='363/4375'
    assert old['true_S2']['true_S2_rejection_gate_passed'] is True
    assert old['group']['negative_representatives']==1920 and old['group']['representative_multiplicity']==24
    assert old['normalization']['T_exact']=='32/825'
    assert old['normalization']['H_normalizer']==old['normalization']['core_G0_normalizer']=='64*pi^6'
    point=list(map(Q,[1,1,1]))+[Q(9,10),Q(27,50),Q(0)]
    assert all(-1<=x<=1 for x in point)
    theory=Path(cfg['continuity_theorem']).read_text()
    assert 'also proves that $G$ is continuous' in theory
    interface=Path(cfg['fixed_time_interface']).read_text()
    assert 'm_T=ess inf_[-1,1]^6 K_T^beta>0' in interface and 'cost T/(64m_T)' in interface
    source=(C/'sparse4_probe.cpp').read_text()
    assert '!(y<neg)' in source and 'arb_sub(K.a,G.a,H.a,PREC)' in source
    # Digest every frozen external object again to exclude edits during this audit.
    for group,path,h in checks:assert sha(path)==h,(group,path,'changed during saved audit')
    assert sha(C/'full_manifest.json')==reg['producer_full_manifest_sha256']
    assert sha(B/'manifest.json')==reg['precheck_manifest_sha256']
    intervals={name:enc_interval(*bounds) for name,bounds in iv.items()}
    intervals['threshold_minus_K_T']=enc_interval(*margin)
    dec={name:decimals(*bounds) for name,bounds in {**iv,'threshold_minus_K_T':margin}.items()}
    return {'status':'PASS_SAVED_ONLY_FIXED_FRAME_TIME_OBSTRUCTION','configurations':6,'hash_checks':len(checks),'hash_check_groups':dict(sorted(Counter(g for g,_,_ in checks).items())),'native_execution':{'attempt':1,'native_pid':59613,'native_pid_observed_absent':True,'session_id':70248,'launch_chunk_id':s['launch_chunk_id'],'terminal_chunk_id':s['terminal_chunk_id'],'exit_code':0,'process_wall_seconds':e['process_wall_seconds'],'kernel_seconds':n['kernel_seconds'],'historical_execution_evidence':'frozen cross-consistent receipts; no historical OS reconstruction'},'checkpoint':{'lines':61,'running_records':60,'computed_records':1,'completed_images':1920,'cones_per_image':13229,'paired_terms':25399680,'logical_terms_each_field':609592320,'G_DD5':25399680,'H_DD5':25399680,'total_DD5':50799360,'zero_node_paired_terms':259937,'confluent_paired_terms':84987,'unique_nodes':303744,'monotone_saved_progress':True},'exact_intervals':intervals,'outward_decimals_50_places':dec,'integer_certificate':{'common_unit':'2^-512','G_endpoints':list(map(str,integer_iv['G'])),'H_T_endpoints':list(map(str,integer_iv['H_T'])),'K_T_endpoints':list(map(str,integer_iv['K_T'])),'K_lower_padding_from_G_minus_H_units':str(sublo-kl),'K_upper_padding_from_G_minus_H_units':str(ku-subhi),'strict_comparison':'105600*K_upper_integer < 2^512','positive_integer_gap':str(gap_numerator)},'K_encloses_saved_G_minus_H':True,'threshold_exact':'1/105600','K_upper_strictly_below_threshold':True,'off_true_S2_squared_distance':'363/4375','manual_continuity_interface':'K=G-H_T continuous near x0; cube-relative neighborhoods have positive six-dimensional measure, including boundary x0. Hence essinf_cube K <= K(x0) <= certified K_upper < T/4096.','negative_representatives':'Inherited complete precheck; odd scalar data and DD order five give reflection factor (-1)^6=1; weight24 once.','scope':{'beta':'1/80','T':'32/825','all_T_excluded':False,'all_beta_excluded':False,'original_Rc_excluded':False,'new_Cn_bound':False,'global_floor_certified':False,'fixed_time_essential_upper_obstruction':True},'new_scalar_native':0,'old_precheck_configurations_rerun':0,'group_pair_reenumerations':0,'new_matrix_operations':0,'author_writer_runs':0,'MATLAB':0}

def main():
    assert sys.argv[1:]==['--run']
    assert not (P/'execution_receipt.json').exists(),'one saved-only registered attempt; do not rerun'
    receipt={'status':'RUNNING','pid':os.getpid(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_manifest_sha256':sha(P/'input_manifest.json'),'new_scalar_native':0,'old_precheck_configurations_rerun':0,'group_pair_reenumerations':0}
    (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    tic=time.monotonic()
    try:
        result=check();elapsed=time.monotonic()-tic
        assert elapsed<=600
        (P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
        receipt.update(status=result['status'],exit_code=0,configurations=6,actual_seconds=elapsed,results_sha256=sha(P/'results.json'))
    except BaseException as err:
        receipt.update(status='STOP_SAVED_AUDIT_FAILURE_NO_RETRY',exit_code=2,actual_seconds=time.monotonic()-tic,error=repr(err));raise
    finally:
        receipt['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
