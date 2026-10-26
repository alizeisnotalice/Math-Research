#!/usr/bin/env python3
"""UNRUN until exact primary release; three fixed pure-exact transfer groups."""
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,177)
from pathlib import Path
import hashlib
import json
import os
import resource
import sys
import time

resource.setrlimit(resource.RLIMIT_CPU,(20,20))
signal.signal(signal.SIGXCPU,signal.SIG_DFL)
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
load=lambda p:json.loads(Path(p).read_text())
release_path=H/'parent_release.json'
if not release_path.exists():raise RuntimeError('UNRELEASED105: source/proof/static review only')
release=load(release_path)
if not (release.get('authorized') is True and str(release.get('round'))=='105'
        and release.get('script_sha256')==sha(Path(__file__))
        and release.get('scope_sha256')==sha(H/'scope_v2.json')
        and release.get('source_manifest_sha256')==sha(H/'source_manifest.json')):
    raise RuntimeError('105 release is not bound to this exact source/scope/script')
if any((H/n).exists() for n in ('started.json','results.json','run_receipt.json')):
    raise RuntimeError('105 exclusive once-only output already exists')

# Mathematical modules and arithmetic begin only after the primary release gate.
from fractions import Fraction as Q
from math import comb
START,CPU=time.monotonic(),time.process_time()
with (H/'started.json').open('x') as f:
    json.dump({'pid':os.getpid(),'round':'105','release_sha256':sha(release_path),
      'hard_fuse_seconds':177,'SIGALRM_is_default':signal.getsignal(signal.SIGALRM)==signal.SIG_DFL},f,indent=2);f.write('\n')
checks=0
out={'status':'RUNNING','round':'105','groups':[],'root_certificates':[],
     'kernel_native_old_grid_calls':0,'old104_100_executions':0,'children':0,
     'MATLAB_calls':0,'lower_bound_project_work':0,'automatic_retries':0,
     'new_Cn_claim':False,'108_proof_accepted':False}
def ck(ok,message):
    global checks
    checks+=1
    if not ok:raise AssertionError(message)
    if time.process_time()-CPU>=18:raise TimeoutError('105 CPU prelimit')
def group(i,name):
    ck(i==len(out['groups'])+1<=3,'fixed group cap/order')
    out['groups'].append({'id':i,'name':name,'pass':True})
def enclosure(a,b,formula):
    ck(0<a<=b,'positive scalar enclosure')
    return {'formula':formula,'lower':str(a),'upper':str(b),'display_only':[float(a),float(b)]}
def root(rad,n,label):
    rad=Q(rad);D=10**60
    ck(rad>0 and n in (2,5),'fixed positive root interface')
    target=(rad.numerator*D**n)//rad.denominator
    lo,hi=0,1<<((target.bit_length()+n-1)//n+1)
    while hi-lo>1:
        m=(lo+hi)//2
        if m**n<=target:lo=m
        else:hi=m
    a,b=Q(lo,D),Q(lo+1,D)
    ck(a**n<=rad<=b**n,'exact outward root certificate '+label)
    out['root_certificates'].append({'label':label,'radicand':str(rad),'degree':n,'lower':str(a),'upper':str(b)})
    return a,b
try:
    manifest=load(H/'source_manifest.json');scope=load(H/'scope_v2.json');inp=load(H/'coefficient_inputs.json')
    for ent in manifest['local_source_files']:
        ck(sha(ent['path'])==ent['sha256'],'source hash '+ent['path'])
    for name,expected in manifest['prepared_artifact_pins'].items():
        ck(sha(H/name)==expected,'prepared artifact hash '+name)
    ck(scope['maximum_new_exact_groups']==3 and inp['root_decimal_digits']==60,'fixed scope')
    old104=load(inp['source104_results_path'])
    ck(sha(inp['source104_results_path'])==inp['source104_results_sha256'],'104 exact saved result pin')
    ck(old104['configurations'][0]['actual_full_symbol_coercivity_reused_from098']==inp['old_coercivity'],'104 old coefficient identity only')
    pi_l,pi_u=map(Q,inp['pi_interval'])
    ck(Q(3)<pi_l<pi_u<Q(22,7),'inherited exact pi enclosure')

    # Group1: exact constants in the prepared analytic cap proof.
    delta=Q(inp['delta']);q=Q(inp['q'])
    ck(delta==Q(1,8) and q==Q(6,5),'fixed delta/q')
    count=4*comb(6,2)
    ck(count==inp['equality_vector_count']==60,'all signed two-coordinate equality directions')
    cap=delta**2/4;union=count*cap;outside=1-union
    ck((cap,union,outside)==(Q(1,256),Q(15,64),Q(49,64)),'arbitrary3D chordal-cap and union fractions')
    ck(248**2>2*175**2,'sqrt2 upper rational comparison for cstar')
    ck(Q(20)-14*Q(248,175)==Q(4,25),'cstar squared exceeds1/400')
    ck(2*114**2>160**2,'relative gain sqrt2/160 exceeds1/114')
    eta=union+outside*Q(114,115)**3
    ck(Q(31,32)<eta<1,'eta positive and improves uniform peak')
    r2=Q(8,27)
    ck((Q(1024,3)/512)**3==r2,'old/sharp Q3 coefficient ratio squared')
    ck(Q(3,128**2)*r2==Q(2,192**2),'principal K coefficient exact simplification')
    ck(Q(9*6,64**2)*r2==Q(1,16**2),'104 D coefficient exact simplification')
    out['cap_certificate']={'delta':str(delta),'single_cap_fraction_upper':str(cap),
      'union_fraction_upper':str(union),'outside_fraction_lower':str(outside),
      'eta':str(eta),'cstar_strict_lower':'1/20','relative_gain_strict_lower':'1/114',
      'arbitrary3D_cap_proof':'PROOF.md section3; projected equality vector may lie outside subspace',
      'Q3_105_upper':'eta/[pi^2*(16sqrt2)^3]',
      'principal_K_sharp':'sqrt2/(192*pi^4)','principal_K_105':'eta*sqrt2/(192*pi^4)',
      'D104_sharp':'Sigma/(16*pi^4)','D104_105':'eta*Sigma/(16*pi^4)'}
    group(1,'rational_square_constants_cap_fractions_and_exact_coefficient_identities')

    # Group2: exact scalar rescaling factors; no old budget is regenerated.
    rk=root(r2,2,'sharp_norm_factor')
    ik=root(r2**3,5,'sharp_I_factor')
    ia=root(eta**6*r2**3,5,'105_cap_I_factor')
    ra=(eta*rk[0],eta*rk[1])
    ck(ra[1]<rk[0]<rk[1]<1 and ia[1]<ik[0]<ik[1]<1,'strict scalar improvements')
    out['old104_formal_budget_scaling']={
      'sharp_norm_factor':enclosure(*rk,'(2/3)^(3/2)'),
      '105_cap_norm_factor':enclosure(*ra,'eta*(2/3)^(3/2)'),
      'sharp_I_factor':enclosure(*ik,'(2/3)^(9/5)'),
      '105_cap_I_factor':enclosure(*ia,'eta^(6/5)*(2/3)^(9/5)'),
      'meaning':'fixed104 P/B/full norm upper expressions multiply by norm factor; fixed I allocations multiply by I factor',
      'actual_G_or_I_lower_bound':False,'numeric_complete104_budget_recalculated':False,
      '104_acceptance_dependency':'separate;105 does not certify104'}
    group(2,'fixed_old104_to105_formal_norm_and_I_rescaling_factors_only')

    # Group3: one conditional BL coefficient comparison, not a108 proof.
    sq3=root(Q(3),2,'sqrt3_for_conditional_BL_factor')
    rb=(8/(3*sq3[1]*pi_u),8/(3*sq3[0]*pi_l))
    ck(eta>Q(31,32) and 93**2>8*32**2,'conditional BL analytic strict comparison')
    ck(rb[1]<ra[0],'conditional BL coefficient strictly beats105 cap coefficient')
    ratio=(rb[0]/ra[1],rb[1]/ra[0])
    ck(ratio[1]<1,'conditional dominance separated interval')
    out['conditional108_comparison']={
      'premise':'IF108 proves Q3<=1/(4096*pi^3) uniformly on all actual3D subspaces',
      'premise_accepted_or_proved_here':False,
      'BL_over_old098_factor':enclosure(*rb,'8/(3sqrt3*pi)'),
      'BL_over_105_cap_factor':enclosure(*ratio,'2sqrt2/(pi*eta)'),
      'strictly_smaller_than_all_planned105_coefficients':True,
      'budget_consequence':'by positive homogeneity, conditionally strictly smaller fixed104 norm and I allocations; no old104 recomputation needed',
      'full104_budget_recalculation_omitted':True,
      'scope':'conditional coefficient algebra only; not an independent108 certification'}
    group(3,'one_conditional_unproved108_coefficient_comparison_and_stop')
    for ent in manifest['local_source_files']:
        ck(sha(ent['path'])==ent['sha256'],'post-run source immutability')
    out.update(status='PASS_105_BACKUP_TRANSFER_CONDITIONAL108_DOMINANCE_NO_FULL104_RECALC',
      proof_qualification='external Xie_v1 Theorem1 source input;FW26 theorem1 proof read, classical originals not fully read;needed Kershaw special case proved in PROOF AppendixA',
      remaining_gaps=['external105 theorem/transfer independent review','108 premise not proved here','global g_L','complete same-model I_U','far/complement and complete positive S'],
      no_new_Cn=True)
except Exception as exc:
    out.update(status='FAILED_FIRST_FAILURE_STOP',error=repr(exc))
finally:
    out.update(assertions=checks,CPU_seconds_before_save=time.process_time()-CPU,
      wall_seconds_before_save=time.monotonic()-START,source_manifest_sha256=sha(H/'source_manifest.json'),
      script_sha256=sha(Path(__file__)),release_sha256=sha(release_path))
    with (H/'results.json').open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    receipt={'status':out['status'],'pid':os.getpid(),'completed_groups':len(out['groups']),
      'assertions':checks,'CPU_seconds_before_receipt_write':time.process_time()-CPU,
      'wall_seconds_before_receipt_write':time.monotonic()-START,
      'SIGALRM_is_default':signal.getsignal(signal.SIGALRM)==signal.SIG_DFL,
      'fuse_remaining_before_receipt_write':signal.getitimer(signal.ITIMER_REAL)[0],
      'hard_fuse_never_cancelled':True,'CPU_hard_seconds':20,'children':0,
      'automatic_retries':0,'results_sha256':sha(H/'results.json')}
    with (H/'run_receipt.json').open('x') as f:json.dump(receipt,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(receipt,ensure_ascii=False),flush=True)
# The default177-second fuse remains armed through final writes/stdout/exit.
sys.stdout.flush();sys.stderr.flush()
os._exit(0 if out['status']=='PASS_105_BACKUP_TRANSFER_CONDITIONAL108_DOMINANCE_NO_FULL104_RECALC' else 1)
