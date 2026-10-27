"""C160 exact preparation/gating helpers; importing performs no computation."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
RESULT156_SHA='f9a0689da14fb34d48c3ffb65b6abaef06193906ff71cc3cab77420074bee21f'
FIELD156=['budgets','finite_best','valid_true_beta2_upper']
CODE_FILES=['PROOF.md','common.py','prepare_for_terra.py','run_phaseA.py','run_phaseB.py','HANDOFF_TO_TERRA.md','STATIC_SOURCE_PINS.json']
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def req(ok,msg):
    if not ok:raise ArithmeticError(msg)
def new(p,name,x):
    with(Path(p)/name).open('x')as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def rows(paths):return [{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}for p in paths]
def check_rows(rs):
    for r in rs:req(sha(r['path'])==r['sha256'],'SOURCE_HASH:'+r['path'])
def check_package():
    for r in read(P/'PACKAGE_MANIFEST.json')['files']:req(sha(P/r['name'])==r['sha256'],'PACKAGE_HASH:'+r['name'])
    check_rows(read(P/'STATIC_SOURCE_PINS.json')['sources'])
def source_roles():return {r['role']:Path(r['path'])for r in read(P/'STATIC_SOURCE_PINS.json')['sources']}
def coeffs(model,candidate):
    x=[F(v,int(candidate['denominator']))for v in candidate['numerators']]
    req(len(x)==49 and min(x)>=0 and x[0]>0 and sum(x)==1 and int(candidate['denominator'])==1<<44,'EXACT_SIMPLEX')
    caps=[F(model['beta0'])]+[F(v)for row in model['beta_by_h']for v in row]
    req(len(caps)==49 and min(caps)>0,'POSITIVE49CAPS')
    a=[v/b for v,b in zip(x,caps)]
    req(sum(v*b for v,b in zip(a,caps))==1,'EXACT_NORMALIZED_CAP_COST')
    return x,caps,a
def model_preflight():
    check_package();s=source_roles();old=read(s['128_input']);r=read(s['156_result']);q=read(P/'quantization_certificate.json');inp=read(P/'inputs.json');m=inp['model_input']
    req(sha(s['156_result'])==RESULT156_SHA,'FULL156RESULT')
    b=F(r['budgets']['finite_best']['valid_true_beta2_upper']);D=1<<96;N=(b.numerator*D+b.denominator-1)//b.denominator
    req(r['h']=='1/128'and read(s['156_inputs'])['h']=='1/128'and b==F(r['budgets']['finite_best']['budget_expression_interval'][1]),'FULL156_H_AND_UPPER_ENDPOINT_IDENTITY')
    req(q['complete_source_upper_exact']==str(b)and q['source_result_sha256']==RESULT156_SHA and q['source_field_selector']==FIELD156,'EXACT156UPPER_PROVENANCE')
    req(int(q['D'])==D and int(q['N'])==N and (N-1)*b.denominator<D*b.numerator<=N*b.denominator,'EXACT_UPWARD_DYADIC96')
    req(F(q['betaD_exact'])==F(N,D)and 0<b<=F(N,D)<b+F(1,D),'CONSERVATIVE_CAP')
    req(old['h_values'][5]=='1/128'and all(old['Phi_table'][n][2]==n*(n-1)//2 for n in range(33))and inp['changed_flat_index']==32,'TRUE_BETA2_IDENTITY')
    for key in ['model','h_values','beta0','p_intervals','Phi_table','threshold_M_power6_exact','prefix_Kmax','prefix_state_count','tail_K_values']:req(m[key]==old[key],'UNCHANGED_'+key)
    req(len(m['beta_by_h'])==8 and all(len(row)==6 for row in m['beta_by_h']),'COMPLETE48CAPTABLE')
    for j,row in enumerate(m['beta_by_h']):
        for k,v in enumerate(row):
            if(j,k)==(5,1):req(F(v)==F(N,D)<F(old['beta_by_h'][j][k]),'SOLE_IMPROVED_GENUINE_UPPER_CAP')
            else:req(v==old['beta_by_h'][j][k],'NO_OTHER_CAP_CHANGED')
    phi=m['Phi_table'];req(len(phi)==33 and all(len(row)==7 and row[0]==1 and all(isinstance(v,int)and v>=0 for v in row)for row in phi),'NONNEGATIVE_COMPLETE_FEATURES')
    req(all(phi[n+1][k]>=phi[n][k]for n in range(32)for k in range(7)),'MONOTONE_PHI')
    for side in (0,1):
        ps=[F(t[side])for t in m['p_intervals']];req(len(ps)==8 and 0<ps[0]and ps[-1]<1 and all(a<b for a,b in zip(ps,ps[1:])),'WHOLE_CDF_LAW')
        req(all((p*(1<<192)).denominator==1 for p in ps),'UNCHANGED192GRID_CDF')
    req(all(F(a)<=F(b)for a,b in m['p_intervals']),'CDF_COUPLING')
    req(inp['warm_start_49']==read(s['128_candidate']),'ORIGINAL_FULL128_WARM_START')
    coeffs(m,inp['warm_start_49'])
    req(m['prefix_Kmax']==16 and m['prefix_state_count']==735471 and m['tail_K_values']==list(range(17,33)),'UNCHANGED_COMPLETE_PREFIX_AND_TAIL')
    theta=F(8,9)**6*24576**5/65000**6
    req(theta==F(m['threshold_M_power6_exact'])and F(inp['threshold_gstar_power6_exact'])==F(8,9)**6*24576**5/63500**6,'COMPLETE_THRESHOLDS')
    return inp,s
def check_source156_acceptance(release):
    audit=release.get('source156_saved_audit',{})
    req(audit.get('result_sha256')==RESULT156_SHA and audit.get('independent_saved_review_passed')is True and audit.get('main_accepted')is True,'SOURCE156_SAVED_AUDIT_NOT_ACCEPTED')
    req(isinstance(audit.get('evidence_files'),list)and len(audit['evidence_files'])>=1,'SOURCE156_SAVED_AUDIT_EVIDENCE_REQUIRED')
    check_rows(audit['evidence_files'])
def phase_requirements(phase):
    d=P/phase;reg=read(d/'registration.json')
    result={'round':'160','phase':phase,'authorized':True,'max_science_runs':1,'workers':1 if phase=='phaseA'else 15,'nice_minimum':10,'time_or_CPU_limit':None,'automatic_restarts':0,'source156_result_sha256':RESULT156_SHA,'source156_saved_audit_required':True,'base128_complete_caps_sha256':sha(source_roles()['128_input']),'source_base_complete_caps_accepted':True,'package_manifest_sha256':sha(P/'PACKAGE_MANIFEST.json'),'inputs_sha256':sha(P/'inputs.json'),'quantization_certificate_sha256':sha(P/'quantization_certificate.json'),'registration_sha256':sha(d/'registration.json'),'run_sha256':sha(P/('run_phaseA.py'if phase=='phaseA'else'run_phaseB.py')),'new_kernel_erfc_native_density_or_old_functional_calls':0,'new_C6_bound_permitted':False}
    result['scope']='ONE_NEW_GENUINE_CAP_NORMALIZED_FLOAT_SELECTOR'if phase=='phaseA'else'ONE_FIXED_NEW_NORMALIZED_FULL_H_CERTIFICATE'
    if phase=='phaseA':result.update({'max_objective_gradient_calls':150,'optimizer_invocations':1})
    else:result.update({'selection_sha256':sha(d/'selection.json'),'phaseA_result_sha256':reg['phaseA_result_sha256'],'candidate_sha256':sha(d/'candidate.json'),'jobs_sha256':sha(d/'jobs.json'),'logical_root_requests':735487,'arithmetic_bits':192,'raw_a_quantized':False})
    return result
def gated(phase):
    d=P/phase;release=read(d/'parent_release.json')
    for k,v in phase_requirements(phase).items():req(release.get(k)==v,'RELEASE:'+k)
    req(release.get('main_approved')is True and release.get('independent_review_approved')is True,'SEPARATE_EXACT_HASH_DUAL_GATE')
    check_source156_acceptance(release)
    req(not(d/'started.json').exists(),'NO_AUTOMATIC_RESTART')
    req(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']),'ONE_THREAD_ENV')
    os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
    inp,s=model_preflight();reg=read(d/'registration.json')
    for k,path in [('package_manifest_sha256',P/'PACKAGE_MANIFEST.json'),('inputs_sha256',P/'inputs.json'),('quantization_certificate_sha256',P/'quantization_certificate.json')]:req(reg[k]==sha(path),'FIXED_PREPARED_'+k)
    return inp,s,reg,release
