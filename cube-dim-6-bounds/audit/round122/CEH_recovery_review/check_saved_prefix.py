"""One E122 independent saved/new-prefix verification; no old program."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import datetime,hashlib,json,os,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
start,cpu=time.monotonic(),time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
if any((P/n).exists() for n in ['started.json','review_results.json','run_receipt.json']):raise SystemExit('one attempt only')
reg=read(P/'registration.json')
assert sha(Path(__file__))==reg['script_sha256']
assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
save('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'registration_sha256':sha(P/'registration.json')})
count=0;status='FIRST_FAILURE_STOP';stage='source';failure=None
def ck(b,label):
    global count
    assert b,label
    count+=1
def saved(d):return tuple(Q(d[s+'_mantissa'])*Q(2)**int(d[s+'_binary_exponent']) for s in ['lower','upper'])
def contains(o,v):
    if o==v:return True
    if isinstance(o,dict):return any(contains(x,v) for x in o.values())
    if isinstance(o,list):return any(contains(x,v) for x in o)
    return False
try:
    src=read(P/'source_manifest.json');pins={}
    for row in src['files']:
        p=Path(row['path']);ck(sha(p)==row['sha256'],'source pin');ck(p.stat().st_size==row['bytes'],'source size');pins[str(p)]=sha(p)
    A=Path(reg['author_directory']);am=read(A/'manifest.json')
    ck(sha(A/'manifest.json')==reg['author_manifest_sha256'],'author final seal')
    for row in am['whitelist']:ck(sha(A/row['path'])==row['sha256'] and (A/row['path']).stat().st_size==row['bytes'],'full whitelist')
    inp=read(A/'upstream.json');out=read(A/'results.json');peaks,cap,h,l=out['configurations'];cert=read(A/'certificate.json')
    old62=read(Path(reg['source062']));oldreg=read(Path(reg['registration062']));old33=read(Path(reg['source033']));oldCE=read(Path(reg['source059account']));oldCEH=read(Path(reg['source062account']));old111=read(Path(reg['source111']));old111c=read(Path(reg['source111constants']))
    for path,digest in oldreg['source_sha256'].items():ck(sha(Path(path))==digest,'original062 ten source pins')
    fh=inp['high062']
    for k,v in fh['constants'].items():ck(old62['constants'][k]==v,'unchanged062 constant '+k)
    for k in ['Pshort','Plong','Bfull','epsilon']:ck(old62[k]==fh[k],'complete saved062 '+k)
    ck(old62['union']['full_time_high_norm_upper']==fh['union_norm'],'old full union norm')
    ck(old62['union']['simple_rational_upper']==fh['simple_high_B_upper'],'old062 simple H bound')
    for name in ['C_root','E_root']:ck(contains(oldCE,inp['ledger'][name]),'original059 C/E root copy')
    ck(contains(oldCEH,inp['ledger']['H_root']),'original062 H-account root')
    ck(oldCEH['paid_C_E_H_I_upper']==inp['ledger']['paid_CEH']==inp['paid_CEH_111'],'original062 paidCEH copy')
    ck(contains(old111c,inp['paid_CEH_111']),'original111 paidCEH copy')
    ref=old111['final']['comparisons']['ref108_frozen_reference']
    ck(inp['old_far_allocation']['lower']==ref['S_plan_lower'] and inp['old_far_allocation']['upper']==ref['S_plan_upper'],'full111 both endpoints')
    ck('Kref=1/(8pi5)' in ref['reference_assumption'] and ref['g0']=='unproved globally','111 qualified branch')
    rc=read(A/'run_receipt.json');obs=read(A/'observed_tool_receipt.json');tool=read(A/'actual_tool_receipt.json');rel=read(A/'parent_release.json');areg=read(A/'registration.json')
    ck(out['status']==rc['status']=='COMPLETE_FOUR_FIXED_CONFIGURATIONS' and rc['completed_configurations']==4 and rc['failure'] is None,'author four complete configurations')
    ck(rc['pid']==tool['science_PID']==73840 and obs['observer_pid']==tool['observer_PID']==73832,'author process lineage')
    ck(obs['exit_code']==tool['science_tool_exit_code']==0 and obs['child_reaped'] and obs['children_created']==1,'author complete observed exit')
    ck(tool['unique_science_run_count']==1 and tool['unique_configurations']==4 and tool['failures']==0,'one producer no continuation')
    ck(rel['script_sha256']==obs['science_entry_sha256']==sha(A/'run_prepared.py') and rel['registration_sha256']==sha(A/'registration.json'),'released science source')
    ck(rel['observer_sha256']==sha(A/'observe_once.py') and out['release_sha256']==obs['release_sha256']==sha(A/'parent_release.json'),'release observer chain')
    ck(read(A/'started.json')['pid']==73840 and read(A/'observer_started.json')['pid']==73832,'exclusive started records')
    ck(read(A/'checkpoint.json')['configurations']==out['configurations'],'four complete checkpoint')
    for n,digest in areg['prepared_pins'].items():ck(sha(A/n)==digest,'prepared pins')
    ck(cert['results_sha256']==sha(A/'results.json') and cert['upstream_sha256']==sha(A/'upstream.json'),'certificate result/input bindings')
    stage='saved_exact'
    for row,got in zip(inp['peaks033'],peaks['nine_cases']):
        old=next(x for x in old33['cases'] if x['id']==row['id'])
        ck(all(old[k]==row[k] for k in ['dimension','Gram_determinant','R_exact']),'033 source case unchanged')
        r=int(row['dimension']);D=Q(row['Gram_determinant']);rr=Q(row['R_exact']);ratio=(2+Q(2*r,3))**r/(16**(2*r)*D*rr**2)
        ck(ratio==Q(got['ratio_BL_over_exact_squared']) and Q(got['strict_residual'])==ratio-1>0,'exact direct-rank ratio')
    ck(len(inp['peaks033'])==len(peaks['nine_cases'])==9,'all nine cases')
    counts=[Q(29,32)*x['unordered_count'] for x in inp['triple_classes']]
    ck(counts==list(map(Q,[435,1740,435,580,1305])) and sum(counts)*6==26970,'no-ones exact class counts')
    C0=24576*sum(n*Q(row['R3']) for n,row in zip(counts,inp['triple_classes']))
    ck(C0==Q(44743079,2970)==Q(peaks['Crem_times_pi6_exact'])==Q(29,32)*Q(inp['B3_pi6_exact']),'restored Crem identity')
    d=Q(31,64);s=list(map(Q,cap['squares']));av=cap['time_coefficients'];p0=Q(cap['p0'])
    ck(s==[6*d*d,Q(10,3)*d*d,2*d*d] and av==[992,960,928] and Q(cap['C0'])==C0 and p0==Q(216,2**30),'unchanged gap and full-law peak constants')
    prodA=av[0]*av[1]*av[2]
    ck(C0>p0*prodA and all(x>0 for x in s),'strict ratio limit above1')
    def poly(u):
        z=Q(1)
        for si,ai in zip(s,av):z*=si+ai*u
        return C0*u**3-p0*z
    ul,uu=map(Q,cap['u_star_interval']);tl,tu=map(Q,cap['T_interval'])
    ck(0<ul<uu<1 and uu-ul==Q(1,2**192),'saved192-step crossing bracket')
    ck(poly(ul)<=0<=poly(uu) and list(map(Q,cap['u_polynomial_endpoint_signs']))==[poly(ul),poly(uu)],'cubic endpoint residuals')
    ck(0<tl and tl**2<=ul<uu<=tu**2 and Q(15,1000)<tl<tu<Q(16,1000),'exact T squares')
    beta=list(map(Q,cap['beta']));hs=list(map(Q,cap['partial_fraction_coefficients']))
    ck(beta==[si/ai for si,ai in zip(s,av)] and len(set(beta))==3 and min(beta)>0,'distinct true resolvent roots')
    for i in range(3):
        den=Q(1)
        for j in range(3):
            if i!=j:den*=beta[j]-beta[i]
        ck(hs[i]*den==1,'signed residue identity')
    ck(sum(hs)==0 and sum(hs[i]*sum(beta[j] for j in range(3) if j!=i) for i in range(3))==0 and sum(hs[i]*beta[(i+1)%3]*beta[(i+2)%3] for i in range(3))==1,'complete partial-fraction polynomial identity')
    ck(cap['B_actual_le_p_before_enlargement'] and not cap['B_bound_le_p_asserted'],'actual subdensity ordering')
    # Only NEW finite prefixes, independently rounded at640bits and using
    # atan half-angle. No full-B or old defining integral is evaluated.
    stage='new_prefix_Arb640'
    import flint
    from flint import arb,fmpq,ctx
    ctx.prec=640
    def A(x):
        x=Q(x);return arb(fmpq(x.numerator,x.denominator))
    def primitive(t):
        total=arb(0)
        for bi,hi in zip(beta,hs):
            root=A(bi).sqrt();x=A(t)/root
            theta=2*(x/(1+(1+x*x).sqrt())).atan()
            total+=A(hi)*theta/root
        return total/prodA
    ftl,ftu=primitive(tl),primitive(tu)
    fl,fu=map(Q,cap['finite_primitive_enclosure'])
    ck(A(fl)<ftl and ftu<A(fu),'independent new prefix verifies saved512 outward endpoints')
    def endpoint(x,upper):
        m,e=(x.upper() if upper else x.lower()).man_exp()
        return Q(int(m))*Q(2)**int(e)
    ftlI=(endpoint(ftl,False),endpoint(ftl,True));ftuI=(endpoint(ftu,False),endpoint(ftu,True))
    ck(0<ftlI[0]<=ftlI[1]<=ftuI[0]<=ftuI[1],'independent prefix intervals monotone')
    piL,piU=saved(fh['constants']['pi']);jl,ju=map(Q,cap['B_new_time_integral'])
    ck(3<piL<piU<Q(22,7),'saved062 pi endpoints')
    ck(jl==C0/piU**6*fl+p0/(piU**6*5*tu**5) and ju==C0/piL**6*fu+p0/(piL**6*5*tl**5),'full prefix+tail endpoint directions')
    independentJL=C0/piU**6*ftlI[0]+p0/(piU**6*5*tu**5)
    independentJU=C0/piL**6*ftuI[1]+p0/(piL**6*5*tl**5)
    ck(jl<=independentJL<independentJU<=ju and Q(4105,10000)<jl<ju<Q(4106,10000),'independent complete newJ enclosure inside saved bounds')
    stage='H_and_account'
    vp=saved(fh['constants']['V_power_5_over6']);ps,pl=saved(fh['Pshort']),saved(fh['Plong']);oldN=saved(fh['union_norm'])
    ck(h['Pshort_frozen']==list(map(str,ps)) and h['Plong_frozen']==list(map(str,pl)) and h['V_power_5_over6_frozen']==list(map(str,vp)),'full frozen P/volume norm inputs')
    bl,bu=map(Q,h['B_new_norm']);nl,nu=map(Q,h['new_H_norm']);il,iu=map(Q,h['new_H_I'])
    ck(bl==vp[0]*jl and bu==vp[1]*ju,'new B norm only')
    ck(nl==ps[0]+pl[0]+bl and nu==ps[1]+pl[1]+bu and 0<nl<=nu,'same-domain Minkowski norm')
    ck(0<il<=iu and il**5<=nl**6<=nu**6<=iu**5,'new H-I exact saved root bounds')
    oldI=Q(inp['ledger']['H_root']['upper'])
    ck(Q(h['strict_norm_improvement_residual'])==oldN[0]-nu>0 and Q(h['strict_I_improvement_residual'])==oldI-iu>0 and h['strict_improvement'],'strict formula and account improvement')
    for name in ['C_root','E_root','H_root']:
        z=inp['ledger'][name];ck(Q(z['lower'])**int(z['degree'])<=Q(z['radicand'])<=Q(z['upper'])**int(z['degree']),'original saved C/E/H exact root')
    ce=Q(inp['ledger']['C_root']['upper'])+Q(inp['ledger']['E_root']['upper']);old=Q(l['old_CEH_I_upper']);new=Q(l['new_CEH_I_upper']);delta=Q(l['I_recovery_exact'])
    ck(ce+oldI==old==Q(inp['ledger']['paid_CEH'])==Q(inp['paid_CEH_111']),'old exact partial account')
    ck(new==ce+min(oldI,iu) and Q(l['new_H_I_upper_used'])==min(oldI,iu) and delta==old-new>0,'only H replaced')
    ck(l['C_I_upper_unchanged']==inp['ledger']['C_root']['upper'] and l['E_I_upper_unchanged']==inp['ledger']['E_root']['upper'],'C/E unchanged exactly')
    rl,ru=map(Q,l['old_far_allocation']);newrl,newru=map(Q,l['new_far_sufficient_allocation'])
    ck(rl==Q(ref['S_plan_lower']) and ru==Q(ref['S_plan_upper']) and newrl==rl+delta and newru==ru+delta,'same exact recovery added to BOTH111 endpoints')
    ck(Q(25,10**8)<delta<Q(251,10**9) and Q(473,10**8)<new<Q(474,10**8),'strict coarse recovery and newCEH displays')
    ck(cert['new_H']['new_B_full_time_integral']==cap['B_new_time_integral'] and cert['new_H']['new_H_norm_interval']==h['new_H_norm'] and cert['new_H']['new_H_I_interval']==h['new_H_I'],'newH certificate copies')
    ck(cert['exact_CEH_account']['new_CEH_I_upper']==l['new_CEH_I_upper'] and cert['exact_CEH_account']['I_recovery_exact']==l['I_recovery_exact'],'account certificate copies')
    ck(cert['remaining_far_sufficient_allocation']['new']==l['new_far_sufficient_allocation'],'full far endpoint certificate copies')
    ck(not l['g0_globally_proved'] and not l['actual_far_paid'] and not l['global_S_closed'] and not out['new_C6_bound'],'partial-only scope')
    for p,digest in pins.items():ck(sha(Path(p))==digest,'source unchanged after independent audit')
    save('review_results.json',{'status':'PASS_ONE_SAVED_AND_NEW_PREFIX_VERIFIER_CONFIGURATION','assertions':count,'new_verifier_configurations':1,'source_paths_pre_post':len(pins),'author_source_pins':31,'original062_pins':10,'arb_precision_bits':640,'flint_version':getattr(flint,'__version__','unavailable'),'new_prefix_method':'atan half-angle, signed exact residues, saved T endpoints only','independent_prefix_TL_interval':list(map(str,ftlI)),'independent_prefix_TU_interval':list(map(str,ftuI)),'independent_Jnew_enclosure_inside_saved512':True,'exact_delta':str(delta),'exact_newCEH':str(new),'same_delta_added_to_both111_endpoints':True,'original_P_volume_tau_C_E_recomputed':False,'author_generator_old_function_peak_MC_native_grids_MATLAB_children':0,'partial_budget_only_no_new_Cn':True,'globalg0_far_S_unpaid':True,'source_hashes':pins,'wall_before_result_write':time.monotonic()-start,'CPU_before_result_write':time.process_time()-cpu})
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    save('run_receipt.json',{'status':status,'pid':os.getpid(),'assertions':count,'stage':stage,'failure':failure,'attempts':1,'automatic_retries':0,'children':0,'current_user_runtime_ceiling':None,'wall_before_final_receipt_write':time.monotonic()-start,'CPU_since_start_before_final_receipt_write':time.process_time()-cpu,'CPU_process_total_before_final_receipt_write':time.process_time(),'timing_scope':'one process before final receipt/stdout/exit; complete tool roundtrip recorded separately','script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json')})
    print(json.dumps({'status':status,'assertions':count,'pid':os.getpid()}))
if status!='COMPLETED':raise SystemExit(1)
