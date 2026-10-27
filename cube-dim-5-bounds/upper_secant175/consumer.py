#!/usr/bin/env python3
"""175 finite two64bit fourth-root witnesses, exact secant lower and total gain.
No H/R/J, factor, mean, old endpoint power or old slope witness execution.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,gcd
import argparse,os,resource,subprocess,sys,time
from common import sha,read,put,lines,inputcheck
HERE=Path(__file__).resolve().parent
N=1<<160
B=1<<64
B4=1<<256

def setup(psha):
    assert sha(HERE/'plan.json')==psha;p=read(HERE/'plan.json')
    for name in ['consumer.py','common.py','sources.json','prepare.py','proof.md','selfcheck.py']:
        assert sha(HERE/name)==p['code_and_source_hashes'][name]
    x=p['preflight'];assert sha(x['path'])==x['sha256'];pre=read(x['path'])
    assert pre['status']=='all_exact_cached_fields_bound' and pre['new_root_slope_gain_DFS_cases']==0
    assert pre['strict_original112_active_chord_bound_by_source_certificate'] is True
    assert pre['sources_sha256']==p['code_and_source_hashes']['sources.json']
    r=read(Path(x['path']).parent/'registration.json')
    assert sha(Path(x['path']).parent/'registration.json')==p['preflight_registration_sha256']
    assert r['prepare_sha256']==p['code_and_source_hashes']['prepare.py'] and r['common_sha256']==p['code_and_source_hashes']['common.py']
    x=p['selfcheck'];assert sha(x['path'])==x['sha256'];sc=read(x['path'])
    assert sc['status']=='syntax_import_scope_and8new_small_root_boundaries_pass' and sc['actual_node_root_gain_DFS_cases']==0
    assert sha(HERE/'authorization.json')==p['authorization_sha256']
    ident={'plan_sha256':psha,'preflight_sha256':x['sha256'],'consumer_sha256':sha(__file__),'sources_sha256':pre['sources_sha256']}
    return p,pre,read(HERE/'sources.json'),ident

def inputs(src,pre):
    inputcheck(src)
    for x in pre['shards']+[pre['topology']]:assert sha(x['path'])==x['sha256']

def new_witness(s):
    L,U=s['range_L_i'],s['range_U_i'];assert 0<=L<U and U>=1<<144
    assert s['strict112_active'] is True and s['self112']<s['certificate112_old_selfU']
    Y=isqrt(isqrt((U*B4)//N));Z=isqrt(isqrt((L*B4)//U))
    yl=U*B4-Y**4*N;yr=(Y+1)**4*N-U*B4
    zl=L*B4-Z**4*U;zr=(Z+1)**4*U-L*B4
    assert yl>=0 and yr>0 and zl>=0 and zr>0 and 0<=Z<=B and Y>0
    # R(Z/B)=(B^4+ZB^3+Z^2B^2+Z^3B+Z^4)/(B*(B^3+ZB^2+Z^2B+Z^3)).
    D=B**3+Z*B**2+Z**2*B+Z**3
    rn=B*D+Z**4
    yn,yd=Y*rn,B*B*D
    g=gcd(yn,yd);yn//=g;yd//=g
    sn,sd=s['slope164_numerator'],s['slope164_denominator'];assert sn>0 and sd>0 and 16*sn>=sd
    use_new=yn*sd>sn*yd
    if use_new:sn,sd=yn,yd
    g=gcd(sn,sd);sn//=g;sd//=g;assert 16*sn>=sd
    delta=s['mean112_i']-s['mean162_i'];assert delta>=0
    num=s['weight']*delta*sn;den=sd<<(4*s['level']);total,rem=divmod(num,den)
    value=s['self112']-total;assert s['selfL']<=value<=s['self166']
    return {'key':s['key'],'Y':Y,'Z':Z,'Y_left_residual':yl,'Y_next_residual':yr,'Z_left_residual':zl,'Z_next_residual':zr,
        'new_secant_numerator':yn,'new_secant_denominator':yd,'selected_slope_numerator':sn,'selected_slope_denominator':sd,
        'selected_new_secant':use_new,'total_gain_relative112_units':total,'division_remainder':rem,'new_selfU':value,
        'baseline166_selfU':s['self166'],'strict_self_gain_units':s['self166']-value}

def worker(out,i,psha,resume=False):
    plan,pre,src,ident=setup(psha);assert read(out/'registration.json')['identity']==ident
    idx=pre['shards'][i];assert sha(idx['path'])==idx['sha256'];dest=out/f'worker{i:02d}'
    if not resume:
        dest.mkdir(exist_ok=False)
        put(dest/'registration.json',{'identity':ident,'worker':i,'pid':os.getpid(),'source_index':idx,'chunk_size':256})
        chunks=[];count=0;improved=0;selected=0
    else:
        state=read(dest/'state.json');assert state['identity']==ident and state['status']=='prefix'
        assert not (dest/'fault_prefix.json').exists(),'Uncommitted arithmetic prefix requires separate explicit disposition'
        chunks=state['chunks'];count=state['completed'];improved=state['self_improved_vs166'];selected=state['selected_new_secant']
        for c in chunks:assert sha(c['path'])==c['sha256']
    put(dest/f'ready{os.getpid()}.json',{'pid':os.getpid(),'worker':i})
    while not (out/'start.json').exists():time.sleep(.01)
    start=time.monotonic();buf=[];skipped=0;skip_initial=count if resume else 0
    def checkpoint():
        nonlocal buf
        if not buf:return
        f=dest/f'chunk{len(chunks):04d}.json';put(f,buf);chunks.append({'path':str(f),'sha256':sha(f),'count':len(buf)});buf=[]
        put(dest/'state.json',{'identity':ident,'completed':count,'chunks':chunks,'status':'prefix',
            'self_improved_vs166':improved,'selected_new_secant':selected},True)
    try:
        for s in lines(idx['path']):
            if skipped<skip_initial:skipped+=1;continue
            r=new_witness(s);buf.append(r);count+=1;improved+=r['strict_self_gain_units']>0;selected+=r['selected_new_secant']
            if len(buf)==256:checkpoint()
        checkpoint();assert count==idx['count'] and sha(idx['path'])==idx['sha256']
        usage=resource.getrusage(resource.RUSAGE_SELF)
        result={'identity':ident,'worker':i,'pid':os.getpid(),'status':'complete','completed':count,'chunks':chunks,
            'new_fourth_roots':2*count,'new_gains':count,'self_improved_vs166':improved,'selected_new_secant':selected,
            'seconds':time.monotonic()-start,'user_cpu_seconds':usage.ru_utime,'system_cpu_seconds':usage.ru_stime,'maxrss_native_units':usage.ru_maxrss}
        put(dest/'result.json',result);put(dest/'state.json',result,True)
    except BaseException as e:
        if buf:put(dest/'fault_prefix.json',buf)
        put(dest/'fault.json',{'pid':os.getpid(),'exception':repr(e),'complete_chunks':chunks,'action':'preserve;no automatic restart'});raise

def aggregate(out,pre,ident):
    tree=read(pre['topology']['path']);nodes={n['key']:n for n in tree['nodes']};assert len(nodes)==97972
    updates={};workers=[]
    for i in range(15):
        w=read(out/f'worker{i:02d}/result.json');assert w['identity']==ident and w['worker']==i and w['status']=='complete';workers.append(w)
        expected=iter(lines(pre['shards'][i]['path']));count=0
        for c in w['chunks']:
            assert sha(c['path'])==c['sha256'];rows=read(c['path']);assert len(rows)==c['count'];count+=len(rows)
            for r in rows:
                s=next(expected);k=r['key'];assert k==s['key'] and k not in updates
                assert r['baseline166_selfU']==s['self166']==nodes[k]['self166']
                assert s['self112']==nodes[k]['self112'] and s['selfL']==nodes[k]['selfL']
                assert s['selfL']<=r['new_selfU']<=s['self166'];updates[k]=r['new_selfU']
        assert next(expected,None) is None and count==w['completed']==pre['shards'][i]['count']
    assert len(updates)==97971 and set(nodes)-set(updates)=={pre['nonactive_key']}
    roots=tree['roots'];assert len(roots)==len(set(roots))==52360
    assert set(roots)=={k for k,n in nodes.items() if n['parent'] is None}
    for k,n in nodes.items():
        assert len(n['children'])==len(set(n['children']))
        for ch in n['children']:assert nodes[ch]['parent']==k and nodes[ch]['level']==n['level']+1
        if n['parent'] is not None:assert k in nodes[n['parent']]['children']
    eff={};rows=[];strict_eff=0;strict_roots=0
    for k in sorted(nodes,key=lambda k:(-nodes[k]['level'],k)):
        n=nodes[k];value=updates.get(k,n['self112']);children=n['children']
        assert n['selfL']<=value<=n['self166']
        e=min(value,sum(eff[c] for c in children)) if children else value
        assert e<=n['eff166'];eff[k]=e;strict_eff+=e<n['eff166']
        rows.append({'key':k,'new_selfU':value,'new_effU':e,'baseline166_selfU':n['self166'],'baseline166_effU':n['eff166'],
            'strict_self_gain_units':n['self166']-value,'strict_eff_gain_units':n['eff166']-e})
    units=sum(eff[k] for k in roots);base=F(tree['baseline166_UH'])*N
    assert base.denominator==1 and units<=base.numerator
    strict_roots=sum(eff[k]<nodes[k]['eff166'] for k in roots)
    put(out/'effective_nodes.json',rows)
    result={'task':'upper_secant175','status':'candidate_pending_independent_audit_primary','pid':os.getpid(),
        'workers':15,'worker_pids':[w['pid'] for w in workers],'new_fourth_roots':195942,'new_gains':97971,
        'active_nodes':97971,'nodes':97972,'roots':52360,'unchanged_self_nodes':1,'nonactive_key':pre['nonactive_key'],
        'self_improved_vs166':sum(w['self_improved_vs166'] for w in workers),'effective_improved_vs166':strict_eff,
        'roots_improved_vs166':strict_roots,'new_secant_selected':sum(w['selected_new_secant'] for w in workers),
        'new_UH':str(F(units,N)),'baseline166_UH':str(F(base.numerator,N)),
        'exact_gain_vs166':str(F(base.numerator-units,N)),'all_node_self_eff_monotone_vs166':True,
        'all_new_root_left_nonnegative_next_positive':True,'gain_baseline':'112 self minus one total floor relative112 mean;not166 self',
        'new_H_R_J_factor_mean_oldpower_old164witness_calls':0,'global_Cn_claim':False}
    put(out/'results.json',result)

def run(out,psha):
    plan,pre,src,ident=setup(psha);inputs(src,pre);out.mkdir(exist_ok=False);procs=[];start=time.monotonic()
    put(out/'registration.json',{'task':'175_unique_fixed_batch','identity':ident,'pid':os.getpid(),'workers':15,'new_fourth_roots':195942,
        'new_gains':97971,'checkpoints':390,'chunk_size':256,'wall_cap':None,'automatic_restart':False,'global_Cn_claim':False})
    try:
        for i in range(15):
            f=(out/f'worker{i:02d}.log').open('xb')
            p=subprocess.Popen([sys.executable,'-B',__file__,'worker','--output',str(out),'--plan-sha',psha,'--worker',str(i)],stdout=f,stderr=subprocess.STDOUT);f.close();procs.append(p)
        put(out/'worker_pids.json',{'pids':[p.pid for p in procs]})
        while not all(list((out/f'worker{i:02d}').glob('ready*.json')) for i in range(15)):
            if any(p.poll() is not None for p in procs):raise RuntimeError('worker exited before common start')
            time.sleep(.01)
        put(out/'start.json',{'pid':os.getpid(),'all15ready':True,'monotonic':time.monotonic()})
        samples=[];next_sample=0.
        while any(p.poll() is None for p in procs):
            if any(p.poll() not in (None,0) for p in procs):raise RuntimeError('worker failure;preserve prefix')
            now=time.monotonic()-start
            if now>=next_sample:
                ps=subprocess.run(['ps','-p',','.join(str(p.pid) for p in procs),'-o','pid=,ppid=,%cpu=,rss=,etime=,state='],capture_output=True,text=True)
                samples.append({'seconds':now,'ps_exitcode':ps.returncode,'output':ps.stdout});next_sample=now+.1
            time.sleep(.02)
        codes=[p.wait() for p in procs]
        for i,p in enumerate(procs):put(out/f'worker_exit{i:02d}.json',{'worker':i,'pid':p.pid,'exitcode':codes[i],'joined':True})
        put(out/'worker_exitcodes.json',{'pids':[p.pid for p in procs],'exitcodes':codes});assert codes==[0]*15
        put(out/'resource_samples.json',{'command':'ps pid,ppid,%cpu,rss,etime,state;only actual15 task workerPIDs','samples':samples})
        aggregate(out,pre,ident);inputs(src,pre);assert setup(psha)[3]==ident
        r=read(out/'results.json');r.update({'seconds':time.monotonic()-start,'all_inputs_before_after_unchanged':True,'registration_sha256':sha(out/'registration.json')})
        put(out/'completion.json',r);print(r)
    except BaseException as e:
        # Fault cleanup only; the normal computation has no time limit.
        events=[]
        for phase,delay in [('grace',.25),('terminate',2.),('kill',2.)]:
            for p in procs:
                if p.poll() is None and phase!='grace':
                    try:getattr(p,phase)();events.append({'pid':p.pid,'signal_phase':phase})
                    except ProcessLookupError:pass
            deadline=time.monotonic()+delay
            while any(p.poll() is None for p in procs) and time.monotonic()<deadline:time.sleep(.02)
        codes=[]
        for i,p in enumerate(procs):
            try:rc=p.wait(timeout=.01);joined=True
            except subprocess.TimeoutExpired:rc=p.poll();joined=False
            codes.append(rc);put(out/f'worker_exit{i:02d}.json',{'worker':i,'pid':p.pid,'exitcode':rc,'joined':joined,'fault_cleanup':True})
        put(out/'fault.json',{'pid':os.getpid(),'exception':repr(e),'exitcodes':codes,'cleanup_events':events,
            'prefix_preserved':True,'action':'stop;no automatic retry/newscope'});raise

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['run','worker']);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--plan-sha',required=True);ap.add_argument('--worker',type=int);ap.add_argument('--resume-authorized',action='store_true');a=ap.parse_args()
    if a.mode=='run':assert not a.resume_authorized;run(a.output,a.plan_sha)
    else:assert 0<=a.worker<15;worker(a.output,a.worker,a.plan_sha,a.resume_authorized)
