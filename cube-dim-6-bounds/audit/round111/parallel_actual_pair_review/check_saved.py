#!/usr/bin/env python3
"""Independent111 saved certificate audit; no author code import or cell replay."""
import sys
sys.dont_write_bytecode=True
import os,json,gzip,hashlib,time,subprocess,resource,traceback
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from datetime import datetime,timezone
OUT=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round111/parallel_actual_pair_geometry')
START=time.monotonic();CPU=time.process_time();HASHES={}
def need(x,msg):
    if not x:raise ValueError(msg)
def sha(p):
    p=Path(p);v=hashlib.sha256(p.read_bytes()).hexdigest();old=HASHES.setdefault(str(p),v);need(old==v,'Changed input '+str(p));return v
def read(p):sha(p);return json.loads(Path(p).read_text())
def gzread(p):sha(p);return json.loads(gzip.decompress(Path(p).read_bytes()))
def save(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
def iv(d):return Q(d['lower']),Q(d['upper'])
def dy(d):
    t=Q(2)**int(d['common_binary_exponent']);return int(d['lower_mantissa'])*t,int(d['upper_mantissa'])*t
def qr(p,a,b):
    pts={a,b}
    if p[2] and a < -p[1]/(2*p[2]) < b:pts.add(-p[1]/(2*p[2]))
    d={x:p[0]+p[1]*x+p[2]*x*x for x in pts};l=min(d.values());r=max(d.values())
    return l,r,{x for x,v in d.items() if v==l},{x for x,v in d.items() if v==r}
def above_gap(p,g,k,r):
    z=p-g*g-k*r*r
    return g>=0 and z>=0 and z*z>=4*g*g*k*r*r
def root_gap_upper(p,g,k,r):
    z=p-g*g-k*r*r
    return g>=0 and (z<=0 or z*z<=4*g*g*k*r*r)
def postcheck():
    for p,v in list(HASHES.items()):need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v,'Postcheck changed evidence')
def arithmetic():
    from flint import arb,fmpq,ctx
    ctx.prec=512
    def A(x):
        x=Q(x);return arb(fmpq(x.numerator,x.denominator))
    def ep(x):
        n,e=x.lower().man_exp();m,f=x.upper().man_exp();return Q(int(n))*Q(2)**int(e),Q(int(m))*Q(2)**int(f)
    def inside(x,bounds,msg):
        a,b=ep(x);l,r=bounds;need(l<=a<=b<=r,msg)
    return arb,A,ep,inside
def origin_poly(k,branch,a,b):
    m=6-k;cap=Q(1,2) if branch=='v' else Q(1,16)
    stars=[Q(k+m*u,6) if branch=='v' else Q(k-m*u,6) for u in (a,b)]
    if min(stars)>=cap:t=cap
    elif max(stars)<=-cap:t=-cap
    else:
        need(min(stars)>=-cap and max(stars)<=cap,'Cell crosses actual clamp')
        c=Q(k*m,6);return c,(-2 if branch=='v' else 2)*c,c
    return k*(1-t)**2+m*t*t,(-2 if branch=='v' else 2)*m*t,Q(m)
def breaks(k,a,b):
    m=6-k;points={a,b}
    for x in (Q(-3-k,m),Q(3-k,m),Q(k*8-3,8*m),Q(k*8+3,8*m),Q(1,2)):
        if a<x<b:points.add(x)
    return sorted(points)
def slab(jobid):
    reg=read(OUT/'registration.json');need(sha(Path(__file__))==reg['checker_sha256'],'Child exact source registration')
    need(sha(AUTHOR/'geometry_certificate.json')==reg['geometry_sha256'],'Child geometry pin')
    job=next(j for j in reg['jobs'] if j['id']==jobid);k=job['k'];m=6-k;delta=Q(1,64)
    index=read(AUTHOR/f'checkpoints/{jobid}/phase1-final/index.json')
    need(sha(AUTHOR/index['header']['path'])==index['header']['sha256'],'Header hash');h=gzread(AUTHOR/index['header']['path'])
    lineage=read(AUTHOR/'phase1_lineage.json')[jobid]
    need(h['job']==job and index['job_id']==jobid and index['tag']=='phase1-final','Exact slab job identity')
    need(h['geometry_sha256']==reg['geometry_sha256'] and h['transport_sha256']==index['transport_sha256']==lineage['transport_sha256'],'Geometry/phase1 transport binding')
    refinements=h['refinements'];prefix=hashlib.sha256(json.dumps(refinements,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    need(prefix==lineage['refinement_prefix_sha256'] and len(refinements)==lineage['refinements'],'Saved complete phase1 tree lineage')
    ua,ub=map(Q,job['u']);br=breaks(k,ua,ub)
    need(list(map(Q,h['initial_u_partition']))==br,'True clamp and line partition')
    expected={(a,b,s,t) for a,b in zip(br,br[1:]) for s,t in [(Q(0),Q(1,2)),(Q(1,2),Q(1))]}
    live={};seen=set()
    for seed in h['initial_seeds']:
        box=tuple(map(Q,seed['u']+seed['s']));cid=seed['id'];need(box in expected and cid not in seen,'Exact seed cover')
        expected.remove(box);live[cid]=(*box,0,0,cid);seen.add(cid)
    n0=job['initial_leaves'];need(not expected and len(live)==h['initial_leaves']==n0,'Full initial cover')
    for split in refinements:
        p=split['parent'];need(p in live,'Existing parent exactly once');a,b,s,t,du,ds,seed=live.pop(p);ids=split['children']
        need(len(ids)==2 and len(set(ids))==2 and all(i not in seen for i in ids),'Fresh two children')
        if split['axis']=='u':
            need(du<=ds,'Alternating normalized axis');mid=(a+b)/2;boxes=[(a,mid,s,t,du+1,ds,seed),(mid,b,s,t,du+1,ds,seed)]
        else:
            need(split['axis']=='s' and du>ds,'Alternating normalized axis');mid=(s+t)/2;boxes=[(a,b,s,mid,du,ds+1,seed),(a,b,mid,t,du,ds+1,seed)]
        for cid,box in zip(ids,boxes):live[cid]=box;seen.add(cid)
    leaves=[]
    need(h['leaf_chunks']==index['leaf_chunks'],'Complete leaf chunk set')
    for block in h['leaf_chunks']:
        p=AUTHOR/block['path'];need(sha(p)==block['sha256'] and p.stat().st_size==block['bytes'],'Exact compressed leaf bytes')
        leaves.extend(gzread(p)['leaves'])
    need(len(leaves)==len(live)==index['leaves']==h['leaf_count']==4096 and {x['id'] for x in leaves}==set(live),'All4096 saved leaves')
    need(len(seen)==index['cell_bounds']==h['cell_bounds']==lineage['cell_bounds']==8192-n0 and len(refinements)==4096-n0,'Phase1 finite cover counts')
    phase=h['completed_phases'];need(len(phase)==1 and phase[0]['phase']==1 and phase[0]['entry_counts']=={'leaves':n0,'cell_bounds':n0,'nextid':n0,'refinements':0} and phase[0]['continued_from_transport_sha256'] is None,'No refreshed/phase2 state')
    arb,A,ep,inside=arithmetic();pi=arb.pi();sqrtk=A(k).sqrt();rootkm=A(Q(k*m,6)).sqrt()
    ap=1/(2*A(k*m).sqrt());mass=A(Q(10,7))*pi*pi/2*A(Q(m,k)).sqrt()*A(delta)**A(Q(14,5))
    inside(ap,iv(h['large_shape_coefficient_1_over_Jvw']),'True pair coarea');inside(mass,iv(h['mass_constant']),'True sqrt(m/k) face/cone mass')
    ml,mh=iv(h['mass_constant']);geometry=read(AUTHOR/'geometry_certificate.json')['by_k'][str(k)]
    groups={b:[(tuple(map(Q,x['P'])),x['J2_squared'],x['n']) for x in geometry['branches'][b]['classes']] for b in ('v','w')}
    rootD={D:A(D).sqrt() for D in (20,32,36)};cache={};L=Q(0);U=Q(0);mingap=None
    for row in leaves:
        a,b,s,t,du,ds,seed=live[row['id']];rlo=delta*s**5;rhi=delta*t**5
        need(tuple(Q(row[x]) for x in ('a','b','s0','s1'))==(a,b,s,t) and row['u_depth']==du and row['s_depth']==ds and row['seed']==seed and du+ds<=28,'Exact leaf rectangle/depth/seed')
        g=Q(row['positive_gap_lower']);need(g>0,'Strict saved positive denominator');mingap=g if mingap is None else min(mingap,g)
        blo=A(0);bhi=A(0)
        for branch in ('v','w'):
            Wlo=A(0);Whi=A(0)
            for p,D,n in groups[branch]:
                key=(p,a,b)
                if key not in cache:
                    pl,pu,_,_=qr(p,a,b);cache[key]=(pl,pu,A(pl).sqrt(),A(pu).sqrt())
                pl,pu,rl,ru=cache[key];need(above_gap(pl,g,k,rhi),'Every saved pair gap by exact rational squaring')
                Wlo+=A(n)/(rootD[D]*(ru-sqrtk*A(rlo)));Whi+=A(n)/(rootD[D]*(rl-sqrtk*A(rhi)))
            fl,fu,_,_=qr(origin_poly(k,branch,a,b),a,b)
            if branch=='v':
                need(b<=Q(1,2) or a>=Q(1,2),'No line switch crossing');lo,hi=(a,b) if b<=Q(1,2) else (1-b,1-a)
            else:lo,hi=Q(7,8)-b,Q(7,8)-a
            need(above_gap(fl,g,k,rhi) and above_gap(Q(k*m,6)*lo*lo,g,k,rhi),'Every saved origin/line gap by exact squaring')
            h0lo=A(fl).sqrt()-sqrtk*A(rhi);h0hi=A(fu).sqrt()-sqrtk*A(rlo)
            h1lo=rootkm*A(lo)-sqrtk*A(rhi);h1hi=rootkm*A(hi)-sqrtk*A(rlo)
            blo+=Wlo/(h0hi*h0hi*h1hi*h1hi);bhi+=Whi/(h0lo*h0lo*h1lo*h1lo)
        tl,tu=Q(row['T_lower']),Q(row['T_upper'])
        need(0<tl<=ep(ap+A(rlo)*blo)[0] and ep(ap+A(rhi)*bhi)[1]<=tu,'Every saved full T interval contains independent512bit geometry')
        pl,pu=Q(row['power_lower']),Q(row['power_upper'])
        need(0<pl<=pu and pl**5<=tl**6 and pu**5>=tu**6,'All full rational fifth-power inequalities')
        weight=(b-a)*(t**14-s**14);need(Q(row['weight'])==weight>0,'Exact full s13 weight')
        l,u=Q(row['lower']),Q(row['upper']);need(l==ml*weight*pl and u==mh*weight*pu,'All full rational leaf contribution endpoints')
        L+=l;U+=u
    need((L,U)==(Q(h['lower']),Q(h['upper']))==(Q(index['shape_lower']),Q(index['shape_upper'])),'Exact4096leaf sums')
    postcheck();result=dict(status='PASS_ALL_SAVED_LEAF_CERTIFICATES',job_id=jobid,pid=os.getpid(),k=k,shape_lower=str(L),shape_upper=str(U),
        leaves=4096,tree_nodes=len(seen),all_rational_endpoints_exact=True,independent_geometry_bits=512,
        bitwise_replay_of_original384_operations=False,minimum_saved_gap=str(mingap),source_hashes=HASHES,
        wall_before_result=time.monotonic()-START,cpu_before_result=time.process_time()-CPU,author_replays=0,native_model_calls=0)
    save(f'slab_{jobid}.json',result);print(json.dumps({'status':'PASS','job':jobid,'pid':os.getpid()}))
def geometry_audit(geo):
    cases=0
    for k in (2,3,4,5):
        m=6-k;nA=k-1
        for branch in ('v','w'):
            table={}
            for a in range(nA+1):
                for b in range(m+1):
                    for c in range(nA+1):
                        for d in range(m+1):
                            omit=(nA,m if branch=='v' else 0)
                            if (a,b)==omit or (c,d)==omit:continue
                            for ta in range(max(0,a+c-nA),min(a,c)+1):
                                for tb in range(max(0,b+d-m),min(b,d)+1):
                                    h=a+c-2*ta+b+d-2*tb
                                    if not h:continue
                                    n=comb(nA,a)*comb(a,ta)*comb(nA-a,c-ta)*comb(m,b)*comb(b,tb)*comb(m-b,d-tb)
                                    ai=1+2*a-nA;bi=2*b-m;aj=1+2*c-nA;bj=2*d-m;dot=6-2*h;D=36-dot*dot
                                    p=(Q(k)-Q(6*ai*ai-2*dot*ai*aj+6*aj*aj,D),-Q(12*ai*bi-2*dot*(ai*bj+aj*bi)+12*aj*bj,D),Q(m)-Q(6*bi*bi-2*dot*bi*bj+6*bj*bj,D))
                                    table[p,D]=table.get((p,D),0)+n;cases+=1
            rows=geo['by_k'][str(k)]['branches'][branch]['classes'];saved={(tuple(map(Q,r['P'])),r['J2_squared']):r['n'] for r in rows}
            need(saved==table and len(rows)==len(table)==(23 if k==5 else 38) and sum(table.values())==930,'Independent two-block hypergeometric actual geometry')
            need({D:sum(n for (p,d),n in table.items() if d==D) for D in (20,32,36)}=={20:180,32:450,36:300},'All930 true coarea terms')
            for row in rows:
                lo,hi,lp,hp=qr(tuple(map(Q,row['P'])),Q(1,4),Q(3,4));gl,gu=dy(row['full_radius_positive_gap'])
                need((lo,hi)==(Q(row['P_min_exact']),Q(row['P_max_exact'])) and lp==set(map(Q,row['min_nodes'])) and hp==set(map(Q,row['max_nodes'])),'Exact full-u pair extrema')
                need(lo>Q(k,4096) and gl>0 and above_gap(lo,gl,k,Q(1,64)) and root_gap_upper(lo,gu,k,Q(1,64)),'Every full original384 root-gap endpoint by rational algebra')
    return cases
def main():
    reg=read(OUT/'registration.json');need(sha(Path(__file__))==reg['checker_sha256'],'Exact registered source')
    need(sha(AUTHOR/'manifest.json')==reg['author_manifest_sha256'],'Final author seal')
    for row in read(AUTHOR/'manifest.json')['files']:
        p=AUTHOR/row['name'];need(p.resolve().is_relative_to(AUTHOR.resolve()),'Input whitelist path')
        need(sha(p)==row['sha256'] and p.stat().st_size==row['bytes'],'Every final author artifact')
    for row in read(AUTHOR/'source_manifest.json')['sources']:need(sha(row['path'])==row['sha256'],'Frozen source/admin-asof bytes')
    source=read(AUTHOR/'registration.json');data=read(AUTHOR/'results.json');geo=read(AUTHOR/'geometry_certificate.json');up=read(AUTHOR/'upstream_constants.json')
    for name,digest in source['code_sha256'].items():need(sha(AUTHOR/name)==digest==data['code_sha256'][name],'Original parent/worker arithmetic source binding')
    qual=read(AUTHOR/'qualification_addendum.json');primary108=read(qual['primary108_path'])
    need(sha(qual['primary108_path'])==qual['primary108_sha256'] and primary108['status']=='PRIMARY_AND_INDEPENDENT_STRICT_AUXILIARY_THEOREM_ACCEPT' and primary108['Q3']=='1/(4096pi^3)' and primary108['independent_manifest_sha256']==qual['independent108_manifest_sha256'],'Formally accepted original-model108 input')
    need(source['jobs']==reg['jobs'] and data['phases_completed']==1 and data['k1_recomputed'] is False and data['new_C6_bound'] is False,'Fixed12 phase1 jobs; k1 readonly')
    for k in (2,3,4,5):need(sorted(tuple(map(Q,j['u'])) for j in reg['jobs'] if j['k']==k)==[(Q(1,4),Q(5,12)),(Q(5,12),Q(7,12)),(Q(7,12),Q(3,4))],'Three disjoint complete slabs')
    need(geo['new_pair_visits']==7440 and geo['k1_generated'] is False,'Original one-time geometry scope')
    cases=geometry_audit(geo)
    old107=read(AUTHOR.parents[1]/'round107/actual_pair_geometry/cell_certificate.json');need(up['k1_shape']==old107['shape_integral'],'Unchanged original107 shape interval')
    prior=OUT.parents[1]/'round107/actual_pair_review'
    need(sha(prior/'manifest.json')==reg['independent107_manifest_sha256'],'Frozen107 independent acceptance')
    oldchecked=read(prior/'checked_saved.json');need(oldchecked['status'].startswith('PASS_FULL107') and oldchecked['exact_shape_interval']==list(up['k1_shape'].values()),'Accepted107 full exact endpoints reused without verification replay')
    paid062=read('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round062/primary_joint_account/actual_run/results.json')
    need(Q(up['paid_CEH_I_upper'])==Q(paid062['paid_C_E_H_I_upper']),'Same actual fullG paidCEH')
    receipt=read(AUTHOR/'run_receipt.json');release=read(AUTHOR/'parent_release.json');tool=read(AUTHOR/'tool_observation.json')
    need(tool['initial_tool']['chunk_id']=='da275d' and tool['initial_tool']['session_id']==99408 and tool['terminal_tool']['chunk_id']=='d2325a' and tool['terminal_tool']['exit_code']==0 and tool['parent_PID']==receipt['pid']==60003 and tool['worker_PIDs']==receipt['owned_worker_pids'] and tool['phase2_invocations']==0,'Original actual session/PID/exit/one-phase binding')
    need(receipt['status']=='COMPLETED' and receipt['first_failure'] is None and receipt['retries']==0 and receipt['all_owned_workers_reaped'] and receipt['phases_completed']==1,'Original completion and owned reaps')
    need(all(x==0 for x in receipt['owned_worker_exit_codes'].values()) and len(receipt['owned_worker_pids'])==12 and all(not x for x in receipt['worker_stderr'].values()),'All12 original workers exit0')
    need(release['authorized'] and release['script_sha256']==sha(AUTHOR/'run_prepared.py') and release['registration_sha256']==sha(AUTHOR/'registration.json')==data['registration_sha256'] and receipt['release_sha256']==sha(AUTHOR/'parent_release.json'),'Original release binding')
    children={};outputs={};actions=[];before=resource.getrusage(resource.RUSAGE_CHILDREN);launch=time.monotonic()
    try:
        env=os.environ.copy()
        for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
        for job in reg['jobs']:children[job['id']]=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--job',job['id']],cwd=OUT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
        done=set();last=time.monotonic()
        while len(done)<12:
            for jid,p in children.items():
                if jid in done or p.poll() is None:continue
                stdout,stderr=p.communicate();outputs[jid]={'stdout':stdout.decode(),'stderr':stderr.decode(),'exit_code':p.returncode,'pid':p.pid};done.add(jid)
                need(p.returncode==0,'First audit slab failure '+jid+' '+stderr.decode())
            if time.monotonic()-last>=10:print(json.dumps({'stage':'saved_audit','done_slabs':len(done)}),flush=True);last=time.monotonic()
            if len(done)<12:time.sleep(.05)
    finally:
        for jid,p in children.items():
            if p.poll() is None:p.terminate();actions.append({'job':jid,'action':'terminate_owned_audit_child'})
        for jid,p in children.items():
            try:p.communicate(timeout=5)
            except subprocess.TimeoutExpired:p.kill();p.communicate(timeout=5);actions.append({'job':jid,'action':'kill_reap_owned_audit_child'})
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    L,U=iv(up['k1_shape']);sumleaves=0;sumnodes=0
    for k,n in ((2,5),(3,10),(4,10),(5,5)):
        kl=Q(0);ku=Q(0)
        for job in [j for j in reg['jobs'] if j['k']==k]:
            r=read(OUT/f'slab_{job["id"]}.json');need(r['status']=='PASS_ALL_SAVED_LEAF_CERTIFICATES','Slab result status')
            kl+=Q(r['shape_lower']);ku+=Q(r['shape_upper']);sumleaves+=r['leaves'];sumnodes+=r['tree_nodes']
        L+=n*kl;U+=n*ku
        target=data['final']['new_k2_to5'][str(k)];need(target['multiplicity']==n and (kl,ku)==(Q(target['shape_lower']),Q(target['shape_upper'])),'Exact multiplicity/slab sum')
    final=data['final'];need((L,U)==(Q(final['total_middle_shape_lower']),Q(final['total_middle_shape_upper'])) and sumleaves==49152 and sumnodes==98268,'Complete31middle exact rational endpoint account')
    arb,A,ep,inside=arithmetic();pi=arb.pi();q=A(Q(6,5));Kold=3*A(3).sqrt()/(64*pi**4);Kref=1/(8*pi**5)
    for x,key in ((Kold,'Kold'),(Kref,'Kref'),(Kold**q,'old_q'),(Kref**q,'ref_q')):inside(x,dy(data['scales'][key]),'Original coefficient scaling enclosure')
    paid=Q(up['paid_CEH_I_upper']);targets=[]
    for name,key in (('old098_valid','old_q'),('ref108_frozen_reference','ref_q')):
        sl,su=dy(data['scales'][key]);row=final['comparisons'][name];il,iu=sl*L,su*U
        need((il,iu)==(Q(row['middle_I_lower_expression']),Q(row['middle_I_upper_certificate'])),'Common coefficient and whole31middle')
        need((paid+il,paid+iu)==(Q(row['CEH_plus_middle_I_lower_expression']),Q(row['CEH_plus_middle_I_upper_certificate'])),'AllCEH/fullmiddle paid')
        tl=Q(row['S_plan_lower'])+paid+iu;tu=Q(row['S_plan_upper'])+paid+il;targets.append((tl,tu))
        inside(64*A(Q(1,390000))**q,(tl,tu),'Conditional g0 threshold')
        need(Q(row['S_plan_lower'])>0 and row['separation']=='STRICTLY_FITS_CONDITIONAL_PLAN' and row['g0']=='unproved globally' and 'UNPAID' in row['far'],'Remaining far allowance only')
    need(targets[0]==targets[1],'Identical target scalar for both comparisons')
    ref=final['comparisons']['ref108_frozen_reference'];need(Q(ref['S_plan_lower'])>=Q(4537105317164041237,10**24),'Exact conservative remaining far allowance')
    summary=read(AUTHOR/'exact_summary.json')
    def outward(item,l,u):need(Q(item['lower'])==Q(item['lower_decimal'])<=l<=u<=Q(item['upper_decimal'])==Q(item['upper']),'Every complete outward decimal/rational endpoint')
    outward(summary['whole_middle_shape'],L,U)
    for name,row in final['comparisons'].items():
        s=summary['coefficients'][name]
        for label,lk,uk in (('middle_I','middle_I_lower_expression','middle_I_upper_certificate'),('CEH_plus_middle_I','CEH_plus_middle_I_lower_expression','CEH_plus_middle_I_upper_certificate'),('remaining_far_allowance','S_plan_lower','S_plan_upper')):outward(s[label],Q(row[lk]),Q(row[uk]))
    need(summary['108_coefficient_accepted'] and summary['globalg0_and_far_unpaid'] and not summary['actual_I_lower_supplied'] and summary['results_sha256']==sha(AUTHOR/'results.json'),'Final strict scope')
    postcheck();result=dict(status='PASS_COMPLETE111_SAVED_EXACT_AND_ALL_LEAF_GEOMETRY',source_hashes=HASHES,hypergeometric_cases=cases,
        author_manifest_sha256=reg['author_manifest_sha256'],all_new_leaves_exact=49152,tree_nodes=98268,slabs=12,
        original384bit_endpoints_verified_by_exact_rational_arithmetic=True,independent_geometry_precision=512,
        original384operation_rounding_reproduced=False,whole_shape_lower=str(L),whole_shape_upper=str(U),
        remaining_far_allowance_lower=str(Q(ref['S_plan_lower'])),remaining_far_allowance_short_lower='4537105317164041237/10^24',
        accepted108_Kref=True,k1_recalculated=False,global_g_and_far_unpaid=True,complete_global_S=False,new_C6_bound=False,
        actual_native_model_calls=0,author_generator_or_phase1_replay=0,owned_audit_children=outputs,cleanup_actions=actions,
        all_own_children_waited=True,child_launch_to_all_reaps_wall=time.monotonic()-launch,
        reaped_audit_child_CPU=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
    save('checked_saved.json',result)
    save('run_receipt.json',dict(status='PASS',pid=os.getpid(),utc=datetime.now(timezone.utc).isoformat(),wall_before_receipt=time.monotonic()-START,
        parent_CPU_before_receipt=time.process_time()-CPU,owned_worker_pids={j:p.pid for j,p in children.items()},owned_worker_exits={j:p.returncode for j,p in children.items()},
        all_owned_children_reaped=True,child_CPU_after_reap=result['reaped_audit_child_CPU'],runtime_limit=None,attempts=1,
        author_phase_replays=0,native_model_calls=0,result_sha256=sha(OUT/'checked_saved.json')))
    print(json.dumps({'status':result['status'],'leaves':49152,'slabs':12,'remaining_far_lower':result['remaining_far_allowance_short_lower'],'complete_global_S':False}))
if __name__=='__main__':
    try:
        if len(sys.argv)==3 and sys.argv[1]=='--job':slab(sys.argv[2])
        else:main()
    except BaseException as error:
        suffix='_'+sys.argv[2] if len(sys.argv)==3 else ''
        save('failure'+suffix+'.json',dict(status='STOP_NO_RETRY',pid=os.getpid(),error=repr(error),traceback=traceback.format_exc(),wall=time.monotonic()-START,cpu=time.process_time()-CPU,author_replays=0))
        raise
