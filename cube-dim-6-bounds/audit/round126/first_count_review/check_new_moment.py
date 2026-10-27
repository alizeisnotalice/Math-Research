"""Independent126 new-U/confluent re-enclosure; no old fan/peak run."""
import sys,os
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS']:os.environ[key]='1'
import ast,datetime,hashlib,itertools,json,math,time,traceback
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
start,cpu=time.monotonic(),time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(p):return json.loads(p.read_text())
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
if any((P/n).exists() for n in ['started.json','review_results.json','run_receipt.json']):raise SystemExit('one attempt only')
reg=rd(P/'registration.json');assert sha(Path(__file__))==reg['script_sha256'];assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
save('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'registration_sha256':sha(P/'registration.json')})
checks=0;completed=0;functional_count=0;stage='source';failure=None;status='FIRST_FAILURE_STOP'
def ck(ok,label):
    global checks
    assert ok,label
    checks+=1
def safe(tree,b,h):
    if isinstance(tree,ast.Expression):return safe(tree.body,b,h)
    if isinstance(tree,ast.Constant) and isinstance(tree.value,int):return Q(tree.value)
    if isinstance(tree,ast.Name):return {'b':b,'h':h}[tree.id]
    if isinstance(tree,ast.UnaryOp):
        v=safe(tree.operand,b,h)
        if isinstance(tree.op,ast.USub):return -v
        if isinstance(tree.op,ast.UAdd):return v
    if isinstance(tree,ast.BinOp):
        l,r=safe(tree.left,b,h),safe(tree.right,b,h)
        if isinstance(tree.op,ast.Add):return l+r
        if isinstance(tree.op,ast.Sub):return l-r
        if isinstance(tree.op,ast.Mult):return l*r
        if isinstance(tree.op,ast.Div):return l/r
        if isinstance(tree.op,ast.Pow) and r.denominator==1:return l**int(r)
    raise ValueError('non-rational coefficient AST')
def rational_coeffs(m,b,h):
    K=m-1;dinv=[Q(2)**j/(1-2*b)**(j+1) for j in range(m)];aa=[];bb=[]
    for k in range(m):
        aa.append(sum(Q(math.comb(m,j))*b**(m-j)*dinv[k-j] for j in range(k+1)))
        bb.append(sum((Q(math.comb(m-1,j))*b**(m-1-j)-Q(math.comb(m,j))*b**(m-j))*dinv[k-j] for j in range(k+1)))
    inv2=[Q((-1)**j*math.comb(j+1,j))*b**(-2-j) for j in range(m)]
    inv3=[Q((-1)**j*math.comb(j+2,j))*b**(-3-j) for j in range(m)]
    pp=[-2*h*(u-v) for u,v in zip(inv3,inv2)];qq=[2*v for v in inv2]
    ta=[Q(1)];te=[Q(0)]
    for k in range(K):
        ta.append(sum(pp[j]*ta[k-j] for j in range(k+1))/(k+1))
        te.append((sum(pp[j]*te[k-j] for j in range(k+1))+qq[k])/(k+1))
    return [[sum(aa[j]*ta[k-j] for j in range(k+1)),bb[k],sum(aa[j]*te[k-j] for j in range(k+1)),-aa[k]] for k in range(m)]
def divided(nodes,jets,denom=lambda x:x):
    vals=[jets(z,0) for z in nodes]
    for order in range(1,len(nodes)):
        new=[]
        for j in range(len(nodes)-order):
            new.append(jets(nodes[j],order) if nodes[j]==nodes[j+order] else (vals[j+1]-vals[j])/denom(nodes[j+order]-nodes[j]))
        vals=new
    return vals[0]
try:
    pins={}
    for row in rd(P/'source_manifest.json')['files']:
        p=Path(row['path']);ck(sha(p)==row['sha256'] and p.stat().st_size==row['bytes'],'source byte pin');pins[str(p)]=sha(p)
    A=Path(reg['author_directory']);spec=rd(A/'inputs.json');out=rd(A/'results.json');old=rd(A/'old_caps_snapshot.json');nt=rd(A/'node_table.json');coeffile=rd(A/'rational_jet_coefficients.json')
    rel=rd(A/'primary_release.json');areg=rd(A/'registration.json')
    ck(out['status']=='COMPLETE_EIGHT_TRUE_AUXILIARY_MOMENT_INTERVALS' and out['exit_code']==0 and len(out['rows'])==8,'terminal8h complete')
    ck(rel['source_sha256']==areg['source_sha256']==out['source_sha256']==sha(A/'certify.py') and rel['input_sha256']==out['input_sha256']==sha(A/'inputs.json'),'source/release lineage')
    ck(out['node_table_sha256']==sha(A/'node_table.json') and out['coefficient_sha256']==sha(A/'rational_jet_coefficients.json'),'terminal derived tables')
    fan=rd(Path(reg['fan_cells']));audit=rd(Path(reg['fan_audit']));R0=Q(spec['old_full_R0'])
    ck(Q(fan['full_exp_minus_symbol_integral'])==Q(audit['full_exp_minus_symbol_integral'])==Q(old['R0'])==R0 and audit['status']=='PASS','old audited R0 read only')
    ck(sha(A/'fan_integer_input.txt')==spec['fan_sha256'],'integer fan bytes')
    tok=list(map(int,(A/'fan_integer_input.txt').read_text().split()));ck(tok[0]==21 and len(fan['cells'])==21,'existing21 source cones')
    V=[(1,)+x for x in itertools.product((-1,1),repeat=5)];cursor=1;alphas=[];pairnodes={}
    for cid,cell in enumerate(fan['cells']):
        det=tok[cursor];cursor+=1;rays=[];rates=[];scales=[]
        for j in range(6):
            rate=tok[cursor];ray=tok[cursor+1:cursor+7];cursor+=7;orig=list(map(Q,cell['rays_a0_equals_one'][j]));scale=Q(ray[0])/orig[0]
            ck(scale>0 and all(Q(a)==scale*b for a,b in zip(ray,orig)) and Q(rate)==scale*Q(cell['rates'][j]),'positive scaled source ray/rate')
            rays.append(ray);rates.append(rate);scales.append(scale)
        ck(det>0 and Q(det)==Q(cell['determinant'])*math.prod(scales),'saved determinant scaling, no fan recomputation')
        alphas.append(Q(det,math.prod(rates)))
        for vid,v in enumerate(V):
            dots=[sum(a*b for a,b in zip(v,r)) for r in rays]
            ck(min(dots)>=0 or max(dots)<=0,'new spectral-label weak signs')
            nodes=tuple(sorted(Q(abs(z),q) for z,q in zip(dots,rates) if z))
            ck(1<=len(nodes)<=6 and all(0<b<Q(1,2) for b in nodes),'positive spectral nodes and zero-variable deletion')
            pairnodes[(cid,vid)]=nodes
    ck(cursor==len(tok) and 46080*sum(alphas)==R0,'old source mass consistency only')
    tuples=nt['tuples'];seen=set();weighttotal=Q(0);repeat_count=0
    ck(len(tuples)==nt['unique_tuples']==out['unique_functional_tuples']==512,'all grouped new tuples')
    for tid,row in enumerate(tuples):
        nodes=tuple(map(Q,row['nodes']));ck(row['id']==tid and list(nodes)==sorted(nodes),'ordered tuple ids/nodes')
        weight=Q(0)
        for pair in row['source_cone_label_pairs']:
            key=tuple(pair);ck(key not in seen and pairnodes[key]==nodes,'all cone/label bindings exactly once');seen.add(key);weight+=46080*alphas[key[0]]
        ck(weight==Q(row['weight'])>0,'positive grouped weight')
        weighttotal+=weight
        if len(set(nodes))<len(nodes):repeat_count+=1
        m=len(nodes)
        ck(divided(nodes,lambda b,k:Q(math.comb(m-1,k))*b**(m-1-k))==1,'exact repeated-node polynomial normalization')
    ck(len(seen)==672 and weighttotal==32*R0,'all32 labels on all21 cones, weight normalization')
    asts={key:[ast.parse(s,mode='eval') for s in values] for key,values in coeffile['coefficients'].items()}
    ck(len(asts)==21 and coeffile['h0_and_hInfinity_symbolic_checks'],'rational coefficient schema')
    pirow=rd(Path(reg['pi_source']))['constants']['pi'];piL,piU=[Q(pirow[s+'_mantissa'])*Q(2)**int(pirow[s+'_binary_exponent']) for s in ['lower','upper']]
    ck(3<piL<piU<Q(22,7),'existing exact pi bounds')
    import flint
    from flint import arb,fmpq,ctx
    ctx.prec=768
    def B(q):q=Q(q);return arb(fmpq(q.numerator,q.denominator))
    def grid(value):
        vals=[];S=1<<512
        for side,x in enumerate([value.lower(),value.upper()]):
            m,e=x.man_exp();q=Q(int(m))*Q(2)**int(e);scaled=q*S
            z=scaled.__floor__() if side==0 else scaled.__ceil__();vals.append(Q(z,S))
        return vals
    previous=None;summary=[];ucount=0
    for index,hs in enumerate(spec['h_values']):
        stage='new_h_'+str(index);h=Q(hs);ah=B(h);eh=(-ah).exp();sqrtpi=arb.pi().sqrt();f=ah.sqrt().erf();E=eh*ah.sqrt()/sqrtpi;Tcache={};jetcache={}
        savedh=rd(A/('h_'+str(index)+'.json'));ck({k:v for k,v in savedh.items() if k!='functional_values'}==out['rows'][index],'h terminal summary copy')
        oldrow=old['rows'][index];source=Path(oldrow['source_path']);ck(sha(source)==oldrow['source_sha256'] and oldrow['old_B1z_expression_interval']==[rd(source)['caps']['b_lower'][1],rd(source)['caps']['beta_upper'][1]],'complete old cap source copy, no scalar rerun')
        def jets(m,b):
            key=(m,b)
            if key not in jetcache:
                if b not in Tcache:
                    c=1/B(b)-1;x=c*c*ah
                    Tu=eh*x.hypgeom_u(B(Q(1,2)),B(Q(1,2)))/sqrtpi
                    Te=((c*c-1)*ah).exp()*(c*ah.sqrt()).erfc()
                    ck(Tu.is_finite() and Te.is_finite() and 0<Te<1 and (Tu-Te).contains(0),'new U identity crosschecked by independent direct erfc')
                    Tcache[b]=Te
                rat=rational_coeffs(m,b,h);vals=[]
                for k,cs in enumerate(rat):
                    original=[safe(t,b,h)/math.factorial(k) for t in asts[str((m,k))]]
                    ck(cs==original,'independent ODE Taylor coefficients equal saved derivative/k!')
                    a,bf,ce,d=cs;vals.append(B(a)*Tcache[b]+B(bf)*f+B(ce)*E+B(d))
                jetcache[key]=vals
            return jetcache[key]
        JL=JU=Q(0);functions=[]
        ck(len(savedh['functional_values'])==512,'all512 new functional intervals')
        for row,sv in zip(tuples,savedh['functional_values']):
            nodes=list(map(Q,row['nodes']));m=len(nodes);value=divided(nodes,lambda b,k:jets(m,b)[k],B);lo,hi=grid(value);sl,su=map(Q,sv['interval'])
            ck(sv['tuple_id']==row['id'] and 0<sl<=lo<=hi<=su<1,'independent confluent new interval inside saved384 endpoints')
            wt=Q(row['weight']);JL+=wt*lo;JU+=wt*hi;functional_count+=1
            functions.append({'tuple_id':row['id'],'independent512_grid_interval':[str(lo),str(hi)]})
        ucount+=len(Tcache)
        sJL,sJU=map(Q,savedh['pi_cubed_EdN_interval']);ml,mu=map(Q,savedh['EdN_interval']);nl,nu=JL/piU**3,JU/piL**3
        ck(0<sJL<=JL<=JU<=sJU<32*R0 and ml<=nl<=nu<=mu,'independent exact new moment sums/division inside saved bounds')
        if previous is not None:ck(previous<JL,'positive threshold monotonicity across fixedh')
        previous=JU
        oldL,oldU=map(Q,oldrow['old_B1z_expression_interval'])
        ck(nu<oldL and mu<oldL and Q(savedh['old_lower_minus_new_upper'])==oldL-mu>0,'strict true-k1 cap improvement')
        ck(list(map(Q,savedh['old_cap_over_new_true_interval']))==[oldL/mu,oldU/ml],'positive complete comparison ratios')
        save('recheck_h_'+str(index)+'.json',{'h':hs,'method':'independent erfc/ODE Taylor Hermite at768bit,512bit outward endpoints','new_functional_intervals':functions,'independent_pi_cubed_EdN_interval':[str(JL),str(JU)],'independent_EdN_interval':[str(nl),str(nu)],'strict_old_cap_improvement':True,'source_h_sha256':sha(A/('h_'+str(index)+'.json'))})
        summary.append({'h':hs,'new_EdN_interval':savedh['EdN_interval'],'strict_improvement':True,'source_improvement_display':savedh['cap_improvement_factor_display']});completed+=1
        print(json.dumps({'completed_h':completed,'functional_checks':functional_count}),flush=True)
    ck(completed==8 and functional_count==4096 and out['MC_samples']==out['old_peak_producer_calls']==out['physical_point_calls']==0 and not out['new_C6_bound'] and not out['M_or_global_g0_certificate'],'finite first-count-only scope')
    for path,digest in pins.items():ck(sha(Path(path))==digest,'source unchanged after new-moment re-enclosure')
    save('review_results.json',{'status':'PASS_TRUE_FIRST_COUNT_AND_INDEPENDENT_NEW_U_CONFLUENT_REENCLOSURE','assertions':checks,'completed_h':8,'new_functional_intervals_reenclosed':4096,'grouped_tuples':512,'repeated_node_tuples':repeat_count,'ungrouped_pairs':672,'new_U_erfc_checks':ucount,'Arb_precision_bits':768,'outward_saved_grid_bits':512,'flint_version':getattr(flint,'__version__','unavailable'),'all_saved384_functionals_contain_independent_intervals':True,'all8_true_k1_caps_strictly_improve_old':True,'summary':summary,'source_paths_pre_post':len(pins),'source_hashes':pins,'old_fan_peak_124scalar_125optimizer_MC_native_MATLAB_children':0,'new_global_M_g0_S_Cn':False,'wall_before_result_write':time.monotonic()-start,'CPU_before_result_write':time.process_time()-cpu})
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    save('run_receipt.json',{'status':status,'pid':os.getpid(),'assertions':checks,'stage':stage,'completed_h':completed,'functional_checks':functional_count,'failure':failure,'attempts':1,'automatic_retries':0,'children':0,'science_CPU_cores':1,'current_user_runtime_ceiling':None,'wall_before_final_receipt_write':time.monotonic()-start,'CPU_since_start_before_final_receipt_write':time.process_time()-cpu,'CPU_process_total_before_final_receipt_write':time.process_time(),'timing_scope':'one new126 verifier before final receipt/stdout/exit; complete tool roundtrip separately observed','script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json')})
    print(json.dumps({'status':status,'assertions':checks,'completed_h':completed,'functional_checks':functional_count,'pid':os.getpid()}),flush=True)
if status!='COMPLETED':raise SystemExit(1)
