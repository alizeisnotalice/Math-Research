"""One saved-only finite117 check: exact geometry/interval algebra/bytes, no exp/native."""
import sys
sys.dont_write_bytecode=True
import collections,datetime,hashlib,importlib.util,json,os,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pairs(items):
    out={}
    for k,v in items:
        if k in out:raise ValueError('duplicate JSON key:'+k)
        out[k]=v
    return out
def read(p):return json.loads(Path(p).read_text(),object_pairs_hook=pairs,parse_constant=lambda s:(_ for _ in()).throw(ValueError(s)))
def load(n,p):
    s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
def require(x,s):
    if not x:raise ValueError(s)
def interval(x):
    a,b=map(Q,x);require(a<=b,'interval order');return a,b
def normalized_box(x):return [[str(Q(a)),str(Q(b))]for a,b in x]
def nk(c,w):return digest({'center':[str(Q(x))for x in c],'a':str(Q(w['a'])),'T':str(Q(w['T'])),'engineSHA':catalog['engineSHA']})
def ck(k,b):return digest({'native_key':k,'box':normalized_box(b),'target':'1/65000','quadratic':'frame32'})
def sk(keys,b):return digest({'keys':keys,'box':normalized_box(b)})
def resolve(r):
    x=read(r['path'])
    for k in r.get('selector',[]):x=x[k]
    return x
def scale(x,a):return tuple(sorted(a*v for v in x))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def new(p,x):
    with Path(p).open('x')as f:json.dump(x,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
require(not any((P/n).exists()for n in['saved_check_started.json','checked_saved.json','saved_check_receipt.json']), 'One saved check only')
t=time.monotonic();new(P/'saved_check_started.json',{'PID':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':sha(Path(__file__)),'scope':'saved exact/byte check only;0 exp/kernel/native'})
checks=[];status='FIRST_FAILURE_STOP';error=None
try:
    inp=read(P/'inputs.json');reg=read(P/'registration.json');catalog=read(P/'cache_catalog.json');ci=read(P/'cache_index.json');cc=read(P/'certificate_cache_index.json');forest=read(P/'forest.json');receipt=read(P/'run_receipt.json')
    fh=sha(P/'forest.json');reg_hash=sha(P/'registration.json')
    require(sha(P/'run_prepared.py')==reg['script_sha256']and sha(P/'cache_catalog.json')==reg['catalog_sha256']and sha(P/'source_manifest.json')==reg['source_manifest_sha256'],'frozen source hashes')
    sources=read(P/'source_manifest.json');require(sources['source_count']==len(sources['files']),'source count')
    for x in sources['files']:require(sha(x['path'])==x['sha256'],'source bytes:'+x['path'])
    rel=read(P/'parent_release.json');require(rel['authorized']is True and rel['script_sha256']==reg['script_sha256']and rel['registration_sha256']==reg_hash,'release binding')
    require(receipt['status']=='COMPLETED_96_FINITE_ROOTS'and receipt['failure']is None and receipt['release_sha256']==sha(P/'parent_release.json'),'completed producer receipt')
    checks.append('frozen117 release/script/catalog/source pins and completed producer receipt')
    base=load('saved117_geometry',inp['base_path']);M=load('saved117_dyadics',inp['math_contract_path']);V=load('saved117_validator',P/'native_validator.py');V.M=M
    cache=ci['entries'];seed={x['key']:x for x in catalog['native_entries']};require(len(seed)==183 and len(cache)==759,'183+576 point-window entries')
    data={};jets={};newkeys=[];pids=set();owner_counts=collections.Counter();window_counts=collections.Counter();native_rows=[]
    for k,v in cache.items():
        e=v['expected'];require(k==nk(e['center'],e['window'])and v['status']=='AVAILABLE','available exact cache key')
        require(e['model']==base.MODEL and e['model_definition']==base.MODEL_DEFINITION and e['window']in inp['windows'],'same law/window')
        r=read(v['receipt_path']);n=read(v['native_path']);i=read(v['intent_path']);pr=read(v['process_path'])
        require(i['input']==e and r['input_sha256']==digest(e),'native input identity')
        require(r['native_exit_code']==0 and Path(v['stderr_path']).read_bytes()==b''and r['native_sha256']==sha(v['native_path'])and r['stderr_sha256']==sha(v['stderr_path']),'native bytes/exit/stderr')
        enc=base.encode_center(e['center'],**inp['fan_bounds']);a=Q(e['window']['a']);T=Q(e['window']['T'])
        cmd=[inp['native_paths']['binary'],'full',inp['native_paths']['fan'],str(a.numerator),str(a.denominator),str(T.numerator),str(T.denominator),str(enc['denominator']),*map(str,enc['numerators'])]
        require(r['command']==cmd,'one exact-window argv')
        pp=r.get('source_binary_fan_runtime_hashes')or{inp['native_paths']['source']:r['source_sha256'],inp['native_paths']['binary']:r['binary_sha256'],inp['native_paths']['fan']:r['fan_sha256'],**r['runtime_sha256']}
        require(pp==reg['native_pins'],'native actual pins')
        d=V.validate_native(n,e);j=resolve(v['jet_ref']);base.validate_jet(j,e,test_only=False)
        require(interval(j['F_interval'])==d[0]and[interval(x)for x in j['face_gradient_intervals']]==d[1][1:],'matching whole native jet')
        data[k]=d;jets[k]=j
        if k in seed:
            require(all(v[a]==b for a,b in seed[k].items()),'old artifact unchanged')
        else:
            require(v['origin']=='117'and r['schema']=='C117-one-window-native-v1'and r['one_window_kernel_evaluations']==1,'new single-window schema')
            require(r['owned_native_pid']==pr['native_pid']==v['native_pid']and pr['wrapper_pid']==receipt['PID'],'new original PID binding')
            require(pr['intent_sha256']==sha(v['intent_path'])and r['intent_sha256']==sha(v['intent_path']),'intent bytes')
            require(r['release_sha256']==sha(P/'parent_release.json')and r['owned_native_pid']not in pids,'unique new PID/release')
            pids.add(r['owned_native_pid']);newkeys.append(k);owner_counts[v['owner_root']]+=1;window_counts[str(a)+'..'+str(T)]+=1
            native_rows.append({'key':k,'owner_root':v['owner_root'],'PID':v['native_pid'],'window':e['window'],'center':e['center'],'receipt_path':v['receipt_path'],'receipt_sha256':sha(v['receipt_path']),'native_path':v['native_path'],'native_sha256':sha(v['native_path']),'exit_code':0,'stderr_bytes':0,'wait_evidence':'producer native_exit_code is the return from direct Popen.wait;no independent ps assertion'})
    require(len(newkeys)==len(pids)==ci['new_started']==ci['new_reserved']==receipt['new_kernel_started']==576,'actual key/evaluation/PID count')
    require(len(list((P/'new_native').iterdir()))==576,'no extra native attempt folder')
    checks.append('759 unique physical point-windows;183 original artifacts preserved;576 unique one-window native exits0/empty stderr')
    certs={x['key']:x for x in cc['entries']};require(len(certs)==757,'183+574 certificates')
    for k,c in certs.items():
        require(k==ck(c['native_key'],c['expected']['face_box']),'certificate exact key')
        out=resolve(c['certificate_ref']);base.check_saved_box(jets[c['native_key']],c['expected'],out)
        require(len(out['corners'])==32 and min(Q(r['certificate_upper'])for r in out['corners'])==Q(out['tangent_bound_upper']),'32 rows/upper minimum')
        for row in out['corners']:require(row['passes_target']==(Q(row['certificate_lower'])>=Q(inp['target'])),'corner flag')
    require(cc['new_certificates']==574 and cc['old_exp_recomputed']==0 and len(list((P/'new_certificates').rglob('*.json')))==574,'new/old exp generation accounting')
    sums={sk(x['keys'],x['box']):x for x in catalog['bandset_certificates115']}
    for x in sums.values():require(x['interval']==read(x['source_path'])['separate_band_wholebox_expression_interval'],'115 cached expression source')
    checks.append('757 saved32-corner formula/geometry records checked with exact arithmetic;574 new/0 old exponential regenerations')
    oldforest=read(inp['source114_forest_path']);require(sha(inp['source114_forest_path'])==inp['source114_forest_sha256'],'114 authoritative forest')
    require(inp['inherited_passed_leaves']==['r000']+[n for s in oldforest['roots'].values()for n in s['passed_leaves']],'42 inherited passed ledger')
    source_by_id={n['id']:(rid,n)for rid,s in oldforest['roots'].items()for n in s['nodes'].values()if n['id']in s['residual_frontier']}
    require(len(source_by_id)==96 and set(forest['roots'])=={x['id']for x in inp['roots']},'96 exact continuation roots')
    allids=set();blocked=[];queued=[];passed=[];root_rows=[];partial_counts=collections.Counter();bandset_hits=0;nodes_processed=0;split_count=0
    for spec in inp['roots']:
        rid=spec['id'];s=forest['roots'][rid];nodes=s['nodes'];srcroot,src=source_by_id[spec['source114_id']]
        require(spec['source_node']==src and spec['source_root']==srcroot and s['source114_id']==spec['source114_id'],'root source mapping')
        root=nodes[rid];require(root['parent']is None and all(root[k]==src[k]for k in['raw_box','box','tightening','depth']),'same full root geometry')
        require(s['new_reserved']==owner_counts[rid]<=8 and s['nodes_processed']<=16 and s['inflight']is None,'root caps/null inflight')
        reached=set();stack=[rid];processed=0
        while stack:
            nid=stack.pop();require(nid not in reached and nid not in allids,'unique reachable geometry node');reached.add(nid);allids.add(nid);n=nodes[nid]
            tight=base.tighten(n['raw_box']);require(tight==n['tightening']and tight['box']==n['box'],'exact tightening')
            keys=n.get('partial_prefix',{}).get('keys',[])
            if keys:
                pp=n['partial_prefix'];require(len(keys)<=5 and len(set(keys))==len(keys),'finite prefix')
                ds=[data[k]for k in keys]
                for z,k in enumerate(keys):require(cache[k]['expected']['window']==inp['windows'][z]and list(map(Q,cache[k]['expected']['center']))==list(map(Q,base.midpoint(n['box']))),'same center/adjacent prefix')
                F=ds[0][0];g=ds[0][1];inter=[]
                for z in range(1,len(ds)):
                    ab=(max(ds[z][3][0],ds[z-1][2][0]),min(ds[z][3][1],ds[z-1][2][1]));require(ab[0]<=ab[1],'shared density')
                    inter.append(list(map(str,ab)));F=add(F,ds[z][0]);g=[add(a,b)for a,b in zip(g,ds[z][1])]
                outer={'kind':'B','a':inp['windows'][len(keys)-1]['a'],'T':'1/2'}
                terms=[scale(F,Q(5)),scale(ds[0][3],Q('1/2')),scale(ds[-1][2],-Q(outer['a']))]+[scale(v,Q(x))for v,x in zip(g,base.midpoint(n['box']))]
                E=tuple(sum(t[z]for t in terms)for z in[0,1]);require(E[0]<=0<=E[1],'prefix Euler0')
                require(interval(pp['F_interval'])==F and[interval(x)for x in pp['six_gradients']]==g and pp['outer_window']==outer and interval(pp['Euler_interval'])==E and pp['internal_density_intersections']==inter,'complete exact saved prefix')
            if n['status']in['passed','split']:
                require(keys and F[0]>Q(inp['target']),'strict center before box certificate');refs=n['certificate_refs']
                if refs['cache']=='REUSED115_BANDSET':
                    key=sk(keys,n['box']);require(key in sums and refs['source_path']==sums[key]['source_path'],'real115 reuse');bound=interval(sums[key]['interval']);bandset_hits+=1
                else:
                    require(refs['keys']==keys and len(refs['certificates'])==len(keys),'per-band certificate list');bound=(Q(0),Q(0))
                    for k,ref in zip(keys,refs['certificates']):
                        ckey=ck(k,n['box']);require(ref['key']==ckey and ref['path']==certs[ckey]['certificate_ref']['path']and ref.get('selector',[])==certs[ckey]['certificate_ref'].get('selector',[]),'same full box band ref')
                        c=resolve(certs[ckey]['certificate_ref']);bound=add(bound,(Q(c['tangent_bound_lower']),Q(c['tangent_bound_upper'])))
                require(interval(n['whole_box_interval'])==bound and(n['status']=='passed')==(bound[0]>=Q(inp['target'])),'correct summed box decision')
            if n['status']=='split':
                split_count+=1;sp=base.split_box(n['box']);require(n['split']==sp and n['children']==[nid+'0',nid+'1'],'closed deterministic split')
                for suffix,raw in[('0',sp['lower_raw_box']),('1',sp['upper_raw_box'])]:
                    c=nodes[nid+suffix];require(c['parent']==nid and c['depth']==n['depth']+1 and c['raw_box']==raw,'complete child geometry');stack.append(c['id'])
            else:require(n['children']==[]and n['status']in['passed','queued','blocked','empty'],'leaf represented')
            if n['status']=='blocked':
                require(n['stop_reason']=='NATIVE_BUDGET_PARTIAL_PREFIX'and s['new_reserved']==8,'actual blocked reason')
                require(len(keys)<5 and n['next_window']==inp['windows'][len(keys)],'next missing window retained')
                if keys:require(F[1]<Q(inp['target']),'only strict low prefix requests next band')
                partial_counts[len(keys)]+=1;blocked.append({'root_id':rid,'node':n,'completed_prefix_keys':keys,'next_window':n['next_window']})
            if n['status']=='queued':queued.append({'root_id':rid,'node':n})
            if n['status']=='passed':passed.append({'root_id':rid,'node':n})
            processed+=n['status']in['passed','split','blocked']
        require(reached==set(nodes),'no orphan/lost node');require(processed==s['nodes_processed'],'processednode ledger');nodes_processed+=processed
        q=sorted([n['id']for n in nodes.values()if n['status']=='queued'],key=lambda k:(nodes[k]['depth'],k));b=[n['id']for n in nodes.values()if n['status']=='blocked'];pa=[n['id']for n in nodes.values()if n['status']=='passed'];em=[n['id']for n in nodes.values()if n['status']=='empty']
        require(s['queue']==q and set(s['blocked_leaves'])==set(b)and set(s['passed_leaves'])==set(pa)and set(s['empty_leaves'])==set(em),'complete disjoint leaf ledgers')
        require(s['residual_frontier']==s['queue']+s['blocked_leaves']and s['root_complete']==(not q and not b),'root completion/residual')
        root_rows.append({'root_id':rid,'source114_id':spec['source114_id'],'status':s['status'],'new_native':s['new_reserved'],'processed_nodes':processed,'passed':len(pa),'queued':len(q),'blocked':len(b),'root_complete':s['root_complete']})
    require(len(allids)==894 and split_count==399 and nodes_processed==623 and len(passed)==170 and len(queued)==271 and len(blocked)==54,'actual finite tree counts')
    require(len(forest['residual_frontier'])==325 and set(forest['residual_frontier'])=={x['node']['id']for x in queued+blocked},'full325 residual list')
    require(forest['whole_root_complete']is False and forest['global_g_certificate']is False and forest['new_C6_bound']is False,'no global theorem promotion')
    require(sha(P/'forest.json')==fh,'original forest bytes unchanged by check')
    checks.append('96 original domains exactly matched;894 nodes/399 closed splits;212 total passed plus325 residual leaves cover original ordered face')
    checks.append('54 budget-only blocked leaves retained:52 empty prefixes,2 complete high prefixes;all next windows and759 cached artifacts retained')
    output={'status':'PASS_FINITE_SAVED_COVERAGE_CACHE_UNIQUE_KERNEL_CHECK','checks':checks,'source_count':len(sources['files']),'native_cache_total':759,'old_native_preserved':183,'new_native':576,'new_unique_PIDs':len(pids),'new_native_windows':dict(window_counts),'new_certificates':574,'new_corner_rows':574*32,'root_rows':root_rows,'complete_roots':sum(x['root_complete']for x in root_rows),'new_passed_leaves':170,'inherited_passed_leaves':42,'queued_residual':271,'blocked_residual':54,'residual_total':325,'blocked_prefix_counts':dict(partial_counts),'blocked':blocked,'queued':queued,'passed':passed,'native_rows':native_rows,'reused115_bandset_uses':bandset_hits,'original_forest_sha256':fh,'new_exp_kernel_native_or_baseline_evaluations':0,'independent_exp_enclosure_recomputation':False,'global_g_certificate':False,'proof_scope':'saved exact algebra/geometry/identity/byte checks;existing A057/engine proof and actual generated enclosures remain dependencies'}
    new(P/'checked_saved.json',output);status='COMPLETED'
except BaseException:error=traceback.format_exc()
finally:
    new(P/'saved_check_receipt.json',{'status':status,'PID':os.getpid(),'elapsed_seconds':time.monotonic()-t,'checks':checks,'failure':error,'new_native':0,'new_exp':0,'source_sha256':sha(Path(__file__))})
    print(json.dumps({'status':status,'PID':os.getpid(),'elapsed_seconds':time.monotonic()-t,'checks':len(checks),'new_native_or_exp':0}))
if status!='COMPLETED':raise SystemExit(1)
