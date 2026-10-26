#!/usr/bin/env python3
"""Frozen089 pre-split, cache-aware journal, target105 certificate/tree audit."""
import argparse,json,gzip,os,time,resource
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import lcm
from hashlib import sha256
from interfaces import REG,HERE,common as C,arith as A


def write(path,obj):
    tmp=Path(str(path)+'.tmp');tmp.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n');tmp.replace(path)


def geom(point,vertices):
    p=list(map(F,point));D=lcm(*(x.denominator for x in p));nums=[int(x*D) for x in p];k=[sum(a*b for a,b in zip(nums,v)) for v in vertices]
    return {'eligible':0 not in k and len(set(k))==280,'zero_count':k.count(0),'distinct_count':len(set(k)),'common_denominator':str(D)}


def choose(n,roots,vertices):
    lo,hi=C.physical_box(roots,n);mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096
    candidates=[mid]+[[x+eta*beta**j for j,x in enumerate(mid)] for beta in (F(1,8),F(1,16),F(1,32))];tests=[]
    for i,tail in enumerate(candidates):
        p=['1',*map(str,tail)];test={'candidate':i,**geom(p,vertices)};tests.append(test)
        if test['eligible']:return p,i,tests
    raise AssertionError('no eligible generic anchor')


def certify(n,cert,record,roots):
    C.require(cert['mechanism']=='full_frame_16' and cert['quadratic']=='Q_full=I4/16' and cert['checked_points']==16,'wrong new-target fullQ schema')
    C.require(cert['anchor']==','.join(record['point']),'certificate/record anchor differs')
    values=[[F(a),F(b)] for a,b in record['values']];C.require(len(values)==5 and all(a<=b for a,b in values) and values[0][0]>F(21,200000),'value interval not eligible for g105')
    for a,b in values:C.require((b-a)*2**40<=max(abs(a),abs(b),F(1,10**6)),'actual numeric width condition fails')
    gs=[[F(a),F(b)] for a,b in cert['loggradient_directed']];C.require(len(gs)==4 and all(a<=b for a,b in gs),'unordered gradient')
    for (gl,gu),(dl,du) in zip(gs,values[1:]):
        for h in values[0]:
            for d in (dl,du):C.require(gl*h<=d<=gu*h,'physical gradient/H interval fails')
    delta=F(cert['Delta_lower']);ratio=values[0][0]/F(21,200000);lo_log,hi_log=A.log_bounds(ratio.numerator,ratio.denominator)
    C.require(delta>0 and delta<=F(lo_log,A.S),'Delta is not a lower bound for fresh g105 log')
    lo,hi=C.physical_box(roots,n);point=list(map(F,record['point']));C.require(point[0]==1 and all(a<=x<=b for a,x,b in zip(lo,point[1:],hi)),'anchor outside independently mapped closed box')
    offsets=[(a-x,b-x) for a,x,b in zip(lo,point[1:],hi)];vd=lcm(*(x.denominator for pair in offsets for x in pair));gd=lcm(delta.denominator,*(x.denominator for pair in gs for x in pair));di=int(delta*gd)
    zz=[tuple(int(x*vd) for x in pair) for pair in offsets];gg=[tuple(int(x*gd) for x in pair) for pair in gs];rows=[]
    for bits in product((0,1),repeat=4):
        z=[zz[j][bits[j]] for j in range(4)];num=di*vd*vd+vd*sum(min(a*x,b*x) for (a,b),x in zip(gg,z))-256*gd*sum(x*x for x in z);rows.append((num,bits))
    minimum=min(x[0] for x in rows);margin=F(minimum,gd*vd*vd);C.require(cert['passed']==(minimum>0) and F(cert['minimum_log_margin'])==margin and (minimum,tuple(cert['worst_corner_bits'])) in rows,'new g105 exact16corner differs')
    return str(margin)


def load_base(wanted):
    r=C.Stream(REG['source086']);roots=None;passed=[];splits=[];pending=[];cache={};keys=set();metadata={};h=sha256();h.update(b'{');first=True
    while True:
        field=r.key()
        if field is None:break
        if field=='cache':
            for key,record in r.items(True):
                keys.add(key)
                if not first:h.update(b',')
                first=False;h.update(C.canonical(key));h.update(b':');h.update(C.canonical(record))
                if key in wanted:cache[key]=record
            h.update(b'}')
        elif field=='passed':passed=[C.node(x['node']) for x in r.items()]
        elif field=='splits':
            for x in r.items():
                n=C.node(x['node']);C.require(x['axis']==len(n[1])%4 and x['children']==[[n[0],n[1]+'0'],[n[0],n[1]+'1']],'source split physical geometry differs');splits.append(n)
        elif field=='roots':roots=r.value()
        elif field=='pending_certificates086':pending=r.value()
        else:metadata[field]=r.value()
    r.close();C.require(len(keys)==180040 and len(passed)==96104 and len(splits)==45600 and len(pending)==1856,'086 source finite counts differ')
    C.require(metadata['target']=='21/200000' and metadata['epsilon']=='1/64' and metadata['model']=='alpha1','wrong source model/target/epsilon')
    tree=C.validate_tree(roots,passed+[C.node(p['node']) for p in pending],set(splits),'qualified086source')
    return {'roots':roots,'passed':passed,'splits':splits,'pending':pending,'cache':cache,'keys':keys,'cache_sha256':h.hexdigest(),'tree':tree}


def numeric_join(batch,point_manifest):
    d=HERE/('numeric_batch%04d'%batch);index=json.loads((d/'index.json').read_text());ir=json.loads((d/'input_registration.json').read_text());C.require(index['status']=='numeric_complete_all_new_H_plus4grad','numeric batch incomplete')
    C.require(C.file_sha(ir['point_manifest'])==ir['point_manifest_sha256'],'numeric source manifest changed')
    jobs=point_manifest['points'];proofs=[]
    for plan in index['blocks']:
        state=json.loads((d/('block%04d.json'%plan['block'])).read_text());C.require(state['status']=='complete' and len(state['points'])==plan['stop_exclusive']-plan['start'],'numeric block prefix incomplete')
        for path,h in state['code_locks'].items():C.require(C.file_sha(path)==h,'numeric successful-prefix code/input changed')
        for offset,proof in enumerate(state['points']):
            i=plan['start']+offset;job=jobs[i];C.require(proof['index']==i and proof['point']==job['point'] and proof['raw_checkpoint_sha256']==job['sha256'],'numeric/source point identity differs')
            for interval,production in zip(proof['independent_values_512_dyadic'],job['record']['values']):
                a,b=map(F,interval);L,U=map(F,production);C.require(L<a<=b<U,'independent value is not strictly inside production ball')
            C.require(len(proof['independent_values_512_dyadic'])==5,'not all5 components independently verified');proofs.append([i,proof['point'],proof['independent_values_512_dyadic']])
    C.require(len(proofs)==len(jobs)==index['verified_points'] and all(x['exitcode']==0 and x['is_alive'] is False for x in index['pool_exits']),'numeric finite count/pool exits differ')
    return {'points':len(proofs),'components':5*len(proofs),'results_sha256':C.digest(proofs),'master_pid':index['master_pid'],'workers':[x['pid'] for x in index['pool_exits']],'wall_seconds':index['wall_seconds'],'index_sha256':C.file_sha(d/'index.json')}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--registration',default='structural_through0004_registration.json');parser.add_argument('--resume');args=parser.parse_args();began=time.monotonic();pid=os.getpid()
    registration=HERE/args.registration;snap=json.loads(registration.read_text());out=HERE/(args.registration.replace('_registration.json','_result.json'))
    C.require(C.file_sha(snap['repack086_result'])==snap['repack086_result_sha256'],'fresh086 qualification receipt changed');rep=json.loads(Path(snap['repack086_result']).read_text());C.require(rep['status']=='086_new_target_all97960_qualified_partial' and rep['passed']==96104 and rep['pending']==1856,'086 has not qualified current target')
    inputs={REG['source086']:REG['source086_sha256'],REG['manifest086']:REG['manifest086_sha256'],str(registration):C.file_sha(registration),snap['repack086_result']:snap['repack086_result_sha256'],REG['fan']:REG['fan_sha256']}
    bundles=[];wanted=set()
    for item in snap['frozen_batches']:
        C.require(C.file_sha(item['manifest'])==item['manifest_sha256'],'batch manifest changed');mf=json.loads(Path(item['manifest']).read_text());inputs[item['manifest']]=item['manifest_sha256'];inputs.update(mf['files'])
        root=Path(item['manifest']).parent;events=[json.loads(line) for line in (root/'journal.jsonl').read_text().splitlines()]
        for e in events:
            if e['kind']=='dispatch' and e['body']['cached']:wanted.add(','.join(e['body']['point']))
        bundles.append({'number':item['batch'],'root':root,'manifest':mf,'events':events})
    pending_order=Path(REG['source086']).parent.parent/'pending_order089.json';m086=json.loads(Path(REG['manifest086']).read_text());inputs[str(pending_order)]=m086['external_files'][str(pending_order)]
    for row in json.loads(pending_order.read_text())['records']:wanted.add(row['anchor'])
    for path,h in {**inputs,**REG['protected_old_locks']}.items():C.require(C.file_sha(path)==h,'frozen/protected input changed: '+path)
    fan=json.loads(gzip.decompress(Path(REG['fan']).read_bytes()));base=load_base(wanted)
    result={'task':'audit_g105090_089_structure','pid':pid,'status':'verified_finite_prefix_only','target':'21/200000','epsilon':'1/64','inputs':inputs,'code_sha256':C.file_sha(Path(__file__)),'snapshot_registration_sha256':C.file_sha(registration),'batches':[],'certificates':[],'source_tree':base['tree'],'old180040_recomputed':False,'global_g105_qualified':False,'Cn_C5_I_updated':False}
    reused={}
    if args.resume:
        prior=json.loads(Path(args.resume).read_text());C.require(prior['code_sha256']==result['code_sha256'] and prior['target']==result['target'],'resume implementation/target differs')
        for proof in prior['certificates']:reused[(tuple(proof['node']),proof['certificate_sha256'],proof['record_sha256'])]=proof
        result['reused_certificate_receipt_sha256']=C.file_sha(args.resume)
    previous=None;allnew=set();allpids=set();total_events=0
    for bundle in bundles:
        root=bundle['root'];num=bundle['number'];br=json.loads((root/'registration.json').read_text());initial=json.loads((root/'initial_overlay.json').read_text());final=json.loads((root/'state.json').read_text());mf=bundle['manifest']
        C.require(br['source_state_sha256']==REG['source086_sha256'] and initial['source_sha256']==REG['source086_sha256'] and br['selected_nodes']==initial['queue'][:4096] and C.digest(br['selected_nodes'])==br['selected_sha256'],'batch source/selection differs')
        C.require(br['base_cache_canonical_sha256']==base['cache_sha256'] and br['wall_deadline'] is None and br['workers']==12,'batch cache/resource qualification differs')
        if previous is None:
            pre=[];children=[];hard={};mapper=json.loads((root.parent/'pre_split_mapper.json').read_text());pr=json.loads((root.parent/'pre_split_registration.json').read_text())
            C.require(pr['source_sha256']==REG['source086_sha256'] and pr['source_pending_certificates_sha256']==C.digest(base['pending']) and pr['target']=='21/200000' and pr['epsilon']=='1/64' and pr['new_H_gradient_kernel_log_calls']==0,'pre-split source/target/math scope differs')
            for i,pending in enumerate(base['pending']):
                n=C.node(pending['node']);cert=pending['certificate'];C.require(cert['passed'] is False,'pre-split consumes a passed source parent')
                kids=[[n[0],n[1]+'0'],[n[0],n[1]+'1']];pre.append({'node':list(n),'axis':len(n[1])%4,'children':kids,'failed_certificate':cert,'source086_pending_certificate':True});children.extend(kids)
                expected_map={'children':kids,'failed_certificate_sha256':C.digest(cert),'source070_leaf_index':pending['source070_leaf_index'],'source086_pending_node':list(n)};C.require(mapper[i]==expected_map,'pre-split mapper/certificate identity differs')
                record=base['cache'][cert['anchor']]
                for child in kids:hard[str(child[0])+':'+child[1]]=[record['values'][0][0],cert['minimum_log_margin']]
            children.sort(key=lambda n:(*map(F,hard[str(n[0])+':'+n[1]]),*n))
            C.require(len(pre)==len(mapper)==1856 and initial['new_splits']==pre and initial['queue']==children and len(children)==3712 and initial['hardness']==hard,'zero-H presplit is not exact closed bisection of all1856failed parents')
            C.require(initial['new_cache']=={} and initial['new_passed']==[] and initial['pre_split_count']==1856 and initial['counts']=={'cache_hits':0,'kernel_segment_calls':0,'new_anchor_starts':0,'new_complete':0},'presplit incorrectly creates new values or passed leaves')
            result['pre_split_tree']=C.validate_tree(base['roots'],base['passed']+list(map(C.node,children)),set(base['splits'])|{C.node(x['node']) for x in pre},'zeroHpresplit')
            result['zero_H_presplit_failed_parents']=1856;result['zero_H_presplit_closed_children']=3712
        else:
            C.require(br['parent_overlay_sha256']==C.file_sha(bundles[num-2]['root']/'state.json'),'overlay parent SHA differs')
            for key in ['new_cache','new_passed','new_splits','queue','hardness','counts','pre_split_count']:C.require(initial[key]==previous[key],'cumulative overlay prefix changed: '+key)
        expected=json.loads(json.dumps(initial));queue=list(initial['queue']);hard=dict(initial['hardness']);dispatch={};received={};merged=[];prefix='0'*64;seq=0;newcount=hits=calls=0;available=base['keys']|set(initial['new_cache'])
        pm=json.loads((root/'point_manifest.json').read_text());points={p['path']:p for p in pm['points']};uncached_paths=set();proofs=[]
        for event in bundle['events']:
            body=event['body'];kind=event['kind'];payload={k:v for k,v in event.items() if k!='hash'};C.require(event['seq']==seq+1 and event['previous']==prefix and C.digest(payload)==event['hash'],'journal hash chain differs');seq+=1;prefix=event['hash']
            if kind=='dispatch':
                ordinal=body['ordinal'];C.require(ordinal==len(dispatch) and body['node']==br['selected_nodes'][ordinal] and queue.pop(0)==body['node'],'dispatch does not match finite frontier')
                point,candidate,tests=choose(C.node(body['node']),base['roots'],fan['vertices']);C.require(body['point']==point and body['candidate']==candidate and body['geometry_tests']==tests,'generic anchor/guard differs')
                key=','.join(point);C.require(body['cached']==(key in available),'cache hit/new dispatch classification differs');dispatch[ordinal]=body
            elif kind=='received':
                ordinal=body['ordinal'];got=body['result'];task=dispatch[ordinal];C.require(ordinal not in received and got['kind']=='complete','received fault/counter/repeated ordinal')
                record=got['record'];key=','.join(task['point']);C.require(record['point']==task['point'],'returned point differs')
                if task['cached']:
                    old=expected['new_cache'].get(key,base['cache'].get(key));C.require(old is not None and record==old and got['calls']==0 and got['cache_reused'] is True,'cache reuse changed accepted source value or evaluated it again');hits+=1
                else:
                    entry=points[task['point_checkpoint']];C.require(C.file_sha(entry['path'])==entry['sha256'],'raw point SHA differs');raw=json.loads(Path(entry['path']).read_text())
                    C.require(raw['stage']=='complete' and raw['partial'] is None and raw['record']==record==entry['record'] and raw['task_sha256']==C.digest(task) and raw['point']==task['point'],'raw/task/record association differs')
                    C.require(raw['pid']==got['pid']==entry['pid'] and raw['kernel_segment_calls']==got['calls'] and record['candidate']==task['candidate'] and record['geometry_tests']==task['geometry_tests'],'new point precision/worker/geometry record differs')
                    C.require(key not in base['keys'] and key not in allnew and key not in expected['new_cache'],'old/successful anchor recomputed');allnew.add(key);expected['new_cache'][key]=record;newcount+=1;uncached_paths.add(entry['path'])
                calls+=got['calls'];received[ordinal]=got
            elif kind=='merge':
                ordinal=body['ordinal'];C.require(ordinal==len(merged) and ordinal in received,'merge order differs');task=dispatch[ordinal];record=received[ordinal]['record'];n=C.node(task['node']);cert=body['certificate'];certsha=C.digest(cert);recsha=C.digest(record);identity=(n,certsha,recsha)
                if identity in reused:proof=dict(reused[identity]);proof['reused_without_recompute']=True
                else:proof={'node':list(n),'certificate_sha256':certsha,'record_sha256':recsha,'minimum_log_margin':certify(n,cert,record,base['roots']),'passed':cert['passed'],'cached':task['cached']}
                proofs.append(proof);result['certificates'].append(proof);merged.append(ordinal);available.add(','.join(record['point']));hard.pop(str(n[0])+':'+n[1],None)
                if cert['passed']:expected['new_passed'].append({'node':list(n),'certificate':cert})
                else:
                    kids=[[n[0],n[1]+'0'],[n[0],n[1]+'1']];expected['new_splits'].append({'node':list(n),'axis':len(n[1])%4,'children':kids,'failed_certificate':cert});queue.extend(kids)
                    for child in kids:hard[str(child[0])+':'+child[1]]=[record['values'][0][0],cert['minimum_log_margin']]
            elif kind=='finish':
                C.require(len(merged)==len(dispatch)==len(received)==len(br['selected_nodes']),'finish has incomplete task sets');pool=body['pool_evidence'];C.require(pool['close_called'] is True and pool['join_returned'] is True and pool['master_pid']==br['master_pid'] and len(pool['worker_pids'])==12 and all(p['exitcode']==0 and p['is_alive'] is False for p in pool['exits']),'pool did not exit cleanly')
                C.require({g['pid'] for i,g in received.items() if not dispatch[i]['cached']}<=set(pool['worker_pids']),'actual new point worker not in pool');allpids.update(pool['worker_pids']);allpids.add(pool['master_pid']);expected['pool_evidence']=pool;expected['status']=body['status']
            else:raise AssertionError('unknown journal event kind')
        C.require(uncached_paths==set(points) and newcount==pm['new_points'],'point manifest/new record set differs');queue.sort(key=lambda n:(*map(F,hard[str(n[0])+':'+n[1]]),*n))
        expected.update({'queue':queue,'hardness':hard,'inflight':{},'resume_jobs':{},'next_merge':len(merged),'journal_seq':seq,'journal_hash':prefix,'counts':{'cache_hits':initial['counts']['cache_hits']+hits,'kernel_segment_calls':initial['counts']['kernel_segment_calls']+calls,'new_anchor_starts':initial['counts']['new_anchor_starts']+newcount,'new_complete':initial['counts']['new_complete']+newcount}})
        C.require(expected==final,'independently reconstructed overlay differs');C.require(final['counterexample'] is None and final['faults']==[] and not final['inflight'] and not final['resume_jobs'],'counter/fault/active work remains')
        tree=C.validate_tree(base['roots'],base['passed']+[C.node(r['node']) for r in final['new_passed']]+list(map(C.node,queue)),set(base['splits'])|{C.node(r['node']) for r in final['new_splits']},'089batch'+str(num))
        numeric=numeric_join(num,pm);C.require(numeric['points']==newcount,'independent numeric set does not match new dispatch set')
        result['batches'].append({'batch':num,'selected':len(merged),'new_values':newcount,'cache_hits':hits,'journal_events':seq,'independent_overlay_exact':True,'certificate_results_sha256':C.digest(proofs),'tree':tree,'numeric':numeric,'cumulative_passed':len(base['passed'])+len(final['new_passed']),'cumulative_splits':len(base['splits'])+len(final['new_splits']),'queue':len(queue),'cache':len(base['keys'])+len(final['new_cache'])})
        total_events+=seq;previous=final;write(out,result);print(json.dumps({'pid':pid,'structural089_batch':num,'new_values':newcount,'cache_hits':hits,'queue':len(queue)}),flush=True)
    C.require(len(previous['queue'])==snap['expected_queue'],'snapshot residual queue count differs')
    for path,h in {**inputs,**REG['protected_old_locks']}.items():C.require(C.file_sha(path)==h,'frozen/protected input changed during audit')
    absent=[]
    for p in sorted(allpids):
        try:os.kill(p,0)
        except ProcessLookupError:absent.append(p)
    C.require(len(absent)==len(allpids),'production PID still present or reused')
    result.update({'status':'all_frozen089_structure_values_qualified_partial_g105','through_batch':len(bundles),'new_values_verified':len(allnew),'independent_value_components':5*len(allnew),'certificate_count':len(result['certificates']),'corner_count':16*len(result['certificates']),'journal_events':total_events,'final_counts':{'passed':len(base['passed'])+len(previous['new_passed']),'splits':len(base['splits'])+len(previous['new_splits']),'queue':len(previous['queue']),'cache':len(base['keys'])+len(previous['new_cache']),'roots':52360},'production_pids_os_absent':absent,'global_g105_qualified':False,'old_locks_unchanged':True,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    write(out,result);print(json.dumps({k:result[k] for k in ['pid','status','through_batch','new_values_verified','certificate_count','corner_count','final_counts','wall_seconds']}),flush=True)


if __name__=='__main__':main()
