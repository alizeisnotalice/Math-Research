#!/usr/bin/env python3
"""After value approval: new Delta/gradient/corner, journal, closed-tree audit."""
import sys
sys.dont_write_bytecode=True
import argparse,copy,os,time
from fractions import Fraction as F
from itertools import product
from math import lcm
import common as C
import pf128 as P
from bind import raw_gate

def certificate(n,c,rec,roots,task,ancestry):
    assert c['lower_bound_object']=='G' and c['component_function']=='H_1/128' and c['leaf_epsilon']==rec['epsilon']=='1/128'
    assert c['target']==rec['qualification_target']=='11/100000' and c['point_key']==rec['pair_cache_key'][1]==','.join(rec['point'])
    assert c['mechanism']=='129_anchor_Jensen_g105_U128' and c['checked_points']==16
    proof=task['H64_g105_ancestor'];assert c['H64_g105_ancestor_proof']==proof==ancestry['proof']
    assert ancestry['task_node']==task['node']==list(n) and ancestry['point']==rec['point'] and ancestry['ancestor_packet']['epsilon']=='1/64'
    anc=C.node(proof['node']);assert anc==C.node(ancestry['ancestor_packet']['node']) and anc[0]==n[0] and n[1].startswith(anc[1])
    assert c['H64_lower_source_node']==list(anc) and c['H64_lower']==proof['lower']==ancestry['ancestor_packet']['qualified_lower']=='21/200000'
    vals=[list(map(F,x)) for x in rec['values']];assert len(vals)==5 and vals[0][0]>0
    hu=vals[0][1];assert F(c['H128_upper'])==hu and hu>=F(21,200000)
    q=1024-768*F(21,200000)/hu;assert F(c['penalty_coefficient'])==q and 256<=q<1024
    gs=[list(map(F,x)) for x in c['loggradient_directed']];assert len(gs)==4 and all(a<=b for a,b in gs)
    for (gl,gu),(dl,du) in zip(gs,vals[1:]):
        for h in vals[0]:
            for d in (dl,du):assert gl*h<=d<=gu*h
    delta=F(c['Delta_lower']);ratio=vals[0][0]/C.TARGET;loglo,loghi=P.logarithm(ratio.numerator,ratio.denominator)
    assert delta<=F(loglo,P.S)
    lo,hi=C.box(roots,n);point=list(map(F,rec['point']));assert point[0]==1
    al,ah=C.box(roots,list(anc));assert all(a<=x<=b for a,x,b in zip(al,point[1:],ah))
    offsets=[(a-x,b-x) for a,x,b in zip(lo,point[1:],hi)]
    vd=lcm(*(v.denominator for pair in offsets for v in pair));gd=lcm(delta.denominator,q.denominator,*(v.denominator for pair in gs for v in pair))
    zz=[[int(v*vd) for v in pair] for pair in offsets];gg=[[int(v*gd) for v in pair] for pair in gs];di=int(delta*gd);qq=int(q*gd);rows=[]
    assert F(qq,gd)==q
    for bits in product((0,1),repeat=4):
        z=[zz[j][bits[j]] for j in range(4)]
        num=di*vd*vd+vd*sum(min(a*v,b*v) for (a,b),v in zip(gg,z))-qq*sum(v*v for v in z)
        rows.append((num,bits))
    minimum,bits=min(rows);margin=F(minimum,gd*vd*vd)
    assert c['passed']==(minimum>0) and F(c['minimum_log_margin'])==margin and c['worst_corner_bits']==list(bits)
    return dict(node=list(n),certificate_sha256=C.digest(c),record_sha256=C.digest(rec),margin=str(margin),passed=minimum>0,corners=16,
        q=str(q),H64_ancestor_packet=ancestry['ancestor_packet'],anchor_inside_box=all(a<=x<=b for a,x,b in zip(lo,point[1:],hi)),Delta_negative=delta<0,HL_above_target=vals[0][0]>C.TARGET)


def numeric_gate(batch,entries):
    out=C.SCRATCH/('numeric_batch%04d'%batch);index=C.read(out/'index.json')
    assert index['status']=='all_five_components_strictly_verified' and index['verified_points']==index['global_prefix']==len(entries)
    C.check_locks(index['input_locks']);seen=set()
    for plan,block in zip(index['plans'],index['blocks']):
        path=out/('block%04d.json'%plan['block']);assert C.sha(path)==block['sha256'];s=C.read(path)
        assert s['status']=='complete' and s['input_locks']==index['input_locks'] and len(s['points'])==plan['end']-plan['start']
        for proof in s['points']:
            i=proof['index'];assert i not in seen;seen.add(i);e=entries[i]
            assert proof['point']==e['point'] and proof['raw_sha256']==e['point_checkpoint_sha256'] and proof['pair_cache_key']==e['record']['pair_cache_key']
            assert len(proof['independent_values'])==5
            for ball,original in zip(proof['independent_values'],e['record']['values']):
                lo,hi=map(F,original);a,b=map(F,ball);assert lo<a<=b<hi
    assert seen==set(range(len(entries))) and all(x['exitcode']==0 and x['alive'] is False for x in index['pool_exits'])
    return C.sha(out/'index.json')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--registration',required=True);ap.add_argument('--reviewed-registration-sha',required=True);a=ap.parse_args();began=time.monotonic()
    rp=__import__('pathlib').Path(a.registration);assert C.sha(rp)==a.reviewed_registration_sha;r=C.read(rp);batch=r['batch']
    resultpath=C.HERE/('structure_batch%04d.json'%batch);assert not resultpath.exists(),'do not rerun a completed structural/Delta gate'
    lock=C.read(C.HERE/'source_lock.json');assert C.sha(C.HERE/'source_lock.json')==r['source_lock_sha256'];C.check_locks(lock['code_locks']);C.check_locks(lock['input_locks']);C.check_locks(r['input_locks'])
    tasks=C.read(r['tasks']);pm=C.read(r['point_manifest']);entries=pm['records'];numeric_sha=numeric_gate(batch,entries)
    geo=C.read(C.SCRATCH/'source_geometry.json');roots=geo['roots'];b=C.PROD/('future_batch%04d'%batch);initial=C.read(b/'initial_overlay.json');final=C.read(b/'state.json');reg=C.read(b/'registration.json')
    if batch==1:
        assert initial['queue']==[t['node'] for t in C.read(C.PROD/'tasks.json')]
        assert initial['new_cache_by_epsilon']=={'1/128':{}} and initial['new_passed']==initial['new_splits']==[]
    else:
        prior=C.read(C.HERE/'structure_batch0001.json');assert prior['status']=='new_values_corners_journal_and_partial_tree_qualified'
        old=C.read(C.PROD/'future_batch0001/state.json')
        for key in old:
            if key not in ('next_merge','journal_seq','journal_hash','batch','status'):assert initial[key]==old[key]
        assert initial['queue'][:len(tasks)]==[t['node'] for t in tasks]
    assert initial['target']=='11/100000' and initial['seed_source_manifest_sha256']==lock['input_locks'][str(C.AUTHOR/'result_manifest.json')]
    expected=copy.deepcopy(initial);seq=0;previous='0'*64;dispatch={};received={};certs=[]
    for event in C.lines(b/'journal.jsonl'):
        body=event['body'];kind=event['kind'];assert event['seq']==seq+1 and event['previous']==previous
        assert C.digest({k:v for k,v in event.items() if k!='hash'})==event['hash'];seq+=1;previous=event['hash']
        if kind=='dispatch':
            task=body['task'];i=task['batch_ordinal'];assert i==len(dispatch) and task==tasks[i] and expected['queue'].pop(0)==task['node']
            dispatch[i]=task;assert str(i) not in expected['inflight'];expected['inflight'][str(i)]={'task':task};expected['counts']['new_anchor_starts']+=1
        elif kind=='received':
            i=body['ordinal'];got=body['returned'];task=dispatch[i];e=entries[i];raw=raw_gate(e,task)
            assert i not in received and got['kind']=='complete' and got['epsilon']=='1/128' and got['target']=='11/100000' and got['pair_cache_key']==task['pair_cache_key']
            assert got['record']==e['record']==raw['record'] and got['attempts']==raw['attempts'] and got['calls']==raw['calls'] and got['pid']==raw['pid']
            received[i]=got;expected['inflight'][str(i)]['returned']=got;expected['counts']['kernel_segment_calls']+=got['calls']
        elif kind=='merge':
            i=body['ordinal'];assert i==expected['next_merge'];expected['next_merge']+=1;ent=expected['inflight'].pop(str(i));rec=ent['returned']['record'];task=ent['task'];n=C.node(task['node']);c=body['certificate']
            key=task['pair_cache_key'][1];assert key not in expected['new_cache_by_epsilon']['1/128'];expected['new_cache_by_epsilon']['1/128'][key]=rec;expected['counts']['new_complete']+=1
            certs.append(certificate(n,c,rec,roots,task,geo['ancestry'][str(task['ordinal'])]))
            if c['passed']:expected['new_passed'].append(dict(node=list(n),leaf_epsilon='1/128',lower_bound_object='G',target='11/100000',certificate=c))
            else:
                kids=[[n[0],n[1]+'0'],[n[0],n[1]+'1']];expected['new_splits'].append(dict(node=list(n),axis=len(n[1])%4,children=kids,leaf_epsilon='1/128',target='11/100000',failed_certificate=c));expected['queue'].extend(kids)
        elif kind=='finish':
            assert body['status']=='finite_original_prefix_complete' and len(dispatch)==len(received)==len(certs)==len(tasks)
            pool=body['pool_evidence'];assert len(pool)==15 and all(x['exitcode']==0 and x['alive'] is False for x in pool)
            assert {x['pid'] for x in received.values()}<=set(x['pid'] for x in pool);expected.update(body)
        else:raise AssertionError(('unexpected/fault journal event',kind))
        expected['journal_seq']=seq;expected['journal_hash']=previous
    assert seq==3*len(tasks)+1 and expected==final and final['faults']==final['obstructions_H128']==[] and not final['inflight']
    assert [t['node'] for t in tasks]==[dispatch[i]['node'] for i in range(len(tasks))]
    oldpass=geo['passed'];oldsplits=geo['splits'];tree=C.tree(roots,oldpass+[x['node'] for x in final['new_passed']]+final['queue'],oldsplits+[x['node'] for x in final['new_splits']])
    # Original task set is the frozen seed frontier; no geometry child occurs
    # in the dispatched set, even though its closed children stay in queue.
    all_tasks=C.read(C.PROD/'tasks.json');assert set(map(lambda x:C.node(x['node']),all_tasks))==set(map(C.node,geo['queue']))
    C.check_locks(lock['code_locks']);C.check_locks(lock['input_locks']);C.check_locks(r['input_locks'])
    result=dict(status='new_values_corners_journal_and_partial_tree_qualified',batch=batch,pid=os.getpid(),seconds=time.monotonic()-began,
        numerical_index_sha256=numeric_sha,registration_sha256=a.reviewed_registration_sha,new_points=len(tasks),new_components=5*len(tasks),new_Delta_checks=len(tasks),integer_corners=16*len(tasks),
        passed_this_batch=sum(x['passed'] for x in certs),geometric_splits_this_batch=sum(not x['passed'] for x in certs),certificates=certs,journal_events=seq,exact_overlay_reconstructed=True,tree=tree,
        old_values_recomputed=False,old_Delta_recomputed=False,original_seed_frontier_only=True,no_geometry_child_dispatched=True,global_g110_qualified=False,S_credit=0,formal_C5_update=False,main_final_acceptance_pending=True)
    C.write(resultpath,result);print(__import__('json').dumps({k:v for k,v in result.items() if k!='certificates'}),flush=True)

if __name__=='__main__':main()
