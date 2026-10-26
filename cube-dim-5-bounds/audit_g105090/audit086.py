#!/usr/bin/env python3
"""Fresh g105 check of every086 leaf; no H/gradient value is recomputed."""
import json,time,os,resource
from pathlib import Path
from fractions import Fraction as F
from itertools import product,zip_longest
from math import lcm
from hashlib import sha256
from interfaces import REG,HERE,common as C,arith as A


def write(obj):
    p=HERE/'repack086_result.json';q=Path(str(p)+'.tmp');q.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n');q.replace(p)


def state_record(row):
    return {'node':row['node'],'certificate':row['certificate'],'source070_leaf_index':row['source_leaf_index'],
            'source070_leaf_type':row['source_leaf_type'],'source070_certificate_sha256':row['source_certificate_sha256']}


def verify(row,info,record,roots):
    C.require(row['node']==list(info['node']) and row['source_certificate_sha256']==info['sha'] and row['source_leaf_type']==info['type'],'source leaf identity/type/certificate changed')
    cert=row['certificate'];C.require(cert['target']=='21/200000' and cert['epsilon']=='1/64' and cert['mechanism']=='full_frame_16_new_target_g105' and cert['quadratic']=='Q_full=I4/16' and cert['checked_points']==16,'wrong new-target/fullQ schema')
    C.require(cert['anchor']==info['anchor']==','.join(record['point']) and row['H_interval']==record['values'][0] and row['cached_bits']==record['bits'],'new-target certificate uses another source anchor/value')
    lo,hi=C.physical_box(roots,info['node']);C.require(list(map(F,row['box_lower']))==lo and list(map(F,row['box_upper']))==hi,'independent physical box differs')
    vls=[[F(a),F(b)] for a,b in record['values']];C.require(len(vls)==5 and all(a<=b for a,b in vls) and vls[0][0]>F(21,200000),'source numerical intervals not eligible for g105')
    gs=[[F(a),F(b)] for a,b in cert['loggradient_directed']];C.require(len(gs)==4 and all(a<=b for a,b in gs),'unordered gradient ratio')
    for (gl,gu),(dl,du) in zip(gs,vls[1:]):
        for h in vls[0]:
            for d in (dl,du):C.require(gl*h<=d<=gu*h,'new-target physical gradient/H interval fails')
    delta=F(cert['Delta_lower']);ratio=vls[0][0]/F(21,200000);L,U=A.log_bounds(ratio.numerator,ratio.denominator)
    C.require(delta>0 and delta<=F(L,A.S),'fresh g105 Delta is not below independent integer-log lower bound')
    point=list(map(F,record['point']));C.require(point[0]==1,'not physical x1=1')
    offsets=[(a-x,b-x) for a,x,b in zip(lo,point[1:],hi)]
    vd=lcm(*(x.denominator for pair in offsets for x in pair));gd=lcm(delta.denominator,*(x.denominator for pair in gs for x in pair));di=int(delta*gd)
    zz=[tuple(int(x*vd) for x in pair) for pair in offsets];gg=[tuple(int(x*gd) for x in pair) for pair in gs];rows=[]
    for bits in product((0,1),repeat=4):
        z=[zz[j][bits[j]] for j in range(4)];n=di*vd*vd+vd*sum(min(a*x,b*x) for (a,b),x in zip(gg,z))-256*gd*sum(x*x for x in z);rows.append((n,bits))
    minimum=min(x[0] for x in rows);margin=F(minimum,gd*vd*vd)
    C.require(F(cert['minimum_log_margin'])==margin and cert['passed']==(minimum>0) and (minimum,tuple(cert['worst_corner_bits'])) in rows,'new-target exact16corner mismatch')
    return margin


def main():
    began=time.monotonic();pid=os.getpid();print(json.dumps({'pid':pid,'stage':'086_frozen_source_and_new_target'}),flush=True)
    C.require(C.file_sha(REG['manifest086'])==REG['manifest086_sha256'],'086 manifest changed');m=json.loads(Path(REG['manifest086']).read_text());base=Path(REG['source086']).parent.parent
    inputs={REG['manifest086']:REG['manifest086_sha256'],REG['source086']:REG['source086_sha256'],REG['source070']:REG['source070_sha256']};inputs.update(m['external_files']);inputs.update(m['source_math_locks'])
    for b in m['blocks']:
        for name in ['registration','certificates','checkpoint','summary']:inputs[str(Path(b['directory'])/(name+('.jsonl' if name=='certificates' else '.json')))]=b[name+'_sha256']
    for p,h in {**inputs,**REG['protected_old_locks']}.items():C.require(C.file_sha(p)==h,'input/protected lock mismatch: '+p)
    wanted=set();mapper_hashes=[]
    for line in (base/'active_leaf_mapper.jsonl').open():
        r=json.loads(line);wanted.add(r['anchor']);mapper_hashes.append(C.digest(r))
    C.require(len(mapper_hashes)==97960,'active mapper count differs')
    old,new=C.Stream(REG['source070']),C.Stream(REG['source086'])
    for reader in (old,new):
        while True:
            key=reader.key();C.require(key is not None,'source cache absent')
            if key=='cache':break
            reader.value()
    cache={};cache_count=0;cache_hash=sha256();cache_hash.update(b'{');missing=object()
    for a,b in zip_longest(old.items(True),new.items(True),fillvalue=missing):
        C.require(a is not missing and b is not missing and a==b,'old180040 cache object changed')
        key,record=a
        if cache_count:cache_hash.update(b',')
        cache_hash.update(C.canonical(key));cache_hash.update(b':');cache_hash.update(C.canonical(record));cache_count+=1
        if key in wanted:cache[key]={'point':record['point'],'values':record['values'],'bits':record['bits']}
    cache_hash.update(b'}');C.require(cache_count==180040 and set(cache)==wanted,'source cache/selected anchor count differs')
    source_info=[];source_splits=[];source_roots=None
    while True:
        key=old.key()
        if key is None:break
        if key=='passed':
            for leaf in old.items():
                cert=leaf['certificate'];source_info.append({'node':C.node(leaf['node']),'sha':C.digest(cert),'type':cert.get('mechanism',cert.get('quadratic')),'anchor':cert['anchor']})
        elif key=='splits':source_splits=list(old.items())
        elif key=='roots':source_roots=old.value()
        else:old.value()
    old.close();C.require(len(source_info)==97960 and len(source_splits)==45600 and len(source_roots)==52360,'source070 finite counts changed')
    for i,line in enumerate((base/'active_leaf_mapper.jsonl').open()):
        row=json.loads(line);info=source_info[i];rec=cache[info['anchor']];lo,hi=C.physical_box(source_roots,info['node'])
        C.require(row['source_leaf_index']==i and row['node']==list(info['node']) and row['source_certificate_sha256']==info['sha'] and row['source_leaf_type']==info['type'] and row['anchor']==info['anchor'],'active source mapper mismatch')
        C.require(row['point']==rec['point'] and row['H_interval']==rec['values'][0] and row['cached_bits']==rec['bits'] and row['root']==source_roots[info['node'][0]] and list(map(F,row['box_lower']))==lo and list(map(F,row['box_upper']))==hi,'active mapper numeric/physical fields mismatch')
    risk=json.loads((base/'risk_order.json').read_text());order=risk['source_leaf_indices'];C.require(sorted(order)==list(range(97960)) and risk['source_sha256']==REG['source070_sha256'],'risk order is not complete source permutation')
    result={'task':'audit_g105090_repack086','pid':pid,'status':'finite_verified_prefix_only','target':'21/200000','epsilon':'1/64','inputs':inputs,'code_locks':{'audit086.py':C.file_sha(Path(__file__)),'registration.json':C.file_sha(HERE/'registration.json'),'interfaces.py':C.file_sha(HERE/'interfaces.py')},'old180040_H_recomputed':False,'cache_entries_directly_compared':cache_count,'cache_canonical_sha256':cache_hash.hexdigest(),'blocks':[],'new_target_all_leaf_count':97960,'new_target_geometry_and_integer_log':True}
    hashes=[None]*97960;statuses=[None]*97960;margins=[None]*97960;failed_rows={};seen=set();rank=0;previous=None
    status_lines=(base/'all_leaf_status.jsonl').open();failed_lines=(base/'uncertified_boxes.jsonl').open()
    for block in m['blocks']:
        directory=Path(block['directory']);br=json.loads((directory/'registration.json').read_text());cp=json.loads((directory/'checkpoint.json').read_text())
        C.require(br['parent_checkpoint_sha256']==cp['parent_checkpoint_sha256']==previous,'086 checkpoint chain differs')
        C.require(br['selected_source_leaf_indices']==order[rank:rank+block['tested']] and br['target']=='21/200000' and br['max_cases']==2048,'086 registered finite prefix differs')
        proofs=[];accepted=failed=0
        for line in (directory/'certificates.jsonl').open():
            row=json.loads(line);i=row['source_leaf_index'];C.require(i not in seen and i==order[rank],'duplicated/omitted/reordered source leaf');seen.add(i);rank+=1
            margin=verify(row,source_info[i],cache[source_info[i]['anchor']],source_roots);passed=margin>0;statuses[i]=passed;margins[i]=str(margin);hashes[i]=C.digest(state_record(row));proofs.append([i,row['node'],str(margin),passed])
            expected_status={'certificate_block':block['block'],'minimum_log_margin':str(margin),'node':row['node'],'passed':passed,'source_leaf_index':i,'source_leaf_type':row['source_leaf_type']}
            C.require(json.loads(status_lines.readline())==expected_status,'all-leaf status mapper differs')
            if passed:accepted+=1
            else:
                failed+=1;failed_rows[i]=row;C.require(json.loads(failed_lines.readline())==row,'failed-box ledger differs from exact test')
        C.require(accepted==block['passed']==cp['passed'] and failed==block['uncertified']==cp['uncertified'] and rank==cp['next_rank']==block['next_rank'],'block verdict/count checkpoint differs')
        C.require(cp['certificates_sha256']==block['certificates_sha256'] and cp['registration_sha256']==block['registration_sha256'] and cp['summary_sha256']==block['summary_sha256'],'checkpoint artifact links differ')
        previous=block['checkpoint_sha256'];result['blocks'].append({'block':block['block'],'tested':len(proofs),'passed':accepted,'uncertified':failed,'independent_exact_results_sha256':C.digest(proofs),'checkpoint_sha256':previous})
        write(result);print(json.dumps({'pid':pid,'086_block':block['block'],'verified_leaves':rank,'passed':sum(status is True for status in statuses)}),flush=True)
    C.require(rank==97960 and len(seen)==97960 and sum(statuses)==96104 and len(failed_rows)==1856 and not status_lines.readline() and not failed_lines.readline(),'complete g105 target verdict counts differ');status_lines.close();failed_lines.close()
    locations={};passed_nodes=[];queue_nodes=[];metadata={};pending=[];new_roots=None
    while True:
        key=new.key()
        if key is None:break
        if key=='passed':
            prior_index=-1
            for j,record in enumerate(new.items()):
                i=record['source070_leaf_index'];C.require(statuses[i] is True and i>prior_index and C.digest(record)==hashes[i],'086 accepted state mapping/certificate changed');prior_index=i;locations[i]={'kind':'passed','index':j};passed_nodes.append(C.node(record['node']))
        elif key=='pending_certificates086':
            pending=new.value()
            for j,record in enumerate(pending):
                i=record['source070_leaf_index'];C.require(statuses[i] is False and i not in locations and C.digest(record)==hashes[i],'086 pending source/certificate mapping changed');locations[i]={'kind':'queue','index':j}
        elif key=='queue':queue_nodes=list(map(C.node,new.value()))
        elif key=='roots':new_roots=new.value()
        elif key=='splits':C.require(list(new.items())==source_splits,'old-target split history records changed')
        else:metadata[key]=new.value()
    new.close();C.require(len(locations)==97960 and len(passed_nodes)==96104 and len(queue_nodes)==1856 and source_roots==new_roots,'partial state cardinality/root mapping differs')
    sorted_failed=sorted(failed_rows.values(),key=lambda r:(F(r['H_interval'][0]),F(r['certificate']['minimum_log_margin']),*r['node']))
    C.require([C.node(r['node']) for r in sorted_failed]==queue_nodes and [C.node(r['node']) for r in pending]==queue_nodes,'pending physical order differs')
    for i,line in enumerate((base/'state_leaf_mapper.jsonl').open()):
        row=json.loads(line);C.require(row['source_leaf_index']==i and row['new_target']=='21/200000' and row['new_state_location']==locations[i],'new state source-index location differs');row.pop('new_target');row.pop('new_state_location');C.require(C.digest(row)==mapper_hashes[i],'source leaf type/mapper fields were rewritten')
    pending_order=json.loads((base/'pending_order089.json').read_text());C.require(pending_order['target']=='21/200000' and pending_order['epsilon']=='1/64','wrong089pending target')
    C.require([C.node(r['node']) for r in pending_order['records']]==queue_nodes,'089 pending order differs')
    for record,row in zip(pending_order['records'],sorted_failed):
        i=row['source_leaf_index'];C.require(record['source_leaf_index']==i and record['source_leaf_type']==source_info[i]['type'] and record['anchor']==row['certificate']['anchor'] and record['H_interval']==row['H_interval'] and record['parent_cached_H_lower']==row['H_interval'][0] and record['parent_full_frame_minimum_log_margin']==row['certificate']['minimum_log_margin'],'089 pending cache/failed-certificate source fields differ')
        C.require(record['root']==source_roots[record['node'][0]] and record['box_lower']==row['box_lower'] and record['box_upper']==row['box_upper'] and record['point']==cache[source_info[i]['anchor']]['point'],'089 pending physical point mapping differs')
    tree=C.validate_tree(new_roots,passed_nodes+queue_nodes,{C.node(r['node']) for r in source_splits},'g105086partial')
    C.require(metadata['target']=='21/200000' and metadata['epsilon']=='1/64' and metadata['model']=='alpha1' and metadata['active'] is None and metadata['inflight']=={} and metadata['resume_jobs']=={} and metadata['counterexample'] is None,'partialstate model/live/counter fields differ')
    C.require(metadata['source070']['path']==REG['source070'] and metadata['source070']['sha256']==REG['source070_sha256'] and metadata['source070']['leaf_mapper_sha256']==m['external_files'][str(base/'state_leaf_mapper.jsonl')],'partial source link differs')
    expected_fp=C.digest({'requalification_registration_sha256':m['external_files'][str(base/'requalification_registration.json')],'last_checkpoint_sha256':previous,'state_leaf_mapper_sha256':m['external_files'][str(base/'state_leaf_mapper.jsonl')],'target':'21/200000','source070_sha256':REG['source070_sha256']})
    C.require(metadata['fingerprint086']==expected_fp and metadata['qualification086']['accepted']==96104 and metadata['qualification086']['pending']==1856 and metadata['qualification086']['g105_closed_cover'] is False,'new-target fingerprint/qualification differs')
    for p,h in {**inputs,**REG['protected_old_locks']}.items():C.require(C.file_sha(p)==h,'input/protected byte lock changed')
    result.update({'status':'086_new_target_all97960_qualified_partial','verified_leaves':rank,'verified_corners':16*rank,'passed':96104,'pending':1856,'tree':tree,'old_splits_directly_preserved':45600,'all_new_state_locations_verified':97960,'fingerprint086_independently_verified':True,'minimum_passed_exact_margin':str(min(F(v) for v,s in zip(margins,statuses) if s)),'g105_closed_cover':False,'global_g105_qualified':False,'input_and_old_locks_unchanged':True,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    write(result);print(json.dumps({k:result[k] for k in ['pid','status','verified_leaves','verified_corners','passed','pending','wall_seconds','peak_rss_bytes']}),flush=True)


if __name__=='__main__':main()
