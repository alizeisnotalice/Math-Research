"""Fixed shards, durable per-case ACKs; successful prefix uses no repeated math."""
from pathlib import Path
from fractions import Fraction as F
import copy,os,resource,time,traceback
from common import read,sha,digest,chain,append,atomic
from math174 import source_q,corners
def initial(receiver):
    b=copy.deepcopy(receiver['base_best158']);assert b['first_pass115'] is None and F(b['best_margin110'])>0 and F(b['best_margin115'])<=0
    b['new173_tested_pairs']=0;b['new173_status']='base158_preserved';return b
def update(b,c):
    assert b['first_pass115'] is None;b['new173_tested_pairs']+=1
    if F(c['margin110'])>F(b['best_margin110']):b.update({'best_margin110':c['margin110'],'best_margin115':c['margin115'],'best_source_kind':'new173_band172_single_anchor','best_canonical_id':c['canonical_id'],'best_certificate_ref':c})
    if c['passed115']:b['first_pass115']=c;b['new173_status']='first115PASS_stop_receiver'
    else:b['new173_status']='best_base158_and_all_tested_newwholebox_including_failed'
def mass_case(cid,s,jet,c,times,rowhash):
    old=s['old155'];ref=s['original_g151_jet_ref'];p=old['point'];g=jet['values'][0]
    # Identity truth comes from accepted455 gate; these strict interface predicates remain mandatory.
    assert s['canonical_id']==cid and old['epsilon']==ref['epsilon']==jet['epsilon']=='1/64' and list(map(F,p))==list(map(F,ref['point']))==list(map(F,jet['point'])) and F(p[0])==1
    assert jet['minimal_projection_sha256']==digest({k:v for k,v in jet.items() if k!='minimal_projection_sha256'}) and jet['raw_record_sha256']==old['raw_record_sha256']==ref['raw_record_sha256']
    assert len(old['loggradient_directed'])==4 and all(F(a)<=F(b) for a,b in old['loggradient_directed'])
    packet=read(s['source171_metadata']['new_mass_donor']['packet_path']);end=packet['endpoints'];assert [e['t'] for e in end]==times
    q=source_q(g,[[e['H']['lower'],e['H']['upper']] for e in end],times,old['oldqualified155q'])
    assert c['canonical_id']==cid and c['epsilon']=='1/64' and c['point']==p and c['oldqualified155q']==old['oldqualified155q'] and list(map(F,c['original_g64_interval']))==list(map(F,g))
    for key,x in [('K_upper',q['K']),('L64_denominator',q['lower64']),('Q_band',q['band']),('Q_153_realH32_split',q['split']),('q_new',q['q']),('Q_153_uniform',F(512,3)),('original_g64_U_unchanged',F(g[1])),('real_H32_lower',F(end[32]['H']['lower']))]:assert F(c[key])==x,(cid,key)
    assert list(map(F,c['donor_H64_interval']))==[F(end[24]['H']['lower']),F(end[24]['H']['upper'])] and c['q_strictly_decreased']==(q['q']<F(old['oldqualified155q']))
    assert c['same_cutoff_exact_nonempty_overlap'] and c['denominator_uses_lower_not_upper'] and c['old_Delta_gamma_unchanged'] and c['infinite_tail_included'] and c['band_count']==49 and len(c['bands'])==49
    for j,(r,(m,w,v)) in enumerate(zip(c['bands'],q['terms']),24):
        assert r['a']==times[j] and r['b']==(times[j+1] if j<72 else 'infinity')
        assert F(r['H_upper_a'])==F(end[j]['H']['upper']) and F(r['H_lower_b'])==(F(end[j+1]['H']['lower']) if j<72 else 0)
        assert (F(r['mass_upper']),F(r['weight']),F(r['weighted_upper']))==(m,w,v)
    assert c['source155_body_sha256']==old['source155_body_sha256'] and c['original_g_raw_record_sha256']==old['raw_record_sha256'] and c['original_jet_projection_sha256']==jet['minimal_projection_sha256']
    return {'canonical_id':cid,'q':str(q['q']),'oldq':old['oldqualified155q'],'strict':q['q']<F(old['oldqualified155q']),'K':str(q['K']),'producer_row_hash':rowhash,'all49_products_and_tail_match':True}
def corner_case(pair,row,receiver,b,s,q,C,calculate):
    assert all(row[k]==pair[k] for k in ('pair_index','receiver_index','canonical_id'))
    if pair['prior_test_same_old155_q'] is not None:assert pair['prior_test_same_old155_q']['new_q_strict_decrease_required']
    if b['first_pass115'] is not None:kind='skip_after_firstPASS'
    elif pair['prior_test_same_old155_q'] is not None and not q['strict']:kind='skip_oldpair_q_unchanged'
    else:kind='evaluated'
    assert row['kind']==kind
    if kind=='evaluated':
        c=row['certificate'];old=s['old155'];assert c['node']==pair['node']==receiver['node'] and c['canonical_id']==pair['canonical_id'] and c['receiver_index']==pair['receiver_index'] and c['leaf_index']==pair['leaf_index']
        assert c['point']==old['point'] and c['epsilon']=='1/64' and F(c['q_new'])==F(q['q']) and c['source_mass_row_hash']==q['producer_row_hash'] and c['single_anchor_all16'] and c['checked_points']==16
        assert F(old['Delta115_lower'])+C==F(old['Delta110_lower'])
        if calculate:
            m,bits=corners(old['point'],old['Delta115_lower'],old['loggradient_directed'],q['q'],receiver['box_lower'],receiver['box_upper'])
            assert F(c['margin115'])==m and c['worst_corner_bits']==list(bits)
        assert F(c['margin110'])==F(c['margin115'])+C and c['passed115']==(F(c['margin115'])>0);update(b,c)
    else:assert not calculate or 'certificate' not in row
    return kind
def execute(phase,spec,plan,O,ready,result,start,release,pilot_sid,stop,resume):
    sid=spec['audit_shard'];dest=Path(O)/phase/f'shard{sid:02d}';dest.mkdir(parents=True,exist_ok=True);began=time.monotonic();state=None
    try:
        os.nice(5);fp={'plan_sha256':plan['self_sha256'],'phase':phase,'audit_shard':sid,'spec_sha256':digest(spec)};reg=dest/'registration.json'
        if reg.exists():assert resume and read(reg)['fingerprint']==fp
        else:atomic(reg,{'pid':os.getpid(),'fingerprint':fp})
        source=read(plan['registry_path']);C=F(read(plan['common_log_path'])['upper']);data=read(spec['input_path']);events=chain(spec['event_path']);best={};answers={};counts={'mass':0,'evaluated':0,'skip_after_firstPASS':0,'skip_oldpair_q_unchanged':0}
        if phase=='mass':
            jobs=data['source_ids'];jets=read(spec['g_input_path']);certs=read(spec['certificate_path']);assert len(events)==len(jobs) and set(certs)==set(jobs)
            for cid,r in zip(jobs,events):assert r['kind']=='mass_certificate' and r['canonical_id']==cid and certs[cid]=={'certificate':r['certificate'],'row_hash':r['hash']}
        else:
            jobs=data['pairs'];assert len(events)==len(jobs);best={k:initial(r) for k,r in data['receivers'].items()};answers=read(Path(O)/'independent_q_registry.json')
        journal=dest/'acks.jsonl';prior=chain(journal) if journal.exists() else [];assert not prior or resume;cp=read(dest/'checkpoint.json') if (dest/'checkpoint.json').exists() else None
        atcp={'last_hash':'0'*64,'counts':dict(counts),'best_sha256':digest(best)} if cp and cp['next_index']==0 else None
        for n,a in enumerate(prior):
            assert a['case_index']==n and a['producer_row_hash']==events[n]['hash']
            if phase=='mass':assert a['answer']['canonical_id']==jobs[n];answers[jobs[n]]=a['answer'];counts['mass']+=1
            else:
                pair=jobs[n];kind=corner_case(pair,events[n],data['receivers'][str(pair['receiver_index'])],best[str(pair['receiver_index'])],source[pair['canonical_id']],answers[pair['canonical_id']],C,False);assert a['kind']==kind;counts[kind]+=1
            if cp and n+1==cp['next_index']:atcp={'last_hash':a['hash'],'counts':dict(counts),'best_sha256':digest(best)}
        n=len(prior);previous=prior[-1]['hash'] if prior else '0'*64;assert n<=len(jobs)
        if cp:assert cp['fingerprint']==fp and cp['next_index']<=n and atcp=={k:cp[k] for k in ('last_hash','counts','best_sha256')}
        state={'fingerprint':fp,'next_index':n,'last_hash':previous,'counts':dict(counts),'best_sha256':digest(best)}
        ready.put({'sid':sid,'pid':os.getpid(),'next_index':n,'total':len(jobs)});start.wait();pilot=pilot_sid.value==sid and n<len(jobs);released=False;before=n
        with journal.open('ab') as f:
            while n<len(jobs) and not stop.is_set():
                if not released and not pilot:release.wait();released=True
                if stop.is_set():break
                r=events[n];ack={'seq':n+1,'previous':previous,'case_index':n,'producer_row_hash':r['hash']}
                if phase=='mass':
                    cid=jobs[n];a=mass_case(cid,source[cid],jets[cid],r['certificate'],plan['times'],r['hash']);answers[cid]=a;ack.update({'kind':'mass','answer':a});counts['mass']+=1
                else:
                    pair=jobs[n];kind=corner_case(pair,r,data['receivers'][str(pair['receiver_index'])],best[str(pair['receiver_index'])],source[pair['canonical_id']],answers[pair['canonical_id']],C,True);ack['kind']=kind;counts[kind]+=1
                previous=append(f,ack);n+=1;state.update({'next_index':n,'last_hash':previous,'counts':dict(counts),'best_sha256':digest(best)})
                if phase=='mass' or n%128==0 or pilot:atomic(dest/'checkpoint.json',state)
                if pilot:result.put({'type':'pilot_pass','sid':sid,'next_index':n,'phase':phase});pilot=False;release.wait();released=True
        state.update({'status':'complete' if n==len(jobs) else 'stopped_successful_prefix'});atomic(dest/'checkpoint.json',state)
        # Producer checkpoint/best are checked against the independently reconstructed entire stream.
        if n==len(jobs):
            pc=read(spec['checkpoint_path']);ps=read(spec['summary_path']);pr=read(spec['registration_path']);assert pc['status']==ps['status']=='complete' and pc['next_index']==len(jobs) and pc['last_hash']==(events[-1]['hash'] if events else '0'*64)
            assert pc['fingerprint']==ps['fingerprint']==pr['fingerprint'] and pc['counts']==ps['counts']==counts and pc['best_sha256']==ps['best_sha256']==digest(best)
            if phase=='corners':assert read(spec['best_path'])==best
        atomic(dest/('q_answers.json' if phase=='mass' else 'receiver_best.json'),answers if phase=='mass' else best)
        u=resource.getrusage(resource.RUSAGE_SELF);out={**state,'type':'summary','sid':sid,'pid':os.getpid(),'new_cases':n-before,'cpu_user':u.ru_utime,'cpu_system':u.ru_stime,'peak_rss_bytes':u.ru_maxrss,'seconds':time.monotonic()-began};atomic(dest/'summary.json',out);result.put(out)
    except BaseException as e:
        stop.set();start.set();release.set();out={'type':'fault','sid':sid,'pid':os.getpid(),'exception':repr(e),'traceback':traceback.format_exc(),'prefix':state};atomic(dest/f'fault{os.getpid()}.json',out);result.put(out)
