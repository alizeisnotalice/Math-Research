#!/usr/bin/env python3
"""PREPARED148: inherited certificates + one new exact whole-root integration.

No H/kernel/PF/log/gradient/corner re-evaluation. JSON selector is the frozen
non-numerical137 mmap reader, copied byte-for-byte as json_sections.py.
"""
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from fractions import Fraction as Q
from pathlib import Path
import argparse,gzip,hashlib,json,os,time
from json_sections import Document

HERE=Path(__file__).resolve().parent
TARGET='11/100000'
MODEL='n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
STAGES=('base','137','138','141','143','145','147')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        while b:=f.read(8*1024*1024):h.update(b)
    return h.hexdigest()
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def read(path):return json.loads(Path(path).read_bytes())
def emit(path,value):
    with Path(path).open('x') as f:json.dump(value,f,sort_keys=True,indent=2);f.write('\n')
def key(node):
    r,p=node
    assert type(r) is int and 0<=r<52360 and isinstance(p,str) and not set(p)-{'0','1'}
    return r,p
def sections(path,names):
    with Document(path) as doc:
        spans=dict(doc.object_spans())
        return {n:doc.load_span(spans[n]) for n in names}
def lines(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        for row in f:yield json.loads(row)
def lock_pairs(value):
    """Both historic absolute-path maps and named {path,sha256} rows."""
    out=set()
    if isinstance(value,dict):
        if isinstance(value.get('path'),str) and isinstance(value.get('sha256'),str):out.add((value['path'],value['sha256']))
        for k,v in value.items():
            if isinstance(k,str) and k.startswith('/') and isinstance(v,str) and len(v)==64:out.add((k,v))
            out.update(lock_pairs(v))
    elif isinstance(value,list):
        for v in value:out.update(lock_pairs(v))
    return out
def select(value,pointer):
    for p in pointer:value=value[p]
    return value

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reviewed-plan-sha',required=True)
    ap.add_argument('--index',type=Path,required=True)
    ap.add_argument('--reviewed-index-sha',required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    assert sha(HERE/'plan.json')==a.reviewed_plan_sha and sha(a.index)==a.reviewed_index_sha
    plan=read(HERE/'plan.json');index=read(a.index)
    assert a.reviewed_index_sha==plan['index_sha256']
    assert sha(__file__)==plan['code_sha256'] and sha(HERE/'json_sections.py')==plan['reader_sha256']
    assert plan['maximum_global_integrations']==1 and index['all_components_final'] is True
    assert set(index['components'])==set(STAGES)
    paths={name:Path(row['path']) for name,row in index['sources'].items()}
    # Independent I/O locks may hash in parallel; mathematical integration is exact and serial.
    def filecheck(item):
        name,row=item;got=sha(row['path']);assert got==row['sha256'],name;return name,got
    with ThreadPoolExecutor(max_workers=15) as pool:before=dict(pool.map(filecheck,index['sources'].items()))
    a.output.mkdir(exist_ok=False)
    emit(a.output/'registration.json',{'task':'audit_global148_one_final_forest',
        'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
        'plan_sha256':a.reviewed_plan_sha,'index_sha256':a.reviewed_index_sha,
        'code_sha256':sha(__file__),'reader_sha256':sha(HERE/'json_sections.py'),
        'source_locks':index['sources'],'hash_workers':15,'math_workers':1,
        'new_H_PF_log_gradient_or_corner_calls':0,'old_baseline_recomputations':0,
        'maximum_global_integrations':1,'wall_deadline':None})
    begin=time.monotonic();component_receipts={};lineage=set()
    required_sources={
        'base':{'seed','mapper','state123','curvature128'},
        '137':{'cross137'},'138':{'numeric138a','numeric138'},
        '141':{'cross141'},'143':{'numeric143'},
        '145':{'cross145'},'147':{'numeric147'}}
    for stage in STAGES:
        cfg=index['components'][stage]
        assert set(cfg['production_sources'])==required_sources[stage]
        inherited=set()
        for source in cfg['lock_documents']:inherited.update(lock_pairs(read(paths[source])))
        for source in cfg['production_sources']:
            row=index['sources'][source]
            assert (row['path'],row['sha256']) in inherited,(stage,source,'not bound by inherited independent audit')
            lineage.add(source)
        for requirement in cfg['qualification_requirements']:
            assert select(read(paths[requirement['source']]),requirement['pointer'])==requirement['equals']
        assert cfg['qualification_requirements'],stage
        component_receipts[stage]={'production_sources':cfg['production_sources'],
            'qualification_requirements':cfg['qualification_requirements'],
            'lock_documents':cfg['lock_documents']}
    assert lineage==set().union(*required_sources.values())
    roots=[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
    seed=sections(paths['seed'],['roots','target','model','passed','splits','queue','inflight','resume_jobs','source_passed_mapper_sha256'])
    assert seed['roots']==roots and len(roots)==52360 and seed['target']==TARGET and seed['model']==MODEL
    assert not seed['inflight'] and not seed['resume_jobs']
    mapper=read(paths['mapper']);assert len(mapper)==len(seed['passed'])==98007
    assert seed['source_passed_mapper_sha256']==index['sources']['mapper']['sha256']
    accepted={};internal={};epsilon_count=Counter();mechanisms=Counter();stage_count=Counter();stage_snapshots=[]
    minmargin=None;min_witnesses=[]
    proofraw=(a.output/'leaf_lineage.jsonl.gz').open('xb')
    proofstream=gzip.GzipFile(filename='',mode='wb',fileobj=proofraw,mtime=0)
    def cert_type(c,epsilon):
        assert epsilon in ('1/64','1/128') and c['target']==TARGET and c['passed'] is True
        margin=Q(c['minimum_log_margin']);assert c['checked_points']==16 and margin>0
        q=Q(c['penalty_coefficient']);assert 0<q<=1024
        rule=c.get('curvature_rule');mechanism=c.get('mechanism')
        if rule is not None:
            if epsilon=='1/64':assert rule=='075_H64_Q256' and q==256
            elif rule=='129_anchor_Jensen_g105_U128':assert 256<=q<1024
            else:assert rule=='075_H128_Q1024_no_H64g105_ancestor' and q==1024
            assert c['scope']=='single source allclosedcorners;not cornerwise switching;not integral mean code'
            return rule,margin
        assert c['component_function'] in ('H_epsilon','H_1/128')
        if c['component_function']=='H_1/128':assert epsilon=='1/128'
        if mechanism in ('075_full_frame_16','075_saved_parent_anchor_16','075_saved_global_fulljet_16'):
            assert c['component_function']=='H_epsilon' and q==(256 if epsilon=='1/64' else 1024)
        else:
            assert mechanism in ('fixed-t075_then_anchor-Jensen_time-mass','129_anchor_Jensen_g105_U128')
            assert epsilon=='1/128' and 256<=q<1024
        return mechanism,margin
    def accept(node,c,epsilon,stage,source):
        nonlocal minmargin,min_witnesses
        n=key(node);assert n not in accepted and n not in internal,('duplicate or internal leaf',n)
        rule,margin=cert_type(c,epsilon)
        declared=c.get('leaf_epsilon',c.get('epsilon',epsilon));assert declared==epsilon
        accepted[n]=(epsilon,stage);epsilon_count[epsilon]+=1;mechanisms[rule]+=1;stage_count[stage]+=1
        row={'node':list(n),'epsilon':epsilon,'proved_component':'H_'+epsilon,
            'target':TARGET,'proof_rule':rule,'stage':stage,'production_source':source,
            'source_sha256':index['sources'][source]['sha256'],'certificate_sha256':digest(c),
            'minimum_log_margin':str(margin)}
        if minmargin is None or margin<minmargin:minmargin=margin;min_witnesses=[row]
        elif margin==minmargin:min_witnesses.append(row)
        proofstream.write((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode())
    def add_splits(rows):
        out=set()
        for row in rows:
            n=key(row['node']);r,p=n
            assert n not in internal and n not in accepted and n not in out
            assert row['axis']==len(p)%4 and row['children']==[[r,p+'0'],[r,p+'1']]
            internal[n]=True;out.add(n)
        return out
    def unique_nodes(rows):
        result=[key(n) for n in rows];assert len(set(result))==len(result);return set(result)
    pending=unique_nodes(seed['queue']);add_splits(seed['splits'])
    for i,(row,m) in enumerate(zip(seed['passed'],mapper)):
        assert m['seed_passed_index']==i and row['node']==m['node'] and digest(row)==m['leaf_sha256']
        e=row['leaf_epsilon'];assert e==m['epsilon']==m['pair_cache_key'][0]
        c=row['certificate'];p=c.get('point_key',c.get('pair_cache_key',[None,None])[1]);assert m['pair_cache_key']==[e,p]
        accept(row['node'],c,e,'base','seed')
    del seed['passed']
    def overlay(source,stage,expected_passes,expected_splits):
        nonlocal pending
        obj=sections(paths[source],['target','new_passed','new_splits','queue','inflight','faults','obstructions_H128','pending_evaluation_epsilon'])
        assert obj['target']==TARGET and obj['pending_evaluation_epsilon']=='1/128'
        assert not obj['inflight'] and not obj['faults'] and not obj['obstructions_H128']
        p=unique_nodes([x['node'] for x in obj['new_passed']]);s=unique_nodes([x['node'] for x in obj['new_splits']])
        assert p.isdisjoint(s) and p|s==pending and len(p)==expected_passes and len(s)==expected_splits
        children=unique_nodes([child for x in obj['new_splits'] for child in x['children']])
        assert children==unique_nodes(obj['queue'])
        for row in obj['new_passed']:accept(row['node'],row['certificate'],row['leaf_epsilon'],stage,source)
        add_splits(obj['new_splits']);pending=children
        stage_snapshots.append({'stage':stage,'source':source,'passed':len(accepted),'splits':len(internal),'pending':len(pending)})
    overlay('state123','base',572,5665)
    rows=list(lines(paths['curvature128']));assert len(rows)==11330
    assert unique_nodes([x['child_node'] for x in rows])==pending
    pending=set()
    for row in rows:
        if row['passed']:accept(row['child_node'],row,row['epsilon'],'base','curvature128')
        else:pending.add(key(row['child_node']))
    assert stage_count['base']==104482 and len(pending)==5427 and len(internal)==57549
    stage_snapshots.append({'stage':'base014','passed':len(accepted),'splits':len(internal),'pending':len(pending)})
    def cross(source,stage,count):
        nonlocal pending
        obj=read(paths[source]);assert not obj['new_splits'] and not obj['faults']
        new=list(obj['qualified_children'].values());P=unique_nodes([c['child_node'] for c in new]);R=unique_nodes(obj['pending_children'])
        assert len(P)==count and P.isdisjoint(R) and P|R==pending
        for c in new:accept(c['child_node'],c,c['epsilon'],stage,source)
        pending=R;stage_snapshots.append({'stage':stage,'passed':len(accepted),'splits':len(internal),'pending':len(pending)})
    cross('cross137','137',1153)
    overlay('numeric138','138',3099,1175)
    cross('cross141','141',671)
    overlay('numeric143','143',1562,117)
    cross('cross145','145',108)
    overlay('numeric147','147',126,0)
    assert not pending and len(accepted)==111201 and len(internal)==58841
    assert set(accepted).isdisjoint(internal)
    # Reconstruct every strict leaf prefix, independently of producer's heap/Kraft code.
    derived=set()
    by=defaultdict(list)
    for r,p in accepted:
        by[r].append(p)
        for j in range(len(p)):derived.add((r,p[:j]))
    assert derived==set(internal) and set(by)==set(range(52360))
    rootproofs=[];maxdepth=0
    for r in range(52360):
        ps=by[r];D=max(map(len,ps));cursor=0
        intervals=sorted((int(p,2) if p else 0,len(p),p) for p in ps)
        intervals=sorted((v<<(D-d),(v+1)<<(D-d),p) for v,d,p in intervals)
        for left,right,p in intervals:assert left==cursor and right>left;cursor=right
        assert cursor==1<<D
        maxdepth=max(maxdepth,D)
        rootproofs.append({'root_id':r,'root_index':roots[r],'leaves':len(ps),'depth':D,'exact_tiling_end':cursor})
    proofstream.close();proofraw.close()
    with (a.output/'root_coverage.jsonl.gz').open('xb') as raw, gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as f:
        for row in rootproofs:f.write((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode())
    # All admitted rules bound log(H_epsilon/target) by the saved concave function.
    # Its saved all-corner minimum is a whole-box lower bound; exp(m)>=1+m.
    assert minmargin is not None and minmargin>0 and min_witnesses
    gstar=Q(TARGET)*(1+minmargin)
    with ThreadPoolExecutor(max_workers=15) as pool:after=dict(pool.map(filecheck,index['sources'].items()))
    assert before==after
    emit(a.output/'results.json',{'task':'audit_global148','status':'complete_component_lineage_and_full_cover_qualified_pending_primary',
        'pid':os.getpid(),'seconds':time.monotonic()-begin,'roots':52360,'passed':len(accepted),'splits':len(internal),'pending':0,
        'max_depth':maxdepth,'epsilon_leaf_counts':dict(epsilon_count),'proof_rule_counts':dict(mechanisms),
        'stage_leaf_counts':dict(stage_count),'stage_frontiers':stage_snapshots,'components':component_receipts,
        'every_leaf_proves_same_epsilon_H_on_closed_box':True,
        'H128_nominal_lower_on_complete_F':TARGET,'H128_lower_on_Q_claim':False,
        'minimum_saved_log_margin':str(minmargin),'minimum_margin_witnesses':min_witnesses,
        'every_admitted_rule_qualifies_saved_log_margin_as_whole_box_lower':True,
        'optional_rational_face_H128_floor_for_primary':str(gstar),
        'optional_rational_G_cube_floor_for_primary':str(gstar),
        'nominal_G_cube_lower_qualified_for_primary':TARGET,'same029_general_L1_all_positive_common_radii':True,
        'all_S1_S2_and_closed_boundaries_included':True,'root_index_set_exact_and_all_roots_tiled':True,
        'C5_updated_by_this_agent':False,'new_H_PF_log_gradient_or_corner_calls':0,'old_baseline_recomputations':0,
        'all_source_hashes_unchanged':True,'leaf_lineage_sha256':sha(a.output/'leaf_lineage.jsonl.gz'),
        'root_coverage_sha256':sha(a.output/'root_coverage.jsonl.gz'),'registration_sha256':sha(a.output/'registration.json')})
    print('PASS: one global148 integration; exact full forest and every-leaf H128-on-F interface')

if __name__=='__main__':main()
