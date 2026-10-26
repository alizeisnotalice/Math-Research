#!/usr/bin/env python3
"""Prepare fixed metadata/index and source locks only. Never reads bulk H."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
import ast,gzip,hashlib,json,os
H=Path(__file__).resolve().parent
B=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_integral156')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def main():
    assert not (H/'export_plan.json').exists()
    inv=read(H/'inventory.json');cut=read(H/'cut_metadata.json');assert cut['selected_cut_has_no_legacy080_082_088_self']
    locks={**read(H/'inventory_registration.json')['source_locks'],**read(H/'cut_registration.json')['source_locks']}
    for p,h in locks.items():assert sha(p)==h
    rows=read(B/'upper_path112/nodes.json.gz');cert={r['key']:r for r in read(B/'upper_path112/node_certificates.json.gz')};geos=read(B/'upper_path112/geometries.json.gz')
    pre=read(S/'preflight002/registration.json');updates={}
    for i in range(15):
        w=read(S/f'run001/worker{i:02d}/result.json')
        for ch in w['chunks']:
            for r in read(ch['path']):updates[r['key']]=r
    for r in read(pre['reuse150']['path']):updates[r['key']]=r
    special='5:31,31,31,31';assert set(r['key'] for r in rows)-set(updates)=={special}
    r092={r['key']:r for p in sorted((B/'upper_range092').glob('block*_cells.json.gz')) for r in read(p)}
    input092={r['key']:r for r in read(B/'upper_range092/nodes.json.gz')}
    other={}
    for name in ('095','103'):
        other[name]={r['key']:r['bounds'] for r in read(B/f'upper_adaptive{name}/children_certificates.json.gz')}
    times=pre['times'];assert len(times)==73;counts=Counter();families=Counter()
    with (H/'node_index.jsonl.gz').open('xb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,filename='') as z:
            for ident,r in enumerate(r for r in sorted(rows,key=lambda r:r['key']) if r['key']!=special):
                k=r['key'];c,u=cert[k],updates[k];family='092' if '/upper_range092/' in r['packet_path'] else '095' if '/upper_adaptive095/' in r['packet_path'] else '103'
                assert c['new_selfU']<c['old_selfU'] and u['old_selfU']==c['new_selfU'] and u['new_mean_i']<c['new_M_i']
                if family=='092':
                    d,old=r092[k],input092[k]
                    assert r['old_L_i']==d['rangeL_i']>old['old_rangeL_i'] and r['old_U_i']==d['rangeU_i']<old['old_rangeU_i']
                else:
                    d=other[family][k];assert (r['old_L_i'],r['old_U_i'])==(d['L_i'],d['U_i'])
                shard=ident%15;counts[shard]+=1;families[family]+=1
                row={n:r[n] for n in ('key','level','indices','weight','geometry','point','packet_path','packet_sha256')}
                row.update(id=ident,shard=shard,range_source=family,L_i=r['old_L_i'],U_i=r['old_U_i'],m112_i=c['new_M_i'],self112U=c['new_selfU'],
                    new_mean_i=u['new_mean_i'],new_selfU=u['new_selfU'],h_infinity=geos[r['geometry']]['profile']['h'])
                z.write((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode())
    assert sum(counts.values())==97971
    cap=read(H.parent/'upper_join132/portable_components.json')['objects']['analytic_stationary']['M_upper_portable_cap'];assert F(cap)<3
    j=read(B/'upper_increment139/actual/results.json')['auxiliary_J_interval'][1];candidate=read(S/'run001/results.json')['new_UH']
    put(H/'certificate.json',{'schema':'I167H001','times':times,'model':'alpha1_body16_unnormalized_Lebesgue_Fourier2pi^-5','endpoint_packets':97971,
      'unused_legacy_self_key':special,'cap_domain':'receiver_complete_F_H512_only','M132_cap':cap,'qualified_J139_upper':j,'candidate_UH':candidate,
      'status':'prepared_candidate_payload_not_I156_acceptance','new_I_composition_evaluations':0,
      'numeric_H_interval_truth':'shared accepted old numerical dependency;projection hashes certify identity,not enclosure truth',
      'source_field_projection':{'point':'$.point or $.anchor','endpoints':'$.endpoints[j].H.lower/upper','times':'$.endpoints[j].t -> global times[j]',
        '092_ancestry':'$.original_raw_sha256','095_103_ancestry':'packet gzip SHA and source certificates'},
      'canonical_projection':'UTF8 compact sorted JSON of reduced-rational point and 73 {t,H:{lower,upper}};SHA256',
      'portable_verify_modes':'decode-only versus separately invoked reader arithmetic;preparation executes neither'})
    portable={}
    def add(rel,p):
        p=Path(p);portable[rel]={'path':str(p),'sha256':sha(p)};locks[str(p)]=sha(p)
    for n in ('codec.py','verify_packet.py','certificate.json','node_index.jsonl.gz'):add(n,H/n)
    add('tree.json.gz',B/'upper_path112/numeric_hierarchy.json.gz');add('geometries.json.gz',B/'upper_path112/geometries.json.gz')
    add('mean_factors.json',S/'run001/factors.json');add('factors112.json',B/'upper_path112/factor_catalog.json');add('ratios112.json',B/'upper_path112/ratio_catalog.json')
    for n in ('092','095','103'):add(f'range_{n}.json',B/('upper_range092' if n=='092' else f'upper_adaptive{n}')/'ratio_catalog.json')
    for n in ('pathmean.py','secant.py'):add('math/'+n,B/'upper_path112'/n)
    for n in ('report.md','check.py'):add('139/'+n,B/'upper_increment139'/n)
    for n in ('inputs.json','plan.json','registration.json','certificate.json','results.json'):add('139/actual/'+n,B/'upper_increment139/actual'/n)
    for directory,names in {
      'upper_alpha1_029':['report.md'],'upper_mean085':['report.md'],'upper_shift091':['report.md'],
      'upper_pathmean111':['report.md'],'audit_path121':['report.md','whole_protocol.md','whole_results.json'],
      'upper_rearrangement133':['report.md'],'upper_facepeak134':['report.md'],'upper_shorttime076':['report.md']}.items():
        for n in names:
            p=B/directory/n
            if not p.exists() and directory=='upper_facepeak134':p=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-auditcurvature129-sol/cube-dim-5-bounds/upper_facepeak134')/n
            add('proof_sources/'+directory+'/'+n,p)
    for n in ('portable_components.json','final_report.md','report.md','inputs.json','registration.json','results.json','budget.py','arith.py'):
        add('proof_sources/upper_join132/'+n,H.parent/'upper_join132'/n)
    # Preserve exact old source recipes,without their whole raw cache or ledgers.
    for d,n in [('upper_range092','engine.py'),('upper_adaptive095','kernel.py'),('upper_adaptive103','kernel.py'),('upper_path112','production_v2.py')]:add('proof_sources/'+d+'/'+n,B/d/n)
    for p in H.glob('*.py'):ast.parse(p.read_text())
    code={p.name:sha(p) for p in H.glob('*.py')}
    for n in ('inventory_registration.json','inventory.json','cut_registration.json','cut_metadata.json'):
        locks[str(H/n)]=sha(H/n)
    for p,h in locks.items():assert sha(p)==h
    plan={'task':'upper_integral167_fixed_exact_endpoint_export','status':'prepared_wait_primary_source_code_review_and_explicit_export_authorization',
      'endpoint_packets':97971,'endpoint_intervals':97971*73,'endpoint_rational_fields':97971*73*2,'excluded_unused_legacy_node':special,
      'shards':15,'shard_counts':[counts[i] for i in range(15)],'families':dict(families),'shard_hard_limit_bytes':90000000,'chunk_rows':256,
      'times':times,'index_sha256':sha(H/'node_index.jsonl.gz'),'index_bytes':(H/'node_index.jsonl.gz').stat().st_size,
      'code_locks':code,'source_locks':locks,'portable_files':portable,'portable_metadata_bytes':sum(Path(r['path']).stat().st_size for r in portable.values()),
      'export_numeric_evaluators':0,'decode_or_math_checks_during_preparation':0,'new_H_R_factor_mean_oldpower_calls':0,
      'source_interval_truth_boundary':'same accepted old H enclosures;not a new independent H evaluation certificate',
      'math_checker_execution':'prepared for readers only;not authorized or executed in this session','new_I_C5_claim':False,
      'independent161_and_primary_I156_acceptance':'pending;exported data may not change candidate qualification',
      'wall_cap':None,'fault':'preserve transactional completed chunks;stop;explicit same-plan resume only'}
    put(H/'export_plan.json',plan)
    print(json.dumps({'pid':os.getpid(),'plan_sha256':sha(H/'export_plan.json'),'input_packets':97971,'index_bytes':plan['index_bytes'],'source_locks':len(locks),'portable_metadata_bytes':plan['portable_metadata_bytes'],'new_H_R_oldpower_calls':0}))
if __name__=='__main__':main()
