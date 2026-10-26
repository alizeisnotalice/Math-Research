#!/usr/bin/env python3
"""12984 distinct single-anchor whole-box checks; no H/log module imported."""
import sys
sys.dont_write_bytecode=True
import json,os,time,traceback
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import lcm
from interfaces import REG,HERE,common as C


def write(path,value):
    p=Path(path);q=Path(str(p)+'.tmp')
    with q.open('wb') as f:f.write(json.dumps(value,sort_keys=True,indent=2).encode()+b'\n');f.flush();os.fsync(f.fileno())
    q.replace(p)


def rows(path):
    with Path(path).open() as f:
        for line in f:yield json.loads(line,object_pairs_hook=C.no_duplicates)


def candidate(proof,lo,hi):
    cert=proof['source_certificate'];epsilon=F(proof['epsilon']);penalty=F(cert['penalty_coefficient']);C.require(penalty==1/(16*epsilon**2),'mixed epsilon/curvature')
    delta=F(cert['Delta_lower']);gradients=[list(map(F,x)) for x in cert['loggradient_directed']];point=list(map(F,proof['point']));inside=all(a<=x<=b for a,x,b in zip(lo,point[1:],hi))
    C.require(F(proof['H_interval'][0])>0,'source H lower not positive')
    offsets=[(a-x,b-x) for a,x,b in zip(lo,point[1:],hi)];vd=lcm(*(x.denominator for pair in offsets for x in pair));gd=lcm(delta.denominator,*(x.denominator for pair in gradients for x in pair))
    zz=[tuple(int(x*vd) for x in pair) for pair in offsets];gg=[tuple(int(x*gd) for x in pair) for pair in gradients];di=int(delta*gd);allcorners=[]
    for bits in product((0,1),repeat=4):
        v=[zz[j][bits[j]] for j in range(4)];num=di*vd*vd+vd*sum(min(a*x,b*x) for (a,b),x in zip(gg,v))-int(penalty)*gd*sum(x*x for x in v);allcorners.append((num,bits))
    minimum,bits=min(allcorners);margin=F(minimum,gd*vd*vd)
    return {'proof_source_id':proof['id'],'epsilon':proof['epsilon'],'pair_cache_key':proof['pair_cache_key'],'minimum_log_margin':str(margin),'worst_corner_bits':list(bits),'checked_points':16,'passed':minimum>0,'anchor_inside_child':inside,'saved_Delta_negative':delta<0,'anchor_HL_above_target':F(proof['H_interval'][0])>F(11,100000)}


def source_gate():
    root=Path(REG['reuse_directory']);source107=json.loads(Path(REG['source107']).read_text());proofs=json.loads((root/'proof_sources.json').read_text());inputs=json.loads((root/'child_inputs.json').read_text());oldrows={C.node(row['node']):row for row in rows(Path(REG['source105']).parent/'pending_gap_order.jsonl')}
    verified64={}
    for p in sorted(Path(REG['source108_directory']).glob('block*.json')):
        data=json.loads(p.read_text());C.require(data['status']=='complete','108 source block not qualified')
        for row in data['points']:verified64[row['source_leaf_index']]=row
    receipt116=json.loads((Path(REG['source116_directory'])/'structure_result.json').read_text());C.require(receipt116['status']=='first4096_H128_values_source_journal_and_partial_tree_independently_qualified','116 source gate not qualified')
    verified128={C.node(x['node']):x for x in receipt116['certificates']};lookup={};expected_inputs=[];different=0
    C.require(len(source107['new_splits'])==3246 and len(source107['new_passed'])==850 and len(proofs)==6492 and len(inputs)==6492,'finite input counts differ')
    for index,split in enumerate(source107['new_splits']):
        parent=C.node(split['node']);cert128=split['failed_certificate'];record128=source107['new_cache_by_epsilon']['1/128'][cert128['point_key']];old=oldrows[parent];cert64=old['certificate'];certificate_data=[('1/128',cert128,record128['point'],record128['values'][0],C.digest(record128)),('1/64',cert64,old['pair_cache_key'][1].split(','),old['H_interval'],old['source_record_sha256'])]
        qualified128=verified128[parent];qualified64=verified64[old['source_leaf_index']]
        C.require(not qualified128['passed'] and qualified128['certificate_sha256']==C.digest(cert128) and qualified128['record_sha256']==C.digest(record128),'128 source not bound to independent116')
        C.require(not qualified64['passed'] and qualified64['row_sha256']==C.digest({k:v for k,v in old.items() if k!='gap_rank'}) and qualified64['certificate_sha256']==C.digest(cert64) and qualified64['epsilon']=='1/64','64 source not bound to independent108')
        different+=record128['point']!=old['pair_cache_key'][1].split(',')
        for j,(epsilon,cert,point,H,recordsha) in enumerate(certificate_data):
            proof=proofs[2*index+j];C.require(cert['target']=='11/100000' and cert['leaf_epsilon']==epsilon and cert['passed'] is False,'source parent target/epsilon/classification differs')
            expected={'id':str(index)+':'+epsilon,'parent':list(parent),'epsilon':epsilon,'point':point,'pair_cache_key':[epsilon,','.join(point)],'source_record_sha256':recordsha,'source_certificate_sha256':C.digest(cert),'source_certificate':cert,'H_interval':H,'HL_positive':F(H[0])>0,'anchor_HL_above_target':F(H[0])>F(11,100000),'saved_Delta_negative':F(cert['Delta_lower'])<0}
            C.require({k:v for k,v in proof.items() if k!='Delta_and_gradient_qualification'}==expected and expected['HL_positive'],'saved point/Delta/gradient/record binding differs')
            C.require(proof['id'] not in lookup,'duplicate proof id');lookup[proof['id']]=proof
        for child in split['children']:expected_inputs.append({'child_index':len(expected_inputs),'node':child,'parent':list(parent),'proof_sources':[str(index)+':1/128',str(index)+':1/64']})
    reader=C.Stream(REG['source098']);roots=None;splits=[]
    while True:
        field=reader.key()
        if field is None:break
        if field=='cache_by_epsilon':
            reader.char('{')
            for j,epsilon in enumerate(('1/128','1/64')):
                if j:reader.char(',')
                C.require(reader.value()==epsilon,'source cache epsilon differs');reader.char(':')
                for _ in reader.items(True):pass
            reader.char('}')
        elif field=='roots':roots=reader.value()
        elif field=='splits':
            for row in reader.items():
                node=C.node(row['node']);C.require(row['axis']==len(node[1])%4 and row['children']==[[node[0],node[1]+'0'],[node[0],node[1]+'1']],'old closed split differs');splits.append(node)
        elif field=='passed':
            for _ in reader.items():pass
        else:reader.value()
    reader.close();splits.extend(C.node(x['node']) for x in source107['new_splits'])
    for i,(inp,expected) in enumerate(zip(inputs,expected_inputs)):
        lo,hi=C.physical_box(roots,C.node(expected['node']));expected.update(box_lower=list(map(str,lo)),box_upper=list(map(str,hi)));C.require(inp==expected,'child input physical box/source ids differ')
    accepted=[];reader=C.Stream(REG['source105']);meta={}
    while True:
        field=reader.key()
        if field is None:break
        if field=='passed':
            for row in reader.items():accepted.append(C.node(row['node']))
        else:meta[field]=reader.value()
    reader.close();accepted.extend(C.node(x['node']) for x in source107['new_passed']);old877=source107['queue'][:877]
    C.require(old877==meta['queue'][4096:]==json.loads((root/'old877_unchanged.json').read_text()) and len(accepted)==96875 and len(splits)==51884,'source passed/splits/877 differ')
    C.require(source107['queue'][877:]==[x['node'] for x in inputs],'6492child input is not exact source107 tail')
    tree=C.validate_tree(roots,accepted+list(map(C.node,source107['queue'])),set(splits),'reuse118Source107Partial')
    return lookup,inputs,roots,splits,accepted,old877,tree,different


def main():
    began=time.monotonic();os.nice(5);rp=HERE/'audit_registration.json';registration=json.loads(rp.read_text())
    for p,h in {**REG['input_locks'],**REG['protected_old_locks'],**registration['code_locks']}.items():C.require(C.file_sha(p)==h,'frozen input changed: '+p)
    lookup,inputs,roots,splits,accepted,old877,source_tree,different=source_gate();root=Path(REG['reuse_directory']);author=list(rows(root/'child_results.jsonl'));prefix_path=HERE/'successful_prefix.json'
    prefix=json.loads(prefix_path.read_text()) if prefix_path.exists() else {'registration_sha256':C.file_sha(rp),'status':'verified_prefix_only','results':[]}
    C.require(prefix['registration_sha256']==C.file_sha(rp) and len(author)==6492,'prefix/input binding differs')
    try:
        for i,record in enumerate(prefix['results']):C.require(record==author[i],'durable mathematical prefix changed')
        for i in range(len(prefix['results']),6492):
            inp=inputs[i];lo,hi=C.physical_box(roots,C.node(inp['node']));candidates=[candidate(lookup[sourceid],lo,hi) for sourceid in inp['proof_sources']]
            winner=next((x for x in candidates if x['passed']),None);expected={'child_index':i,'node':inp['node'],'parent':inp['parent'],'candidates':candidates,'winner_proof_source_id':winner['proof_source_id'] if winner else None}
            C.require(expected==author[i],'two independent whole-box candidates/winner differ');prefix['results'].append(expected)
            if len(prefix['results'])%256==0:write(prefix_path,prefix)
        prefix['status']='all6492children_both_wholebox_candidates_complete';write(prefix_path,prefix)
    except BaseException:
        prefix.update(status='fault_preserved_prefix',fault=traceback.format_exc());write(prefix_path,prefix);raise
    newpassed=[];pending=[];best_rows=[];winning={'1/64':0,'1/128':0};outside={'1/64':0,'1/128':0};candidate_pass={'1/64':0,'1/128':0};candidate_out={'1/64':0,'1/128':0};negative={'1/64':0,'1/128':0};both=0
    for inp,result in zip(inputs,prefix['results']):
        candidates=result['candidates'];both+=sum(x['passed'] for x in candidates)>1
        for test in candidates:candidate_pass[test['epsilon']]+=test['passed'];candidate_out[test['epsilon']]+=not test['anchor_inside_child'];negative[test['epsilon']]+=test['saved_Delta_negative']
        winner=next((x for x in candidates if x['passed']),None)
        if winner:
            proof=lookup[winner['proof_source_id']];original=proof['source_certificate'];epsilon=proof['epsilon'];winning[epsilon]+=1;outside[epsilon]+=not winner['anchor_inside_child']
            cert={'target':'11/100000','lower_bound_object':'G','component_function':'H_epsilon','leaf_epsilon':epsilon,'mechanism':'075_saved_parent_anchor_16','proof_source_id':proof['id'],'pair_cache_key':proof['pair_cache_key'],'source_record_sha256':proof['source_record_sha256'],'source_parent_certificate_sha256':proof['source_certificate_sha256'],'Delta_lower':original['Delta_lower'],'loggradient_directed':original['loggradient_directed'],'penalty_coefficient':original['penalty_coefficient'],'minimum_log_margin':winner['minimum_log_margin'],'worst_corner_bits':winner['worst_corner_bits'],'checked_points':16,'passed':True,'anchor_inside_child':winner['anchor_inside_child'],'outside_anchor_allowed_by':'global075/099;g certificate only'}
            newpassed.append({'node':inp['node'],'leaf_epsilon':epsilon,'target':'11/100000','lower_bound_object':'G','certificate':cert})
        else:
            pending.append(inp['node']);best=max(candidates,key=lambda x:F(x['minimum_log_margin']))
            best_rows.append({'child_index':inp['child_index'],'node':inp['node'],'parent':inp['parent'],'box_lower':inp['box_lower'],'box_upper':inp['box_upper'],'kind':'both_saved_parent_anchor_certificates_unproved','best_single_anchor_margin':best['minimum_log_margin'],'best_single_anchor_proof_source_id':best['proof_source_id'],'candidate_results':candidates,'ranking_only':'max of separate wholebox minima;no cornerwise anchor mixing'})
    state=json.loads((root/'state.json').read_text());C.require(state['source105_sha256']==REG['source105_sha256'] and state['source107_sha256']==REG['source107_sha256'] and state['target']=='11/100000','state source/target differs')
    C.require(state['new_passed']==newpassed and state['remaining_children']==pending and state['remaining_original877']==old877 and state['queue']==old877+pending and state['new_splits']==[] and state['inflight']=={} and state['new_H_gradient_kernel_log_calls']==0 and state['global_g110_qualified'] is False,'exact newpass/retirement/877/tree state differs')
    C.require(state['proof_sources_sha256']==C.file_sha(root/'proof_sources.json') and state['child_results_sha256']==C.file_sha(root/'child_results.jsonl') and state['retained_passed_source105']==96025 and state['retained_passed_source107']==850,'source proofs/preserved acceptance differs')
    source_mapper=list(rows(Path(REG['source107']).parent.parent/'pending_mapper_batch0001.jsonl'));expected_mapper=[{**x,'pending118_index':i} for i,x in enumerate(source_mapper[:877])]
    expected_mapper.extend({**row,'pending118_index':877+i} for i,row in enumerate(best_rows));C.require(expected_mapper==list(rows(root/'pending_mapper.jsonl')),'full pending mapper or old877 metadata differs')
    ranked=list(rows(root/'remaining_children_gap_order.jsonl'));byindex={x['child_index']:x for x in best_rows};C.require(len(ranked)==len(best_rows)==5360,'remaining candidate group differs')
    last=None
    for rank,row in enumerate(ranked):
        C.require(row.pop('gap_rank')==rank and row==byindex.pop(row['child_index']),'ranked whole-box priority data differs');m=F(row['best_single_anchor_margin']);C.require(last is None or last<=m,'ranked margins not increasing');last=m
    C.require(not byindex and len(newpassed)==1132 and winning=={'1/64':1080,'1/128':52} and outside=={'1/64':996,'1/128':0} and both==0,'winner/epsilon/outside groups differ')
    tree=C.validate_tree(roots,accepted+[C.node(x['node']) for x in newpassed]+list(map(C.node,state['queue'])),set(splits),'reuse118FinalPartial')
    pids=[json.loads((root/'registration.json').read_text())['pid'],json.loads((root/'author_verification.json').read_text())['pid']]
    for pid in pids:
        try:os.kill(pid,0)
        except ProcessLookupError:continue
        raise AssertionError('118 process still present or reused')
    for p,h in {**REG['input_locks'],**REG['protected_old_locks'],**registration['code_locks']}.items():C.require(C.file_sha(p)==h,'input changed after finite zeroH audit')
    result={'task':'audit_reuse122','status':'all6492children_12984_single_anchor_candidates_independently_qualified_partial','pid':os.getpid(),'wall_seconds':time.monotonic()-began,'registration_sha256':C.file_sha(rp),'finite_children':6492,'single_anchor_candidates':12984,'integer_closed_corners':207744,'new_passed':1132,'new_passed_by_epsilon':winning,'successful_outside_anchors':outside,'candidate_pass_counts':candidate_pass,'candidate_outside_counts':candidate_out,'negative_Delta_candidates':negative,'child_both_candidates_pass':both,'different_H64_H128_parent_points':different,'passed':98007,'splits':51884,'queue':6237,'old877_verbatim':True,'remaining_children':5360,'tree':tree,'source_tree':source_tree,'new_H_gradient_kernel_log_calls':0,'new_splits':0,'source_Delta_grad_penalty_fields_unchanged':True,'source_qualification_from108_and116_independent_receipts':True,'new_leaf_single_anchor_all16':True,'anchor_inside_child_not_required':True,'HL_above_target_not_required':True,'negative_Delta_allowed':True,'no_cornerwise_anchor_mixing':True,'actual_118_PIDs_absent':pids,'protected_old_files_unchanged':250,'global_g110_qualified':False,'global_S_credit':0,'Cn_C5_I_updated':False}
    p=HERE/'result.json';write(p,result);print(json.dumps(result),flush=True);print(json.dumps({'result_sha256':C.file_sha(p)}),flush=True)


if __name__=='__main__':main()
