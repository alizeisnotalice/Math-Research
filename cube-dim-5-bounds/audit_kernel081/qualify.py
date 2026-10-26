#!/usr/bin/env python3
"""Join completed independent values, closed tree, and locked075/029 theory."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import json,os,time,resource
from common import require,file_sha,digest

HERE=Path(__file__).resolve().parent


def main():
    started=time.monotonic();pid=os.getpid()
    qr=json.loads((HERE/'qualification_registration.json').read_text())
    for path,h in qr['input_locks'].items():require(file_sha(path)==h,'qualification input lock mismatch: '+path)
    reg=json.loads((HERE/'registration.json').read_text());nr=json.loads((HERE/'numeric_result.json').read_text());sr=json.loads((HERE/'structural_result.json').read_text())
    require(nr['status']=='numeric_complete_all1318_H_plus4grad' and nr['verified_count']==len(nr['points'])==1318,'new values are not all independently verified')
    require(sr['status']=='structural_and_all16corner_passed_pending_values' and sr['new_certificates_verified']==1318 and sr['new_corner_points_verified']==21088,'structural/certificate proof incomplete')
    for name,h in nr['code_locks'].items():require(file_sha(HERE/name)==h,'numeric successful-prefix implementation changed')
    for name,h in sr['code_locks'].items():require(file_sha(HERE/name)==h,'structural implementation changed')
    require(nr['shared_production_code'] is False and nr['shared_Arb_or_log_library'] is False,'new value path is not numerically independent')
    manifest=json.loads(Path(reg['manifest']).read_text());jobs=[]
    for item in manifest['batches']:
        p=Path(item['directory'])/'point_manifest.json';require(file_sha(p)==item['files']['point_manifest.json'],'point manifest changed')
        jobs.extend({'batch':item['batch'],'ordinal':ordinal,**entry} for ordinal,entry in enumerate(json.loads(p.read_text())['point_checkpoints']))
    last=manifest['batches'][-1];overlay_path=Path(last['directory'])/'state.json';require(file_sha(overlay_path)==last['files']['state.json'],'last overlay changed')
    cache=json.loads(overlay_path.read_text())['new_cache'];require(len(cache)==len(jobs)==1318,'new cache/catalogue count differs')
    minimum=None;minimum_job=None;checked_components=0
    for index,(point,job) in enumerate(zip(nr['points'],jobs)):
        require(point['index']==index and point['batch']==job['batch'] and point['ordinal']==job['ordinal'] and point['point']==job['point'] and point['point_checkpoint_sha256']==job['sha256'],'verified numeric prefix identity differs')
        require(file_sha(job['path'])==job['sha256'],'raw successful-point checkpoint changed')
        raw=json.loads(Path(job['path']).read_text());key=','.join(job['point']);record=cache[key]
        require(raw['record']==record and raw['stage']=='complete' and raw['partial'] is None,'raw numeric target differs from new cache')
        values=point['independent_values_512_dyadic'];require(len(values)==len(record['values'])==5,'not all5 components were independently computed')
        for independent,production in zip(values,record['values']):
            a,b=map(F,independent);L,U=map(F,production)
            require(L<a<=b<U,'independent saved enclosure is not strictly inside original ball');checked_components+=1
        require(F(values[0][0])>F(1,10000),'independent anchor H lower is not above target')
        h=F(record['values'][0][0])
        if minimum is None or h<minimum:minimum=h;minimum_job={'index':index,'batch':job['batch'],'ordinal':job['ordinal'],'point':job['point'],'production_H_lower':str(h),'independent_H_and_all4_gradients':values}
    counts=sr['final_counts'];require(counts=={'passed':97960,'splits':45600,'queue':0,'cache':180040,'roots':52360},'closed tree final counts differ')
    require(sr['batches'][-1]['tree']['every_root_Kraft_sum_exactly_one'] is True and sr['batches'][-1]['tree']['leaf_count']==counts['passed'],'final tree is not completely certified')
    require(sr['materialization']['direct_object_comparisons']['old_cache']==178722 and sr['materialization']['direct_object_comparisons']['new_cache']==1318,'source/new cache direct comparison incomplete')
    require(all(F(x['minimum_log_margin'])>0 for x in sr['new_certificates'] if x['passed']),'a new accepted closed-box margin is not strict')
    require(sum(x['passed'] for x in sr['new_certificates'])==1001,'new accepted-leaf count differs')
    # Generators of every signed coordinate permutation preserve the16 body lines.
    V={tuple([1,*tail]) for tail in product((-1,1),repeat=4)}
    def representative(vector):return tuple(vector if vector[0]>0 else [-v for v in vector])
    generators=[]
    for coordinate in range(5):
        transformed=set()
        for v in V:
            w=list(v);w[coordinate]*=-1;transformed.add(representative(w))
        require(transformed==V,'body process signed-coordinate symmetry failed');generators.append('sign_flip_'+str(coordinate))
    for coordinate in range(4):
        transformed=set()
        for v in V:
            w=list(v);w[coordinate],w[coordinate+1]=w[coordinate+1],w[coordinate];transformed.add(representative(w))
        require(transformed==V,'body process permutation symmetry failed');generators.append('adjacent_swap_'+str(coordinate))
    for path,h in reg['protected_old_locks'].items():require(file_sha(path)==h,'old source proof lock changed')
    for path,h in qr['input_locks'].items():require(file_sha(path)==h,'qualification input changed during finite check')
    result={'task':'audit_kernel081_qualified_join','pid':pid,'status':'all_values_structure_and_closed_cube_g_qualified_given_locked075_029',
            'qualification_registration_sha256':file_sha(HERE/'qualification_registration.json'),'code_sha256':file_sha(Path(__file__)),
            'input_locks':qr['input_locks'],'final_state_sha256':reg['final_state_sha256'],'numeric_points':1318,'independent_components':checked_components,
            'new_H_and_all4_physical_gradients_independently_enclosed':True,'new_cones_evaluated':1318*3504,'new_corner_points_verified':21088,
            'old_cache178722_recomputed':False,'old_local_leaf_qualification':'inherited unchanged from locked078/069/062 chain',
            'final_counts':counts,'minimum_new_anchor_H_candidate':minimum_job,'body_line_symmetry_generators_verified':generators,
            'closed_ordered_face_H_lower':'H_(1/64)(1,y)>1/10000 for1>=y2>=y3>=y4>=y5>=0',
            'all_cube_faces_G_lower':'G(y)>1/10000 for||y||infinity=1; includes extended infinite singular values',
            'global_G_lower':'G(x)>(1/10000)*||x||infinity^(-4) for x!=0, by full Green homogeneity only',
            'origin':'canonical G(0)=infinity','essential_closed_unit_cube_g_lower':'1/10000',
            'H_epsilon_homogeneity_used':False,'same_alpha1_L1_all_positive_scale_qualification':'locked029; no alpha1.5 substitution',
            'Cn_or_C5_updated':False,'I_updated_or_recomputed':False,'main_retains_final_acceptance':True,
            'analytical_dependency':'true017 cone representation/normalization and locked075/029 measure-theoretic proofs are inherited, not formally reproved by integer checks',
            'wall_seconds':time.monotonic()-started,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'MATLAB_calls':0,'new_H_evaluations_in_this_join':0}
    output=HERE/'qualification_result.json';temp=HERE/'qualification_result.json.tmp';temp.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');temp.replace(output)
    print(json.dumps({k:result[k] for k in ('pid','status','numeric_points','independent_components','final_counts','wall_seconds','peak_rss_bytes')}),flush=True)
    print(json.dumps({'qualification_result_sha256':file_sha(output)}),flush=True)


if __name__=='__main__':main()
