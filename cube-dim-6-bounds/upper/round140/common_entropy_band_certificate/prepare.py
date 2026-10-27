#!/usr/bin/env python3
"""Read-only C134 selection and exact-hash preparation; zero scientific logs."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,itertools,json
P=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round134/fixed_atlas_local_bisection')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(x):return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def key(k,p):return hashlib.sha256(canon({'semantic_key':k,'point':p})).hexdigest()
def points(box):return [['1']+[box[i][b]for i,b in enumerate(bits)]for bits in itertools.product((0,1),repeat=5)]
def write(name,x):(P/name).write_text(json.dumps(x,indent=2)+'\n')
def main():
    result=read(SRC/'result.json');cfg=read(SRC/'fixed_roots_and_atlas.json');cache=read(SRC/'point_cache.json');exit=read(SRC/'exit_receipt.json')
    assert result['status']==exit['status']=='FINITE_FIXED_BATCH_COMPLETE'
    assert result['target']==cfg['target']=='1/63500'
    assert exit['result_sha256']==sha(SRC/'result.json') and exit['point_cache_sha256']==sha(SRC/'point_cache.json')
    assert len(cfg['roots'])==7 and len(cache['entries'])==38055
    assert cache['initial_entry_count']==16929 and cache['actual_counts']['new_exp']==21126
    assert cache['old_corner_exp_recomputed']==0 and result['native_calls']==0
    frozen=[SRC/n for n in ['result.json','exit_receipt.json','fixed_roots_and_atlas.json','point_cache.json','fixed_source_entries.json','input_hashes.json','PROOF.md','run_batch.py']]
    selected=[];paid=0;total=0;used=set();partition_refs=[]
    for root in cfg['roots']:
        assert [len(b['candidate_keys'])for b in root['bands']]==[16,16,0,0,0]
        rp=SRC/'roots'/(root['original_node_id']+'.json');frozen.append(rp);tree=read(rp)
        assert tree['original_node_id']==root['original_node_id'] and tree['original_box']==root['box']
        assert tree['closed_partition_volume_exact']==tree['original_volume_exact']
        partition_refs.append({'root_id':root['original_node_id'],'path':str(rp),'sha256':sha(rp),'closed_partition_leaf_ids':tree['closed_partition_leaf_ids']})
        for nid in tree['closed_partition_leaf_ids']:
            total+=1;node=tree['nodes'][nid]
            if node['status']=='LOCAL_LEAF_PASS':paid+=1;continue
            assert node['status']=='DEPTH6_UNRESOLVED' and node['depth']==6
            row={'root_id':root['original_node_id'],'leaf_id':nid,'box':node['box'],'source_keys':[b['candidate_keys']for b in root['bands'][:2]],'windows':[b['window']for b in root['bands'][:2]]}
            for keys in row['source_keys']:
                assert len(set(keys))==16
                for k in keys:
                    for point in points(node['box']):
                        ck=key(k,point);r=cache['entries'][ck]
                        assert r['semantic_key']==k and r['point']==point and F(0)<F(r['lower'])<=F(r['upper'])
                        used.add(ck)
            selected.append(row)
    assert len(selected)==11 and total==73 and paid==62 and len(used)==8544
    write('selected_metadata.json',{'round':'140','target':'1/63500','selected_leaf_count':11,'selected_leaves':selected,'old_paid_leaf_count_reference_only':62,'C134_total_leaves':73,'original_closed_partition_references':partition_refs,'unique_required_source_point_entries':len(used),'atlas_enlarged':False})
    write('input_hashes.json',{'frozen_sources':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}for p in frozen],'source_directory_allowlist':str(SRC),'read_live129':False,'source_134_independent_acceptance':'parent assignment states accepted; formal approval must be reconfirmed in dual release'})
    write('prepared_packet.json',{'status':'PREPARED_ONLY_NO_LOG_ENTROPY_SCIENCE','round':'140','target':'1/63500','leaves':11,'pairs_per_leaf':256,'total_pairs':2816,'source_points_per_leaf':1024,'unique_required_source_point_entries':len(used),'max_Arb_source_logs':8544,'max_Arb_entropy_logs':5632,'max_Arb_logs_total':14176,'max_Arb_corner_score_combinations':90112,'float_objective_calls_per_pair':70,'max_float_objective_calls':197120,'weight_denominator':2**24,'precision_bits':512,'workers':1,'nice_minimum':10,'time_limit':None,'native_calls':0,'heat_kernel_calls':0,'source_point_exp_calls':0,'final_exp_calls':0,'old_paid_leaves_retried':0,'global_g_certificate':False,'new_C6_bound':False,'run_batch_sha256':sha(P/'run_batch.py'),'proof_sha256':sha(P/'PROOF.md')})
    import importlib.util
    spec=importlib.util.spec_from_file_location('c140_prepared_runner',P/'run_batch.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    write('release_requirements.json',mod.release_requirements())
    print(json.dumps({'status':'PREPARED_ONLY','selected_leaves':len(selected),'paid_leaves_reference_only':paid,'unique_required_source_points':len(used),'new_log_calls':0,'new_exp_calls':0}))
if __name__=='__main__':main()
