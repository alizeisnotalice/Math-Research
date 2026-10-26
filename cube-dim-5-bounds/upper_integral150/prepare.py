#!/usr/bin/env python3
"""Read-only source inventory and old-field ranking only. No new mean/H/R.

Registers before selecting64 current physical leaves. Filters ensure the
saved112 self upper came from its new convex-chord branch; no old powers
or old range/mean proofs are re-evaluated.
"""
from datetime import datetime,timezone
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,os,time

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
P=MAIN/'upper_path112'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p);return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def write(n,x):
    with (HERE/n).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def main():
    assert not (HERE/'selected64.json').exists(),'Preserve this finite preparation'
    names=['nodes.json.gz','node_certificates.json.gz','numeric_hierarchy.json.gz','geometries.json.gz','factor_catalog.json','ratio_catalog.json','scalar_registration.json','manifest.json']
    locks={n:{'path':str(P/n),'sha256':sha(P/n)} for n in names}
    write('preflight_registration.json',{'task':'upper_integral150_finite_preparation','pid':os.getpid(),
        'utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),'sources':locks,
        'finite_source_nodes':97972,'select_current_physical_leaves':64,
        'ranking':'saved weighted112 selfU-selfL decreasing,thenlevel/indices',
        'eligibility':'unexpanded physical leaf;112new_selfU<old_selfU;current self=112newself;U>=1/65536',
        'new_H_R_gradient_log_exp_scalar_factor_or_mean_calls':0,'old_power_baseline_runs':0})
    began=time.monotonic();rows=read(P/'nodes.json.gz');certs=read(P/'node_certificates.json.gz');tree=read(P/'numeric_hierarchy.json.gz')
    assert len(rows)==len(certs)==len(tree['nodes'])==97972 and len(tree['frontier'])==94580
    rm={r['key']:r for r in rows};cm={r['key']:r for r in certs};nodes=tree['nodes'];eligible=[]
    for key in tree['frontier']:
        row,c,n=rm[key],cm[key],nodes[key]
        assert not n['children'] and n['selfU']==c['new_selfU'] and n['selfL']==row['old_selfL']
        if c['new_selfU']<c['old_selfU'] and row['old_U_i']>=1<<144 and row['old_L_i']<row['old_U_i']:
            eligible.append(key)
    keys=sorted(eligible,key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:64]
    assert len(keys)==64
    factors=read(P/'factor_catalog.json');geos=read(P/'geometries.json.gz');cat=read(P/'ratio_catalog.json')
    times=read(P/'scalar_registration.json')['time_endpoints'];assert len(times)==73 and times[0]=='1/512'
    records=[];gids=set();packets={}
    for rank,key in enumerate(keys):
        row,c,n=rm[key],cm[key],nodes[key];packet=Path(row['packet_path']);assert sha(packet)==row['packet_sha256']
        raw=read(packet);assert raw['key']==key and raw['point']==row['point'] and [x['t'] for x in raw['endpoints']]==times
        gids.add(row['geometry']);packets[key]={'path':str(packet),'sha256':row['packet_sha256']}
        records.append({'rank':rank,'key':key,'source_node':row,'source112_certificate':c,
            'saved_selfL':n['selfL'],'saved_selfU':n['selfU'],'old_weighted_gap_units':n['selfU']-n['selfL'],
            'U_ge_1over65536':True,'new112_chord_branch_strictly_active':True})
    needed=set()
    for gid in gids:
        assert set(factors[gid])==set(times)
        h=F(geos[gid]['profile']['h'])
        for t in times:
            x=str(h/F(t));assert x in cat;needed.add(x)
    result={'task':'upper_integral150','status':'fixed64_source_preparation_complete_no_new_mean',
        'pid':os.getpid(),'seconds':time.monotonic()-began,'eligible_current_physical_leaves':len(eligible),
        'selected':64,'geometry_ids':sorted(gids),'time_endpoints':times,'band_count':73,
        'finite_selected_band_products':4672,'maximum_geometry_band_factors':len(gids)*73,
        'all_range_ratio_endpoint_keys_already_cached':True,'used_existing_range_ratio_keys':len(needed),
        'new_R_H_gradient_log_exp_factor_or_mean_calls':0,'old_power_baseline_runs':0,
        'source_locks':locks,'selected_packet_locks':packets,
        'saved_complete_UH':tree['budget']['UH'],'new_UH_claim':False}
    write('selected64.json',records);result['selected64_sha256']=sha(HERE/'selected64.json')
    write('preflight_results.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['source_locks','selected_packet_locks','time_endpoints']},sort_keys=True))
if __name__=='__main__':main()
