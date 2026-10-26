#!/usr/bin/env python3
"""Only source qualification/sharding. No factors, means, gains or new UH."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,collections,gzip,json,os,resource,subprocess,sys,time
from fractions import Fraction as F
from common import sha,read,put,lines,endpoint_arrays,lockcheck

HERE=Path(__file__).resolve().parent

def inventory_worker(out,i):
    reg=read(out/'registration.json');idx=reg['shards'][i]
    assert sha(idx['path'])==idx['sha256']
    dest=out/f'qualification{i:02d}';dest.mkdir(exist_ok=False)
    put(dest/'registration.json',{'pid':os.getpid(),'shard':idx,'mode':'parse_identity_dyadic_band_completeness_only','new_factors_products_means_gains_UH':0})
    start=time.monotonic();count=0;bytes_=0;maxE=0;maxpacket=0;source=collections.Counter()
    for row in lines(idx['path']):
        p=Path(row['packet_path']);assert sha(p)==row['packet_sha256'],str(p)
        raw=read(p);lo,up,E=endpoint_arrays(raw,reg['times'],row['point'],row['key'])
        count+=1;bytes_+=p.stat().st_size;maxpacket=max(maxpacket,p.stat().st_size);maxE=max(maxE,E);source[row['source_group']]+=1
        if count%256==0:put(dest/'state.json',{'pid':os.getpid(),'qualified_prefix':count,'new_factor_products':0},True)
    assert count==idx['count']
    # Complete source lock after parsing, without any factor/product calculation.
    for row in lines(idx['path']):assert sha(row['packet_path'])==row['packet_sha256']
    result={'pid':os.getpid(),'count':count,'compressed_packet_bytes':bytes_,'largest_packet_bytes':maxpacket,'max_endpoint_binary_exponent':maxE,
        'source_groups':dict(source),'seconds':time.monotonic()-start,'maxrss_native_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'all_packet_hashes_before_after_unchanged':True,'new_factors_products_means_gains_UH':0}
    put(dest/'result.json',result)

def prepare(out):
    spec=read(HERE/'sources.json');lockcheck(spec['source_locks']);out.mkdir(exist_ok=False);start=time.monotonic()
    put(out/'initial_registration.json',{'task':'upper_integral156_source_inventory_only','pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
        'prepare_sha256':sha(__file__),'common_sha256':sha(HERE/'common.py'),'sources_sha256':sha(HERE/'sources.json'),'finite_original_nodes':97972,'workers':15,'new_factor_products_means_gains_UH':0})
    src=spec['source_locks'];rows=read(src['112_nodes']['path']);cert=read(src['112_certificates']['path']);tree=read(src['112_tree']['path'])
    cm={c['key']:c for c in cert};nodes=tree['nodes'];rm={r['key']:r for r in rows};roots=set(tree['roots'])
    assert len(rm)==len(cm)==len(nodes)==97972 and len(tree['roots'])==52360
    times=read(src['112_scalar_registration']['path'])['time_endpoints'];assert len(times)==73 and times[0]=='1/512'
    geos=read(src['112_geometries']['path']);fc=read(src['112_factors']['path']);rc=read(src['112_ratios']['path'])
    assert sorted(geos)==sorted(fc)==['g000','g001','g002','g003','g004']
    for g in geos:
        assert set(fc[g])==set(times)
        for t in times:assert str(F(geos[g]['profile']['h'])/F(t)) in rc
    # Existing150 values are copied, not reconstructed or recalculated.
    reuse=read(src['150_nodes']['path']);oldf=read(src['150_factors']['path']);r150=read(src['150_results']['path']);audit=read(src['152_result']['path'])
    assert len(reuse)==64 and len(oldf)==73 and r150['new_UH']==audit['new_UH']
    assert audit['status']=='finite73factor64root_arithmetic_independently_qualified'
    reused={r['key']:r for r in reuse};assert len(reused)==64
    eligible=[];unchanged=[];counts=collections.Counter();gcounts=collections.Counter();topology=[]
    for key in nodes:
        n=nodes[key];r=rm[key];c=cm[key]
        assert c['packet_sha256']==r['packet_sha256'] and n['level']==r['level'] and n['indices']==r['indices']
        assert n['selfU']==c['new_selfU'] and n['selfL']==r['old_selfL']
        assert r['old_L_i']<=c['new_M_i']<=r['old_U_i'] and c['old_selfU']==r['old_selfU']
        assert (n['parent'] is None)==(key in roots)
        for ch in n['children'] or []:assert nodes[ch]['parent']==key and nodes[ch]['level']==n['level']+1
        topology.append({**{k:n[k] for k in ['key','level','parent','selfL','selfU','effU']},'children':n['children'] or []})
        if key in reused:
            x=reused[key];assert n['parent'] is None and not n['children'] and n['effU']==n['selfU']==x['old_selfU']
            assert r['geometry']=='g002' and x['packet_sha256']==r['packet_sha256'];counts['reused150']+=1;continue
        reason=[]
        if not c['new_selfU']<c['old_selfU']:reason.append('112strict_chord_not_active')
        if not r['old_U_i']>=1<<144:reason.append('U_below_2^-16')
        if not r['old_L_i']<r['old_U_i']:reason.append('range_not_strict')
        if reason:unchanged.append({'key':key,'reasons':reason});counts['unchanged']+=1;continue
        group='range092' if '/upper_range092/' in r['packet_path'] else 'adaptive095' if '/upper_adaptive095/' in r['packet_path'] else 'adaptive103'
        assert f'/upper_{group}/' in r['packet_path'] if group!='range092' else '/upper_range092/' in r['packet_path']
        eligible.append({'key':key,'geometry':r['geometry'],'level':r['level'],'weight':r['weight'],'point':r['point'],
            'old_L_i':r['old_L_i'],'old_U_i':r['old_U_i'],'old_mean_i':c['new_M_i'],'old_selfL':n['selfL'],'old_selfU':n['selfU'],
            'packet_path':r['packet_path'],'packet_sha256':r['packet_sha256'],'source_group':group})
        gcounts[r['geometry']]+=1;counts['new_mean_nodes']+=1
    # Sort by source and file path, then distribute consecutive source blocks evenly.
    eligible.sort(key=lambda r:(r['source_group'],r['packet_path'],r['key']))
    shards=[[] for _ in range(15)]
    for group in sorted({r['source_group'] for r in eligible}):
        batch=[r for r in eligible if r['source_group']==group]
        for j,r in enumerate(batch):shards[j%15].append(r)
    indices=[]
    for i,batch in enumerate(shards):
        p=out/f'index{i:02d}.jsonl.gz'
        with gzip.open(p,'xt') as f:
            for r in batch:f.write(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
        indices.append({'worker':i,'path':str(p),'sha256':sha(p),'count':len(batch),'bytes':p.stat().st_size})
    # Topology is external only; each consumption worker never loads it.
    tp=out/'topology.json.gz'
    with gzip.open(tp,'xt') as f:json.dump({'roots':tree['roots'],'nodes':topology,'old_complete_UH':tree['budget']['UH']},f,separators=(',',':'))
    put(out/'reuse150.json',reuse);put(out/'unchanged.json',unchanged)
    reg={'task':'upper_integral156_source_inventory_only','pid':os.getpid(),'initial_registration_sha256':sha(out/'initial_registration.json'),
        'sources_sha256':sha(HERE/'sources.json'),'times':times,'shards':indices,'topology':{'path':str(tp),'sha256':sha(tp)},
        'reuse150':{'path':str(out/'reuse150.json'),'sha256':sha(out/'reuse150.json')},'unchanged':{'path':str(out/'unchanged.json'),'sha256':sha(out/'unchanged.json')},
        'counts':dict(counts),'new_mean_geometry_counts':dict(gcounts),'new_factors':292,'reused_factors':73,'new_factor_products_means_gains_UH':0}
    put(out/'registration.json',reg)
    procs=[]
    for i in range(15):
        log=(out/f'worker{i:02d}.log').open('xb');p=subprocess.Popen([sys.executable,'-B',__file__,'worker','--output',str(out),'--worker',str(i)],stdout=log,stderr=subprocess.STDOUT);log.close();procs.append(p)
    exits=[p.wait() for p in procs];put(out/'worker_exitcodes.json',{'pids':[p.pid for p in procs],'exitcodes':exits})
    assert exits==[0]*15,exits
    results=[read(out/f'qualification{i:02d}/result.json') for i in range(15)]
    assert sum(r['count'] for r in results)==len(eligible)
    lockcheck(spec['source_locks'])
    result={'task':'upper_integral156','status':'qualified_and_sharded_only_no_new_factor_product_mean_gain_or_UH',
        'pid':os.getpid(),'seconds':time.monotonic()-start,'counts':dict(counts),'new_mean_geometry_counts':dict(gcounts),'new_time_factors':292,'reused_g002_factors':73,
        'fixed_new_band_products':73*len(eligible),'source_original_nodes':len(nodes),'source_roots':len(tree['roots']),'times':times,'registration_sha256':sha(out/'registration.json'),
        'external_files':{p.name:{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(out.iterdir()) if p.is_file()},
        'workers':results,'all_source_hashes_before_after_unchanged':True,'compressed_packet_bytes':sum(x['compressed_packet_bytes'] for x in results),
        'largest_packet_bytes':max(x['largest_packet_bytes'] for x in results),'max_endpoint_binary_exponent':max(x['max_endpoint_binary_exponent'] for x in results),
        'new_factors_products_means_gains_UH':0,'global_Cn_claim':False}
    put(out/'results.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['external_files','workers','times']},sort_keys=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','worker']);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--worker',type=int);a=ap.parse_args()
    if a.mode=='prepare':prepare(a.output)
    else:inventory_worker(a.output,a.worker)
