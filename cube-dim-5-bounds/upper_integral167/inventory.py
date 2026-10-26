#!/usr/bin/env python3
"""One finite metadata/byte-layout inventory. No H or envelope evaluator imports.

Reads old metadata and at most ten fixed representative endpoint packets.
Stats all recorded packets; does not bulk read/hash/export their contents.
"""
from pathlib import Path
from collections import Counter,defaultdict
import gzip,hashlib,json,os,struct,zlib

HERE=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds')
P156=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-curvature128-sol/cube-dim-5-bounds/upper_integral156')
A161=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-audit136-sol/cube-dim-5-bounds/audit_integral161')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def write(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n')
def family(p):
    p=str(p)
    return next(n for n in ('upper_range092','upper_adaptive095','upper_adaptive103') if '/'+n+'/' in p)
def dyad(s):
    a=s.split('/');n=int(a[0]);d=1 if len(a)==1 else int(a[1])
    assert len(a)<=2 and d>0 and d&(d-1)==0
    return n,d.bit_length()-1
def triple(lo,hi):
    n,e=dyad(lo);m,f=dyad(hi);E=max(e,f);L=n<<(E-e);U=m<<(E-f)
    assert 0<L<=U
    return L,U-L,E
def uint(n):
    assert n>=0
    b=n.to_bytes(max(1,(n.bit_length()+7)//8),'big')
    assert len(b)<65536
    return struct.pack('>H',len(b))+b
def payload(raw):
    # Prototype only: lower integer,nonnegative width,common binary exponent.
    out=bytearray();widths=[];exps=[];lowerbits=[]
    for e in raw['endpoints']:
        L,W,E=triple(e['H']['lower'],e['H']['upper'])
        out.extend(struct.pack('>i',E));out.extend(uint(L));out.extend(uint(W))
        widths.append(W.bit_length());exps.append(E);lowerbits.append(L.bit_length())
    return bytes(out),widths,exps,lowerbits
def main():
    assert not (HERE/'inventory_registration.json').exists()
    names=['upper_range092/nodes.json.gz','upper_range092/engine.py','upper_range092/verify.py','upper_range092/prepare.py',
      'upper_adaptive095/jobs.json.gz','upper_adaptive095/children_certificates.json.gz','upper_adaptive095/kernel.py',
      'upper_adaptive103/jobs.json.gz','upper_adaptive103/children_certificates.json.gz','upper_adaptive103/kernel.py',
      'upper_path112/nodes.json.gz','upper_path112/node_certificates.json.gz','upper_path112/numeric_hierarchy.json.gz',
      'upper_path112/geometries.json.gz','upper_path112/factor_catalog.json','upper_path112/production_v2.py',
      'audit_path121/check_whole_v2.py','audit_path121/whole_results.json',
      'upper_global080/engine_v2.py','upper_global080/initial_cells.json.gz',
      'upper_adapt082/kernel.py','upper_reenvelope088/engine.py']
    locks={str(BASE/n):sha(BASE/n) for n in names}
    for p in [P156/'plan.json',P156/'sources.json',P156/'common.py',A161/'binding.json']:
        locks[str(p)]=sha(p)
    rows=read(BASE/'upper_path112/nodes.json.gz');groups=defaultdict(list)
    for r in rows:groups[family(r['packet_path'])].append(r)
    selected={}
    for name,rr in sorted(groups.items()):
        for i in (0,len(rr)//2,len(rr)-1):selected[rr[i]['key']]=rr[i]
    special='5:31,31,31,31';selected[special]=next(r for r in rows if r['key']==special)
    assert len(selected)<=10
    write(HERE/'inventory_registration.json',dict(task='upper_integral167_readonly_metadata_schema',pid=os.getpid(),source_locks=locks,code_sha256=sha(__file__),
        fixed_sample_keys=sorted(selected),sample_packet_limit=10,bulk_packet_hash_or_read=False,bulk_export=False,
        new_H_R_J_factor_mean_range_chord_calls=0,old_successful_baseline_replays=0))
    sizes=Counter();counts=Counter();missing=[];allpaths=set();allshas=set()
    for r in rows:
        p=Path(r['packet_path']);g=family(p);counts[g]+=1;allpaths.add(str(p));allshas.add(r['packet_sha256'])
        try:sizes[g]+=p.stat().st_size
        except FileNotFoundError:missing.append(str(p))
    out=[]
    for k,r in sorted(selected.items()):
        p=Path(r['packet_path']);assert sha(p)==r['packet_sha256'];raw=read(p)
        assert raw['point']==r['point'] and len(raw['endpoints'])==73
        data,widths,exps,bits=payload(raw)
        out.append(dict(key=k,family=family(p),source_packet_sha256=r['packet_sha256'],raw_top_keys=sorted(raw),
          endpoint_keys=sorted(raw['endpoints'][0]),source_gzip_bytes=p.stat().st_size,prototype_binary_bytes=len(data),prototype_zlib_bytes=len(zlib.compress(data,9)),
          width_bits_minmax=[min(widths),max(widths)],exponent_minmax=[min(exps),max(exps)],lower_bits_minmax=[min(bits),max(bits)],
          top_source_sha_fields={n:v for n,v in raw.items() if 'sha256' in n},contains_weights='weights' in raw))
    rangecells=[]
    for p in sorted((BASE/'upper_range092').glob('block*_cells.json.gz')):rangecells.extend(read(p))
    original=read(BASE/'upper_range092/nodes.json.gz');by={r['key']:r for r in original}
    equal=Counter();keys=defaultdict(list)
    for r in rangecells:
        old=by[r['key']]
        for new,prior,label in [('rangeL_i','old_rangeL_i','rangeL_equals_old'),('rangeU_i','old_rangeU_i','rangeU_equals_old'),('mean_i','old_M_i','mean_equals_old'),('new_self_U','old_self_U','selfU_equals_old')]:
            if r[new]==old[prior]:equal[label]+=1;keys[label].append(r['key'])
    tree=read(BASE/'upper_path112/numeric_hierarchy.json.gz');geo=read(BASE/'upper_path112/geometries.json.gz');plan=read(P156/'plan.json')
    result=dict(status='metadata_only_inventory_complete',pid=os.getpid(),nodes=len(rows),unique_packet_paths=len(allpaths),unique_packet_recorded_shas=len(allshas),
      endpoint_intervals=len(rows)*73,endpoint_rational_fields=len(rows)*73*2,packet_family_counts=dict(counts),source_gzip_bytes_by_family=dict(sizes),source_gzip_bytes_total=sum(sizes.values()),missing_recorded_packets=missing,
      physical_frontier=len(tree['frontier']),roots=len(tree['roots']),internal_nodes=sum(bool(n['children']) for n in tree['nodes'].values()),
      geometry_counts=dict(Counter(r['geometry'] for r in rows)),levels=dict(Counter(r['level'] for r in rows)),
      metadata_file_sizes={str(BASE/n):(BASE/n).stat().st_size for n in names if n.endswith('.gz') or n.endswith('.json')},
      range092_nodes=len(rangecells),range092_stored_equalities=dict(equal),range092_equality_keys={k:v[:12] for k,v in keys.items()},
      samples=out,source156_new_mean_nodes=plan['new_mean_nodes'],source156_reused150_nodes=plan['reused150_unexpanded_roots'],
      source156_unchanged_nodes=plan['unchanged_self_nodes'],source156_unchanged_key=plan['unchanged_self_key'],
      original112_mean_unchanged_keys=[r['key'] for r in read(BASE/'upper_path112/node_certificates.json.gz') if r['new_M_i']==r['old_M_i']],
      all_source_locks_unchanged=all(sha(p)==h for p,h in locks.items()),sample_count=len(out),
      numeric_validity_of_sample_intervals_recomputed=False,new_H_R_J_factor_mean_range_chord_calls=0,old_successful_baseline_replays=0)
    write(HERE/'inventory.json',result)
    print(json.dumps({k:result[k] for k in ('status','pid','nodes','endpoint_intervals','packet_family_counts','source_gzip_bytes_total','roots','physical_frontier','internal_nodes','geometry_counts','range092_stored_equalities','sample_count','all_source_locks_unchanged')}))
if __name__=='__main__':main()
