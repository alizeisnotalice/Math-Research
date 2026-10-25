"""Complete 65536-sign facet contacts against the locked280 true vertices."""
from pathlib import Path
from itertools import permutations,product
from collections import Counter
import hashlib,json,os,time,sys
import numpy as np
from exact import rank,mask_of,save,read,sha,dumps

P=Path(__file__).resolve().parent;CACHE=Path('/tmp/astra_vertex017')
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-audit/cube-dim-5-bounds/audit')

def main():
    began=time.monotonic();source={x.name:sha(x) for x in [Path(__file__),P/'exact.py',P/'registration.md',OLD/'005_vertex_rays_exact.py',OLD/'005_vertex_rays_exact_results.json']}
    fingerprint=hashlib.sha256(dumps(source)).hexdigest()
    locked=read(OLD/'005_vertex_rays_exact_results.json');assert locked['source_sha256']==sha(OLD/'005_vertex_rays_exact.py')
    orbits=[(tuple(r['primitive_ray']),r['q_equals_psi_ray']) for r in locked['records'] if r['normalized_ray_is_true_vertex']]
    assert orbits==[((1,1,0,0,0),16),((1,1,1,1,0),24),((2,1,1,1,1),36)]
    vertex_orbit={}
    for orbit,(n,q) in enumerate(orbits):
        for p in set(permutations(n)):
            for s in product((-1,1),repeat=5):
                z=tuple(144//q*a*b for a,b in zip(p,s));assert z not in vertex_orbit or vertex_orbit[z]==orbit
                vertex_orbit[z]=orbit
    vertices=sorted(vertex_orbit);assert len(vertices)==280
    counts=Counter(vertex_orbit.values());assert [counts[i] for i in range(3)]==[40,80,160]
    V=np.array([(1,)+s for s in product((-1,1),repeat=4)],dtype=np.int64)
    Z=np.array(vertices,dtype=np.int64);assert np.max(np.abs(Z))==9
    assert np.all(np.abs(V@Z.T).sum(axis=0)==144)
    # Each entry is an exact sum of16 signs; int64 range is ample.
    bits=(np.arange(65536,dtype=np.int64)[:,None]>>np.arange(16,dtype=np.int64))&1
    all_normals=(2*bits-1)@V
    normals,multiplicities=np.unique(all_normals,axis=0,return_counts=True)
    assert len(normals)==5817 and int(multiplicities.sum())==65536
    assert np.max(np.abs(normals))==16
    checkpoint=CACHE/'contacts_checkpoint.json.gz';records=[];start=0
    if checkpoint.exists():
        prior=read(checkpoint);assert prior['fingerprint']==fingerprint
        records=prior['records'];start=prior['next_slope']
    for begin in range(start,len(normals),128):
        end=min(begin+128,len(normals));values=normals[begin:end]@Z.T
        assert np.max(values)<=144 # all known vertices feasible for every slope
        for off,a in enumerate(normals[begin:end]):
            contact=np.flatnonzero(values[off]==144).tolist();r=rank([vertices[i] for i in contact])
            records.append({'normal':[int(x) for x in a],'multiplicity':int(multiplicities[begin+off]),'contacts':contact,'rank':r})
        if end%512==0 or end==len(normals) or time.monotonic()-began>120:
            save(checkpoint,{'fingerprint':fingerprint,'next_slope':end,'records':records})
            print(json.dumps({'pid':os.getpid(),'stage':'contacts','processed_slopes':end,'facets_so_far':sum(r['rank']==5 for r in records),'seconds':time.monotonic()-began}),flush=True)
        if time.monotonic()-began>120 and end<len(normals):return
    assert len(records)==5817
    facets=[{'normal':r['normal'],'vertices':r['contacts']} for r in records if r['rank']==5]
    assert len({tuple(f['vertices']) for f in facets})==len(facets)
    data={'denominator':144,'vertices':vertices,'orbit_labels':[vertex_orbit[z] for z in vertices],'facets':facets,
          'source_hashes':source,'fingerprint':fingerprint,'all_contact_records_sha256':hashlib.sha256(dumps(records)).hexdigest()}
    out=CACHE/'contacts.json.gz';save(out,data)
    summary={'status':'complete_exact_facet_contacts','sign_choices':65536,'distinct_slopes':5817,'vertices':280,
             'orbit_counts':[40,80,160],'facets':len(facets),'facet_vertex_counts':dict(sorted(Counter(len(f['vertices']) for f in facets).items())),
             'slope_contact_rank_counts':dict(sorted(Counter(r['rank'] for r in records).items())),
             'contacts_path':str(out),'contacts_sha256':sha(out),'source_hashes':source,
             'pid':os.getpid(),'seconds':time.monotonic()-began,'python':sys.executable,'numpy':np.__version__}
    save(P/'contacts_result.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
