"""Compatible least-vertex pulling of every proper face, exact boundary table."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import prod
import hashlib,time,os,json,importlib.util
from exact import save,read,sha,dumps,mask_of,ids,rank,determinant
P=Path(__file__).resolve().parent;CACHE=Path('/tmp/astra_vertex017')
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/numerics/body_fan_geometry.py')

def main():
    began=time.monotonic();contactpath=CACHE/'contacts.json.gz';facepath=CACHE/'faces.json.gz'
    data=read(contactpath);lattice=read(facepath);assert lattice['contacts_sha256']==sha(contactpath)
    source={str(p):sha(p) for p in [Path(__file__),P/'exact.py',contactpath,facepath,OLD]}
    fingerprint=hashlib.sha256(dumps(source)).hexdigest();cp=CACHE/'pulling_checkpoint.json.gz'
    Z=data['vertices'];masks=[int(x,16) for x in lattice['masks']];dims=lattice['dimensions'];children=lattice['children']
    order=sorted((i for i,d in enumerate(dims) if 0<=d<=4),key=lambda i:(dims[i],masks[i]));by_mask={m:i for i,m in enumerate(masks)}
    state=read(cp) if cp.exists() else {'fingerprint':fingerprint,'next':0,'triangulations':{},'restriction_checks':0}
    assert state['fingerprint']==fingerprint
    T={int(i):[tuple(s) for s in tris] for i,tris in state['triangulations'].items()}
    for pos in range(state['next'],len(order)):
        i=order[pos];Fmask=masks[i];d=dims[i];v=next(ids(Fmask))
        if d==0:tri=[(v,)]
        else:
            tri=sorted({tuple(sorted((v,)+s)) for j in children[i] if not(masks[j]>>v)&1 for s in T[j]})
            assert tri and all(len(s)==d+1 and len(set(s))==d+1 for s in tri)
        for s in tri:
            assert mask_of(s)&Fmask==mask_of(s)
            assert rank([Z[j] for j in s])==d+1
        # Check restriction to every codimension-one face, including all
        # lower-dimensional simplex intersections, not only top simplices.
        if d:
            for child in children[i]:
                intersections=[tuple(j for j in s if(masks[child]>>j)&1) for s in tri]
                assert {s for s in intersections if len(s)==d}==set(T[child])
                child_masks=[mask_of(s) for s in T[child]]
                assert all(any(mask_of(s)&g==mask_of(s) for g in child_masks) for s in intersections)
                state['restriction_checks']+=1
        T[i]=tri;state['next']=pos+1
        if (pos+1)%512==0 or time.monotonic()-began>120:
            state['triangulations']={str(k):v for k,v in T.items()};save(cp,state)
        if time.monotonic()-began>120:return
    state['triangulations']={str(k):v for k,v in T.items()};state['phase']='complete';save(cp,state)
    simplices=[];facets=[];seen=set();oriented_ridges=Counter();ridge_count=Counter();volume_num=0
    for facet in data['facets']:
        i=by_mask[mask_of(facet['vertices'])];assert dims[i]==4
        own=[]
        for s in T[i]:
            assert s not in seen;seen.add(s);D=determinant([Z[j] for j in s]);assert D
            for j in s:assert sum(a*b for a,b in zip(facet['normal'],Z[j]))==144
            own.append(len(simplices));simplices.append({'vertices':s,'det_numerator':abs(D)})
            volume_num+=abs(D)
            for j in range(5):
                ridge=s[:j]+s[j+1:];ridge_count[ridge]+=1;oriented_ridges[ridge]+=(-1)**j*(1 if D>0 else -1)
        facets.append({'normal':facet['normal'],'vertices':facet['vertices'],'simplices':own})
    assert all(n==2 for n in ridge_count.values()) and all(n==0 for n in oriented_ridges.values())
    spec=importlib.util.spec_from_file_location('locked017oldgeometry',OLD);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    old_sum=3840*sum((F(D,prod(q)) for _,D,q in old.base_terms()),F(0))
    new_sum=F(volume_num,144**5);assert old_sum==new_sum
    used=set(j for s in simplices for j in s['vertices']);assert used==set(range(280))
    table={'denominator':144,'vertices':Z,'facets':facets,'simplices':simplices,'geometry_source_hashes':source,
           'all_faces_input_sha256':sha(facepath),'pulling_order':'lexicographic 144*r, least vertex recursively'}
    out=P/'boundary_triangulation.json.gz';save(out,table)
    summary={'status':'complete_exact_original_vertex_boundary_pulling','vertices':280,'facets':len(facets),'boundary_4_simplices':len(simplices),
             'simplices_per_facet_histogram':dict(sorted(Counter(len(f['simplices']) for f in facets).items())),
             'boundary_3_simplices':len(ridge_count),'all_boundary_ridges_incidence_two':True,'all_oriented_ridge_sums_zero':True,
             'immediate_face_restriction_checks':state['restriction_checks'],'normalized_determinant_sum':str(new_sum),
             'volume_P':str(new_sum/120),'exact_volume_matches_old_26880_cones':True,'artificial_vertices_used':0,
             'table_sha256':sha(out),'table_bytes':out.stat().st_size,'source_hashes':source,
             'pid':os.getpid(),'seconds_this_process':time.monotonic()-began}
    save(P/'pulling_result.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
