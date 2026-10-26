#!/usr/bin/env python3
"""Verify the saved finite certificates; no solving, vertex search or kernel."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import gzip,hashlib,itertools,json,math,os,time
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from flint import fmpz_mat,fmpq_mat,fmpq
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
def read(p):
    p=Path(p)
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def detcols(xs):return F(str(fmpq_mat([[fmpq(str(x))for x in row]for row in zip(*xs)]).det()))
def write(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')

def verify_parent(fan,ledger,cols,priority,deadline,remaining_active):
    assert time.time()<deadline
    pid=fan['parent_id'];H=fan['boundary_forms_in_u'];m=len(H)-6
    rays=fan['parent_rays'];pd=abs(int(fmpz_mat(list(zip(*rays))).det()))
    assert F(fan['parent_determinant'])==pd
    vv=fan['vertices'];U=[tuple(map(F,v['u']))for v in vv];N=len(U)
    assert [v['id']for v in vv]==list(range(N))
    assert all(sum(u)==1 and min(u)>=0 for u in U)
    rr=[tuple(v['primitive_physical_ray'])for v in vv]
    orders=[(0,priority[r])if r in priority else(1,r)for r in rr]
    assert orders==sorted(orders) and len(set(rr))==N
    nums=[];dens=[]
    for i,u in enumerate(U):
        d=math.lcm(*(x.denominator for x in u));n=tuple(int(x*d)for x in u);nums.append(n);dens.append(d)
        scale=F(vv[i]['slice_scale']);assert scale>0 and math.gcd(*rr[i])==1
        assert all(rr[i][j]>=rr[i][j+1]for j in range(5))and rr[i][-1]>=0
        assert all(sum(F(rays[k][j])*u[k]for k in range(6))==scale*rr[i][j]for j in range(6))
    enum={v['enumeration_id']:i for i,v in enumerate(vv)};assert set(enum)==set(range(N))
    # Read the entire saved ordered classification ledger. No solve or new point.
    expected=iter(itertools.combinations(range(len(H)),5));counts={'S':0,'O':0,'I':0};seen_enum=set();count=0;digest=hashlib.sha256()
    with gzip.open(ledger,'rt')as f:
        for line in f:
            digest.update(line.encode());a,kind,*payload=json.loads(line);assert tuple(a)==next(expected)
            count+=1;assert count<=remaining_active
            if count%1024==0:assert time.time()<deadline
            M=fmpz_mat([[1]*6]+[H[j]for j in a]);d=int(M.det());counts[kind]+=1
            if kind=='S':assert d==0 and not payload
            elif kind=='O':
                assert d!=0 and len(payload)==1;y=list(map(F,payload[0]));assert sum(y)==1 and min(y)<0
                yd=math.lcm(*(x.denominator for x in y));yn=[int(x*yd)for x in y]
                assert all(dot(H[j],yn)==0 for j in a)
            else:
                assert kind=='I' and d!=0 and len(payload)==1;eid=payload[0];i=enum[eid]
                assert all(dot(H[j],nums[i])==0 for j in a)
                if eid not in seen_enum:
                    assert eid==len(seen_enum)and a==vv[i]['active_witness'];seen_enum.add(eid)
    assert next(expected,None)is None and count==math.comb(len(H),5)and seen_enum==set(enum)
    state=fan['active_enumeration']
    assert state['active_completed']==state['expected_active_systems']==count
    assert counts=={'S':state['singular'],'O':state['outside'],'I':state['inside_systems']}
    assert digest.hexdigest()==state['active_ledger_uncompressed_sha256']
    assert sha(ledger)==state['active_ledger_sha256']
    zeros=[frozenset(i for i,n in enumerate(nums)if dot(h,n)==0)for h in H]
    @lru_cache(None)
    def dim(s):return int(fmpz_mat([nums[i]for i in s]).rank())-1 if s else -1
    cells=fan['full_cells'];codes={c['sign_code']:c for c in cells};assert len(codes)==len(cells)
    def eligible(code):
        return tuple(i for i,n in enumerate(nums)if all((1 if(code>>j)&1 else -1)*dot(H[6+j],n)>=0 for j in range(m)))
    for c in cells:
        code=c['sign_code'];s=eligible(code);assert s==tuple(c['vertices'])and dim(s)==5
        assert c['signs']==[1 if(code>>j)&1 else -1 for j in range(m)]
        predicted=[]
        for j in range(m):
            t=tuple(i for i in s if i in zeros[6+j])
            if dim(t)==4:
                other=code^(1<<j);assert other in codes
                assert set(s)&set(codes[other]['vertices'])==set(t)
                predicted.append({'cut_index':j,'neighbor_sign_code':other,'facet_vertices':list(t)})
        assert c['adjacency']==predicted
    seed=fan['bfs_seed_sign_code'];assert seed in codes and fan['bfs_queue_exhausted']
    reached={seed};todo=[seed]
    while todo:
        c=codes[todo.pop()]
        for a in c['adjacency']:
            q=a['neighbor_sign_code']
            if q not in reached:reached.add(q);todo.append(q)
    assert reached==set(codes)
    assert set().union(*(set(c['vertices'])for c in cells))==set(range(N))
    # All facets via maximal proper zero faces, distinct from producer's rank filter.
    nodes={tuple(n['vertices']):n for n in fan['pulling_face_certificates']}
    assert len(nodes)==len(fan['pulling_face_certificates'])
    for s,n in nodes.items():
        assert time.time()<deadline
        d=dim(s);assert d==n['dimension']and n['apex']==min(s)
        if d==0:
            assert len(s)==1 and n['facets']==[]and n['simplices']==[list(s)];continue
        ss=set(s);proper={tuple(i for i in s if i in z)for z in zeros}
        proper.discard(());proper.discard(s)
        facets=sorted(t for t in proper if not any(set(t)<set(u)for u in proper))
        assert all(dim(t)==d-1 for t in facets)
        assert {tuple(t)for t in n['facets']}==set(facets)
        expected_simp=[]
        for t in n['facets']:
            if n['apex']not in t:
                child=nodes[tuple(t)];expected_simp.extend([[n['apex']]+q for q in child['simplices']])
        assert expected_simp==n['simplices']
        assert len({tuple(sorted(q))for q in expected_simp})==len(expected_simp)
    ds=[[dot(a,r)for a in cols]for r in rr]
    qs=[F(sum(abs(x)for x in row),80)for row in ds]
    assert [F(v['q_true_152'])for v in vv]==qs and min(qs)>0
    simplices=fan['simplices'];seen=set();total=F(0)
    for c in cells:
        s=c['vertices'];signs=[]
        for j in range(152):
            a=[ds[i][j]for i in s];assert not(min(a)<0<max(a));signs.append(1 if max(a)>0 else(-1 if min(a)<0 else 0))
        slope=[F(sum(signs[j]*cols[j][k]for j in range(152)),80)for k in range(6)]
        assert list(map(F,c['linear_symbol_coefficient']))==slope
        assert all(dot(slope,rr[i])==qs[i]for i in s)
        expected_simp=nodes[tuple(s)]['simplices']
        assert [simplices[i]['vertices']for i in c['simplex_indices']]==expected_simp
        cv=F(0)
        for sid in c['simplex_indices']:
            assert time.time()<deadline
            child=simplices[sid];sv=child['vertices'];assert child['id']==sid and child['cell']==c['id']
            key=tuple(sorted(sv));assert key not in seen;seen.add(key)
            volume=abs(detcols([U[i]for i in sv]));assert volume>0 and F(child['relative_slice_volume'])==volume
            determinant=int(fmpz_mat(list(zip(*(rr[i]for i in sv)))).det())
            assert F(child['primitive_physical_determinant'])==determinant
            assert abs(determinant)*math.prod(F(vv[i]['slice_scale'])for i in sv)==pd*volume
            for j in range(152):
                a=[ds[i][j]for i in sv];assert not(min(a)<0<max(a))
            cv+=volume
        assert F(c['relative_slice_volume'])==cv;total+=cv
    assert total==1 and len(seen)==len(simplices)
    return {'parent':pid,'active_classifications_verified_without_solves':count,'vertices':N,'full_cells':len(cells),
            'simplices':len(simplices),'pulling_faces':len(nodes),'relative_volume_exact':'1'}

def main():
    assert sys.argv[1:]in([],['--record'])
    start=time.monotonic();reg=read(HERE/'registration.json');run=read(HERE/'run_receipt.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    assert sha(HERE/'build_global.py')==reg['builder_sha256']and sha(HERE/'method.md')==reg['method_sha256']
    assert run['registration_sha256']==sha(HERE/'registration.json')
    assert run['all_handles_closed']and run['status']=='COMPLETE_21_PARENT_EXACT_GLOBAL_FAN'
    assert run['exit_code']==0 and run['completed_parents']==list(range(1,22))
    deadline=(run['actual_start_epoch'] if sys.argv[1:]==['--record'] else time.time())+reg['limits']['actual_seconds']
    assert time.time()<deadline
    frame=read(reg['frame_path']);cols=frame['column_numerators']
    inherited=read(reg['inherited_parent1_fan']);assert sha(reg['inherited_parent1_fan'])==reg['inherited_parent1_sha256']
    priority={tuple(v['primitive_physical_ray']):v['id']for v in inherited['vertices']}
    summaries=[];verified=0;fans={1:inherited}
    for saved in run['new_parent_receipts']:
        pid=saved['parent'];assert pid==len(summaries)+2
        assert sha(saved['fan_path'])==saved['fan_sha256'];fan=read(saved['fan_path']);fans[pid]=fan
        s=verify_parent(fan,HERE/f'parents/{pid:02d}.active.jsonl.gz',cols,priority,deadline,
                        reg['limits']['active_systems']-run['budget']['active_systems']-verified)
        summaries.append(s);verified+=s['active_classifications_verified_without_solves']
        print(json.dumps({'saved_parent_verified':pid,'active_classifications':verified,'seconds':time.monotonic()-start}),flush=True)
    # Independently reconstruct final aggregation and exact shared-face sets.
    global_fan=read(HERE/'global_fan.json.gz');cross=read(HERE/'cross_parent_faces.json.gz')
    assert sha(HERE/'global_fan.json.gz')==run['global_fan_sha256']
    assert sha(HERE/'cross_parent_faces.json.gz')==run['cross_parent_faces_sha256']
    global_rays=global_fan['rays'];rmap={tuple(x['primitive']):x['id']for x in global_rays}
    assert len(rmap)==len(global_rays)
    facets={};flat=[];total=F(0)
    for pid,fan in sorted(fans.items()):
        vv=fan['vertices'];ss=fan['simplices'];pr=[tuple(r)for r in fan['parent_rays']]
        for s in ss:
            rr=[tuple(vv[i]['primitive_physical_ray'])for i in s['vertices']]
            dd=F(s['primitive_physical_determinant']);assert dd.denominator==1
            dd=abs(int(dd));total+=F(dd,math.prod(r[0]for r in rr))
            flat.append({'parent':pid,'local_id':s['id'],'abs_primitive_determinant':dd,'ray_ids':[rmap[r]for r in rr]})
            for i,r in zip(s['vertices'],rr):assert F(global_rays[rmap[r]]['q_true_152'])==F(vv[i]['q_true_152'])
        for j in range(6):
            k=tuple(sorted(pr[:j]+pr[j+1:]));on={i for i,v in enumerate(vv)if F(v['u'][j])==0}
            tri=set()
            for s in ss:
                a=on.intersection(s['vertices'])
                if len(a)==5:tri.add(tuple(sorted(tuple(vv[i]['primitive_physical_ray'])for i in a)))
            assert tri;facets.setdefault(k,[]).append((pid,j,tri))
    assert total==1 and flat==global_fan['cones']
    walls=[tuple(int(i==j)-int(i==j+1)for i in range(6))for j in range(5)]+[(0,0,0,0,0,1)]
    expected_shared=[];expected_boundary=[]
    for k,owners in sorted(facets.items()):
        assert len(owners)in(1,2)
        if len(owners)==2:
            a,b=owners;assert a[2]==b[2]
            expected_shared.append({'parents':[a[0],b[0]],'parent_facets':[a[1],b[1]],'old_facet_rays':[list(x)for x in k],
                                    'triangulation':[[list(r)for r in q]for q in sorted(a[2])]})
        else:
            a=owners[0];which=[j for j,w in enumerate(walls)if all(dot(w,r)==0 for r in k)];assert which
            expected_boundary.append({'parent':a[0],'parent_facet':a[1],'chamber_walls':which,'old_facet_rays':[list(x)for x in k],
                                      'triangulation':[[list(r)for r in q]for q in sorted(a[2])]})
    assert cross=={'shared_facets':expected_shared,'boundary_facets':expected_boundary,'all_matched':True}
    assert global_fan['full_signed_permutation_multiplicity']==46080
    assert global_fan['logical_terms_per_future_full_kernel']==46080*len(flat)
    assert verified==run['budget']['active_systems']==reg['planned_total_new_active_systems']
    assert 2*verified<=reg['limits']['active_systems']and time.time()<deadline
    out={'status':'PASS_COMPLETE_SAVED_CERTIFICATE_NO_NEW_SOLVES','pid':os.getpid(),'wall_seconds':time.monotonic()-start,
         'total_actual_clock_seconds':time.time()-run['actual_start_epoch'],'new_active_solves':0,
         'saved_active_classifications_verified':verified,'construction_plus_verification_system_operations':2*verified,
         'first_parent_reenumerations':0,'parents':summaries,'global_rays':len(global_rays),'global_simplices':len(flat),
         'shared_old_facets':len(expected_shared),'boundary_old_facets':len(expected_boundary),'normalized_xi1_slice_volume':'1',
         'checker_sha256':sha(__file__),'kernel_G_H_K_p_native':0,'new_C6_bound':False}
    if sys.argv[1:]==['--record']:
        assert not(HERE/'checked_results.json').exists();write(HERE/'checked_results.json',out)
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
