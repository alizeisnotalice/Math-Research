#!/usr/bin/env python3
"""B050: exact finite arrangement construction, no kernel or time parameter."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import deque
from functools import lru_cache
import datetime,gzip,hashlib,json,math,os,time,traceback
import flint
from flint import fmpq_mat,fmpz_mat,fmpq
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def write(p,x):
    p=Path(p);tmp=p.with_name(p.name+'.tmp')
    if p.suffix=='.gz':
        with gzip.open(tmp,'wt',compresslevel=6)as f:json.dump(x,f,separators=(',',':'));f.write('\n')
    else:tmp.write_text(json.dumps(x,indent=2)+'\n')
    os.replace(tmp,p)
def qmat(rows):return fmpq_mat([[fmpq(str(x)) for x in row]for row in rows])
def detcols(cols):return Q(str(qmat(list(zip(*cols))).det()))
def primitive(v):
    d=math.lcm(*(x.denominator for x in v));n=[int(x*d)for x in v];g=math.gcd(*n)
    assert g>0
    return tuple(x//g for x in n),Q(g,d)
def bits(mask):
    while mask:
        b=mask&-mask;yield b.bit_length()-1;mask-=b
class StopBatch(RuntimeError):pass
class Budget:
    def __init__(self,reg):
        self.start=time.monotonic();self.limits=reg['limits'];self.active=0;self.cells=12;self.simplices=14
    def check(self):
        if time.monotonic()-self.start>=self.limits['actual_seconds']:raise StopBatch('actual-time cap reached')
        if self.active>self.limits['active_systems']:raise StopBatch('active-system cap reached')
        if self.cells>self.limits['full_cells']:raise StopBatch('full-cell cap reached')
        if self.simplices>self.limits['full_simplices']:raise StopBatch('full-simplex cap reached')
    def state(self):return {'active_systems':self.active,'full_cells_including_inherited':self.cells,'full_simplices_including_inherited':self.simplices,'actual_elapsed_seconds':time.monotonic()-self.start}

def sourcecheck(reg):
    for p,h in reg['source_sha256'].items():assert sha(p)==h,('source mismatch',p)
    assert sha(__file__)==reg['builder_sha256']
    assert sha(HERE/'method.md')==reg['method_sha256']

def parse_old(path):
    it=iter(map(int,Path(path).read_text().split()));n=next(it);assert n==21;out=[]
    for _ in range(n):
        d=next(it);rr=[]
        for _ in range(6):
            next(it);rr.append(tuple(next(it)for _ in range(6))) # old q deliberately discarded
        out.append((d,rr))
    assert next(it,None)is None
    return out

def build_parent(pid,parent,cols,priority,budget,checkpoint):
    start=time.monotonic();pd,rays=parent
    assert abs(int(fmpz_mat(list(zip(*rays))).det()))==pd
    parent_dots=[tuple(dot(w,r)for r in rays)for w in cols[32:]]
    cross=[j for j,v in enumerate(parent_dots)if min(v)<0<max(v)]
    cuts=[]
    for j in cross:
        g=math.gcd(*parent_dots[j]);cuts.append(tuple(a//g for a in parent_dots[j]))
    assert len(set(cuts))==len(cuts)
    assert not any(tuple(-x for x in h)in cuts for h in cuts)
    H=[tuple(int(i==j)for i in range(6))for j in range(6)]+cuts
    m=len(cuts);expected=math.comb(len(H),5)
    state={'parent':pid,'phase':'ACTIVE_SYSTEMS','expected_active_systems':expected,'active_completed':0,
           'singular':0,'outside':0,'inside_systems':0,'unique_vertices':0}
    checkpoint['current_parent']=state;write(HERE/'checkpoint.json',checkpoint)
    ledger_path=HERE/f'parents/{pid:02d}.active.jsonl.gz';digest=hashlib.sha256();vertices={}
    ones=[1]*6;rhs=fmpq_mat([[1],[0],[0],[0],[0],[0]])
    last_feedback=time.monotonic()
    try:
        with gzip.open(ledger_path,'wt',compresslevel=6)as ledger:
            for ai,active in enumerate(combinations(range(len(H)),5),1):
                budget.check();budget.active+=1;budget.check()
                M=fmpq_mat([ones]+[H[j]for j in active])
                try:
                    sol=M.solve(rhs)
                except ZeroDivisionError:
                    state['singular']+=1;rec=[list(active),'S']
                else:
                    vals=[sol[i,0]for i in range(6)]
                    if any(v<0 for v in vals):
                        state['outside']+=1;rec=[list(active),'O',[str(v)for v in vals]]
                    else:
                        u=tuple(Q(str(v))for v in vals)
                        assert sum(u)==1
                        assert all(dot(H[j],u)==0 for j in active)
                        if u not in vertices:vertices[u]=(len(vertices),tuple(active))
                        state['inside_systems']+=1;rec=[list(active),'I',vertices[u][0]]
                line=json.dumps(rec,separators=(',',':'))+'\n';ledger.write(line);digest.update(line.encode())
                state['active_completed']=ai;state['unique_vertices']=len(vertices)
                if time.monotonic()-last_feedback>25:
                    checkpoint['budget']=budget.state();write(HERE/'checkpoint.json',checkpoint)
                    print(json.dumps({'progress':state,'budget':budget.state()}),flush=True);last_feedback=time.monotonic()
        assert state['active_completed']==expected
        assert state['singular']+state['outside']+state['inside_systems']==expected
        state['active_ledger_uncompressed_sha256']=digest.hexdigest()
        # Convert to globally ordered primitive physical rays before any pulling.
        processed=[]
        for u,(eid,active)in vertices.items():
            x=tuple(sum(Q(rays[j][i])*u[j]for j in range(6))for i in range(6))
            r,scale=primitive(x)
            assert all(r[j]>=r[j+1]for j in range(5))and r[-1]>=0
            order=(0,priority[r])if r in priority else(1,r)
            processed.append((order,u,r,scale,eid,active))
        processed.sort(key=lambda x:x[0]);U=[x[1]for x in processed];N=len(U)
        primitive_rays=[x[2]for x in processed];scales=[x[3]for x in processed]
        assert len(set(primitive_rays))==N
        nums=[]
        for u in U:
            den=math.lcm(*(x.denominator for x in u));nums.append(tuple(int(x*den)for x in u))
        zmask=[];positive=[];negative=[];allmask=(1<<N)-1
        for h in H:
            vv=[dot(h,u)for u in nums]
            zmask.append(sum(1<<i for i,v in enumerate(vv)if not v))
            positive.append(sum(1<<i for i,v in enumerate(vv)if v>=0))
            negative.append(sum(1<<i for i,v in enumerate(vv)if v<=0))
        @lru_cache(None)
        def dimension(mask):
            return int(fmpz_mat([nums[i]for i in bits(mask)]).rank())-1 if mask else -1
        @lru_cache(None)
        def cell_mask(code):
            s=allmask
            for j in range(m):s&=positive[6+j]if(code>>j)&1 else negative[6+j]
            return s
        for z in range(1,5*m+2):
            generic=[z**j for j in range(6)];v=[dot(h,generic)for h in cuts]
            if all(v):break
        else:raise AssertionError('finite generic seed bound failed')
        seed=sum((a>0)<<j for j,a in enumerate(v))
        assert dimension(cell_mask(seed))==5
        queue=deque([seed]);discovered={seed};cells=[];code_to_id={}
        state['phase']='BFS_FULL_CELLS';state['generic_seed_positive_integers']=generic
        while queue:
            budget.check();code=queue.popleft();s=cell_mask(code);assert dimension(s)==5
            cid=len(cells);code_to_id[code]=cid;budget.cells+=1;budget.check()
            adjacent=[]
            for j in range(m):
                t=s&zmask[6+j]
                if dimension(t)==4:
                    neighbor=code^(1<<j);ns=cell_mask(neighbor)
                    assert dimension(ns)==5 and (s&ns)==t
                    adjacent.append({'cut_index':j,'neighbor_sign_code':neighbor,'facet_vertices':list(bits(t))})
                    if neighbor not in discovered:discovered.add(neighbor);queue.append(neighbor)
            cells.append({'id':cid,'sign_code':code,'signs':[1 if(code>>j)&1 else -1 for j in range(m)],
                          'vertices':list(bits(s)),'affine_dimension':5,'adjacency':adjacent})
            state['full_cells_discovered']=len(discovered);state['full_cells_settled']=len(cells)
        assert len(cells)==len(discovered)==len(code_to_id)
        assert all(a['neighbor_sign_code']in code_to_id for c in cells for a in c['adjacency'])
        covered_vertices=0
        for c in cells:covered_vertices|=cell_mask(c['sign_code'])
        assert covered_vertices==allmask
        # Complete exact facets and global-order pulling, with a saved face DAG.
        faces={};pull_cache={}
        def facets(mask):
            d=dimension(mask)
            return sorted({mask&z for z in zmask if mask&z and mask&z!=mask and dimension(mask&z)==d-1})
        def pulling(mask):
            budget.check()
            if mask in pull_cache:return pull_cache[mask]
            d=dimension(mask);assert d>=0;apex=next(bits(mask));fs=[]
            if d==0:
                assert mask.bit_count()==1;out=[(apex,)]
            else:
                fs=facets(mask);assert fs;opposite=[t for t in fs if not(t>>apex)&1];assert opposite
                out=[]
                for t in opposite:
                    for q in pulling(t):
                        simplex=(apex,)+q;assert len(simplex)==d+1
                        out.append(simplex)
                        if d==5 and len(out)+budget.simplices>budget.limits['full_simplices']:
                            raise StopBatch('simplex cap during pulling')
                assert len(set(tuple(sorted(q))for q in out))==len(out)
            faces[mask]={'vertices':list(bits(mask)),'dimension':d,'apex':apex,
                         'facets':[list(bits(t))for t in fs],'simplices':[list(q)for q in out]}
            pull_cache[mask]=out;return out
        ds=[tuple(dot(a,r)for a in cols)for r in primitive_rays]
        rates=[Q(sum(abs(v)for v in a),80)for a in ds];assert min(rates)>0
        vr=[]
        for i,(_,u,r,scale,eid,active)in enumerate(processed):
            vr.append({'id':i,'u':list(map(str,u)),'enumeration_id':eid,'active_witness':list(active),
                       'primitive_physical_ray':list(r),'slice_scale':str(scale),'q_true_152':str(rates[i])})
        simplices=[];seen_simp=set();volume=Q(0);state['phase']='PULLING_AND_TRUE_FRAME_VERIFICATION'
        for c in cells:
            budget.check();s=cell_mask(c['sign_code']);vs=c['vertices']
            signs=[]
            for j in range(152):
                vals=[ds[i][j]for i in vs];assert not(min(vals)<0<max(vals))
                signs.append(1 if max(vals)>0 else(-1 if min(vals)<0 else 0))
            slope=[Q(sum(signs[j]*cols[j][a]for j in range(152)),80)for a in range(6)]
            assert all(dot(slope,primitive_rays[i])==rates[i]for i in vs)
            c['linear_symbol_coefficient']=list(map(str,slope));c['simplex_indices']=[];cv=Q(0)
            for simp in pulling(s):
                budget.check();sk=tuple(sorted(simp));assert sk not in seen_simp;seen_simp.add(sk)
                vol=abs(detcols([U[i]for i in simp]));assert vol>0
                physical=int(fmpz_mat(list(zip(*(primitive_rays[i]for i in simp)))).det());assert physical
                assert abs(physical)*math.prod(scales[i]for i in simp)==pd*vol
                for j in range(152):
                    vv=[ds[i][j]for i in simp];assert not(min(vv)<0<max(vv))
                sid=len(simplices);c['simplex_indices'].append(sid)
                simplices.append({'id':sid,'cell':c['id'],'vertices':list(simp),
                                  'relative_slice_volume':str(vol),'primitive_physical_determinant':str(physical)})
                budget.simplices+=1;budget.check();cv+=vol
            c['relative_slice_volume']=str(cv);volume+=cv
        assert volume==1
        fan={'scope':f'COMPLETE_OLD_PARENT_{pid}','parent_id':pid,'parent_rays':[list(r)for r in rays],
             'parent_determinant':str(pd),'normalization':'sum u=1; relative parent volume=1',
             'crossing_new_form_indices_1_based':[j+1 for j in cross],
             'boundary_forms_in_u':[list(h)for h in H],'vertices':vr,'full_cells':cells,
             'pulling_face_certificates':list(faces.values()),'simplices':simplices,
             'active_enumeration':state,'bfs_seed_sign_code':seed,'bfs_queue_exhausted':True,
             'complete_vertex_pool_used':True,'lower_dimensional_patterns':'retained by closure of all full cells; no 2^m enumeration',
             'global_pulling_order':'inherited parent1 primitive-ray IDs first, then remaining primitive rays lexicographic',
             'relative_volume_exact':'1','frame_sha256':checkpoint['frame_sha256']}
        state.update(phase='COMPLETE',vertices=N,full_cells=len(cells),simplices=len(simplices),
                     face_nodes=len(faces),relative_volume_exact='1',elapsed_seconds=time.monotonic()-start,
                     active_ledger_sha256=sha(ledger_path))
        fanpath=HERE/f'parents/{pid:02d}.fan.json.gz';write(fanpath,fan)
        return fan,{'parent':pid,'status':'COMPLETE','fan_path':str(fanpath),'fan_sha256':sha(fanpath),**state}
    except BaseException:
        partial={'state':state,'vertices':[{'u':list(map(str,u)),'enumeration_id':eid,'active_witness':list(a)}for u,(eid,a)in vertices.items()],
                 'uncovered':True,'no_complete_parent_claim':True}
        write(HERE/f'parents/{pid:02d}.partial.json.gz',partial)
        checkpoint['current_parent']=state;checkpoint['budget']=budget.state();write(HERE/'checkpoint.json',checkpoint)
        raise

def global_checks(fans,cols,reg,budget):
    facets={};rayrates={};cones=[];global_volume=Q(0)
    for pid,fan in sorted(fans.items()):
        budget.check();vv=fan['vertices'];simplices=fan['simplices']
        for v in vv:
            r=tuple(v['primitive_physical_ray']);q=Q(v['q_true_152'])
            assert r not in rayrates or rayrates[r]==q;rayrates[r]=q
        for s in simplices:
            rr=[tuple(vv[i]['primitive_physical_ray'])for i in s['vertices']]
            dq=Q(s['primitive_physical_determinant']);assert dq.denominator==1
            d=abs(int(dq));assert d>0
            global_volume+=Q(d,math.prod(r[0]for r in rr))
            cones.append({'parent':pid,'local_id':s['id'],'rays':rr,'abs_primitive_determinant':d})
        pr=[tuple(r)for r in fan['parent_rays']]
        for j in range(6):
            key=tuple(sorted(pr[:j]+pr[j+1:]));on={i for i,v in enumerate(vv)if Q(v['u'][j])==0}
            induced=set()
            for s in simplices:
                face=on.intersection(s['vertices'])
                if len(face)==5:induced.add(tuple(sorted(tuple(vv[i]['primitive_physical_ray'])for i in face)))
            assert induced
            facets.setdefault(key,[]).append({'parent':pid,'parent_facet':j,'triangulation':induced})
    assert global_volume==1
    walls=[tuple(int(i==j)-int(i==j+1)for i in range(6))for j in range(5)]+[(0,0,0,0,0,1)]
    shared=[];boundary=[]
    for key,owners in sorted(facets.items()):
        budget.check();assert len(owners)in(1,2)
        if len(owners)==2:
            a,b=owners;assert a['triangulation']==b['triangulation'],('cross-parent facet mismatch',a['parent'],b['parent'])
            shared.append({'parents':[a['parent'],b['parent']],'parent_facets':[a['parent_facet'],b['parent_facet']],
                           'old_facet_rays':key,'triangulation':[list(q)for q in sorted(a['triangulation'])]})
        else:
            a=owners[0];which=[j for j,w in enumerate(walls)if all(dot(w,r)==0 for r in key)];assert which
            boundary.append({'parent':a['parent'],'parent_facet':a['parent_facet'],'chamber_walls':which,
                             'old_facet_rays':key,'triangulation':[list(q)for q in sorted(a['triangulation'])]})
    first=[tuple(v['primitive_physical_ray'])for v in fans[1]['vertices']];priority={r:i for i,r in enumerate(first)}
    rr=sorted(rayrates,key=lambda r:(0,priority[r])if r in priority else(1,r));rid={r:i for i,r in enumerate(rr)}
    for c in cones:c['ray_ids']=[rid[tuple(r)]for r in c.pop('rays')]
    data={'status':'COMPLETE_21_PARENT_EXACT_SPARSE4_FAN','beta':'1/80','frame_sha256':reg['frame_sha256'],
          'rays':[{'id':i,'primitive':r,'q_true_152':str(rayrates[r])}for i,r in enumerate(rr)],
          'cones':cones,'parent_count':21,'ray_count':len(rr),'cone_count':len(cones),
          'full_signed_permutation_multiplicity':46080,'logical_terms_per_future_full_kernel':46080*len(cones),
          'ordered_xi1_slice_volume':'1/120','ordered_xi1_slice_normalized_volume':'1',
          'global_negative_pair_factor_not_applied':True,'kernel_calls':0,
          'numeric_engine_int64_or_zero_node_compatibility_claimed':False}
    write(HERE/'global_fan.json.gz',data)
    write(HERE/'cross_parent_faces.json.gz',{'shared_facets':shared,'boundary_facets':boundary,'all_matched':True})
    return {'global_fan_sha256':sha(HERE/'global_fan.json.gz'),'cross_parent_faces_sha256':sha(HERE/'cross_parent_faces.json.gz'),
            'global_rays':len(rr),'global_simplices':len(cones),'shared_old_facets':len(shared),'boundary_old_facets':len(boundary),
            'logical_full_group_terms':46080*len(cones),'ordered_slice_volume_exact':'1/120'}

def main():
    assert sys.argv[1:]==['--run']
    assert not(HERE/'run_receipt.json').exists(),'single registered batch; no checkpoint reset or rerun'
    reg=read(HERE/'registration.json');assert reg['release_status']=='EXPLICIT_MAIN_RELEASE_AFTER_F046_FINAL_PASS'
    sourcecheck(reg)
    parent1=read(reg['inherited_parent1_fan']);assert len(parent1['simplices'])==14
    assert sha(reg['inherited_parent1_fan'])==reg['inherited_parent1_sha256']
    frame=read(reg['frame_path']);cols=[tuple(v)for v in frame['column_numerators']]
    assert len(cols)==152 and frame['common_denominator']==80 and frame['beta']=='1/80'
    old=parse_old(reg['old_fan_path']);priority={tuple(v['primitive_physical_ray']):v['id']for v in parent1['vertices']}
    assert len(priority)==12
    (HERE/'parents').mkdir(exist_ok=False);budget=Budget(reg)
    receipt={'status':'RUNNING','pid':os.getpid(),'actual_start_epoch':time.time(),'actual_start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'registration_sha256':sha(HERE/'registration.json'),'completed_parents':[1],'new_parent_receipts':[],
             'inherited_parent1_no_reenumeration':True,'new_kernel_G_H_K_p_native':0}
    checkpoint={'status':'PARTIAL','frame_sha256':reg['frame_sha256'],'completed_parents':[1],
                'uncovered_parents':list(range(2,22)),'budget':budget.state(),
                'unverified_global_obligations':['cross-parent compatibility','global aggregation']}
    write(HERE/'run_receipt.json',receipt);write(HERE/'checkpoint.json',checkpoint);fans={1:parent1};code=1
    try:
        for pid in range(2,22):
            budget.check();sourcecheck(reg)
            fan,pstate=build_parent(pid,old[pid-1],cols,priority,budget,checkpoint);fans[pid]=fan
            receipt['new_parent_receipts'].append(pstate);receipt['completed_parents'].append(pid)
            checkpoint.update(completed_parents=receipt['completed_parents'][:],uncovered_parents=list(range(pid+1,22)),budget=budget.state())
            checkpoint.pop('current_parent',None);write(HERE/'checkpoint.json',checkpoint);write(HERE/'run_receipt.json',receipt)
            print(json.dumps({'parent_completed':pid,'active_systems':pstate['active_completed'],'vertices':pstate['vertices'],
                              'cells':pstate['full_cells'],'simplices':pstate['simplices'],'seconds':pstate['elapsed_seconds'],
                              'budget':budget.state()}),flush=True)
        g=global_checks(fans,cols,reg,budget);sourcecheck(reg);budget.check()
        receipt.update(status='COMPLETE_21_PARENT_EXACT_GLOBAL_FAN',**g);code=0
        checkpoint.update(status='COMPLETE',uncovered_parents=[],global_checks=g,unverified_global_obligations=[])
    except StopBatch as exc:
        receipt.update(status='PARTIAL_REGISTERED_CAP_STOP',stop_reason=str(exc));code=0
        checkpoint['status']='PARTIAL_REGISTERED_CAP_STOP'
    except BaseException as exc:
        receipt.update(status='PARTIAL_FAILURE_NO_RETRY',error=repr(exc),traceback=traceback.format_exc());code=1
        checkpoint['status']='PARTIAL_FAILURE_NO_RETRY'
    finally:
        receipt.update(exit_code=code,budget=budget.state(),all_handles_closed=True)
        checkpoint['budget']=budget.state();write(HERE/'checkpoint.json',checkpoint);write(HERE/'run_receipt.json',receipt)
        print(json.dumps({'final_status':receipt['status'],'exit_code':code,'budget':budget.state(),'completed':receipt['completed_parents'],'uncovered':checkpoint['uncovered_parents']}),flush=True)
    return code

if __name__=='__main__':raise SystemExit(main())
