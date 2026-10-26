#!/usr/bin/env python3
"""Independent finite algebra readback. No producer imports and no union calls."""
import time
START=time.monotonic()
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product,combinations
from math import prod
import json,gzip,hashlib,os,sys,resource

HERE=Path(__file__).resolve().parent
ENV=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
ZERO=()
MB={}


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def load(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else Path(p).open()) as f:return json.load(f)


def save(name,o):
    (HERE/name).write_text(json.dumps(o,indent=2)+'\n')


@lru_cache(maxsize=4096)
def rational(s):return Q(s)


def canonical(d):return tuple(sorted((k,v) for k,v in d.items() if v))


def dec(o):
    d={}
    for k,s in o.items():
        i,j=map(int,k.split(','));assert 0<=i<5 and 0<=j<5
        d[i,j]=rational(s)
    return canonical(d)


def enc(p):return {f'{i},{j}':str(c) for (i,j),c in p}


def mono(i,j,c=1):
    return (((i%5,j%5),Q(c)*2**(i//5)*3**(j//5)),) if c else ZERO


ONE=mono(0,0)
U=mono(1,0)
V=mono(0,1)


def scale(p,c):return tuple((ij,x*c) for ij,x in p if x*c)


def add(*ps):
    d={}
    for p in ps:
        for ij,c in p:d[ij]=d.get(ij,Q(0))+c
    return canonical(d)


def sub(p,q):return add(p,scale(q,-1))


@lru_cache(maxsize=8192)
def mul(p,q):
    d={}
    for (i,j),x in p:
        for (k,l),y in q:
            key=((i+k)%5,(j+l)%5)
            d[key]=d.get(key,Q(0))+x*y*2**((i+k)//5)*3**((j+l)//5)
    return canonical(d)


def power(p,k):
    ans=ONE
    for _ in range(k):ans=mul(ans,p)
    return ans


def root_box(n):
    s=1<<256;lo=s;hi=2*s;target=n*s**5
    for _ in range(256):
        mid=(lo+hi)//2
        if mid**5<=target:lo=mid
        else:hi=mid
    assert hi==lo+1 and lo**5<target<hi**5
    return Q(lo,s),Q(hi,s)


@lru_cache(maxsize=4096)
def bounds(p):
    lo=Q(0);hi=Q(0)
    for ij,c in p:
        a,b=MB[ij]
        lo+=c*(a if c>0 else b);hi+=c*(b if c>0 else a)
    return lo,hi


def cmp(p,q=ZERO):
    d=sub(p,q)
    if not d:return 0
    lo,hi=bounds(d)
    assert lo>0 or hi<0,('unresolved sign',enc(d))
    return 1 if lo>0 else -1


def pmin(p,q):return p if cmp(p,q)<=0 else q
def pmax(p,q):return p if cmp(p,q)>=0 else q


def inside_saved(p,interval):
    lo,hi=bounds(p);assert Q(interval[0])<=lo<=hi<=Q(interval[1])


def radius(n):
    i=j=0;m=n
    while m%2==0:m//=2;i+=1
    while m%3==0:m//=3;j+=1
    assert m==1
    return mono(i,j,Q(1,2))


def stream_histogram(path):
    """Own bounded-buffer JSON tokenizer, no producer streaming helper."""
    decoder=json.JSONDecoder();found=False
    with gzip.open(path,'rt') as f:
        buf='';pos=0;eof=False
        def fill():
            nonlocal buf,pos,eof
            buf=buf[pos:];pos=0
            chunk=f.read(1<<18)
            if not chunk:eof=True
            buf+=chunk
        def whitespace():
            nonlocal pos
            while True:
                while pos<len(buf) and buf[pos].isspace():pos+=1
                if pos<len(buf) or eof:return
                fill()
        def take():
            nonlocal pos
            whitespace();assert pos<len(buf),'unexpected EOF'
            c=buf[pos];pos+=1;return c
        def value():
            nonlocal pos
            whitespace()
            while True:
                try:o,end=decoder.raw_decode(buf,pos);pos=end;return o
                except json.JSONDecodeError:
                    assert not eof,'bad JSON'
                    fill()
        assert take()=='{'
        while True:
            key=value();assert take()==':'
            if key=='histogram':
                assert not found;found=True;assert take()=='{'
                while True:
                    key=value();assert isinstance(key,str) and take()==':'
                    coeff=value();assert isinstance(coeff,dict)
                    yield key,coeff
                    end=take()
                    if end=='}':break
                    assert end==','
            else:value()
            end=take()
            if end=='}':break
            assert end==','
        whitespace();assert eof and pos==len(buf),'trailing JSON data'
        assert found


def main():
    reg=load(HERE/'registration.json')
    assert sha(__file__)==reg['source_sha256']
    assert not (HERE/'results.json').exists()
    assert all(os.environ.get(k)=='1' for k in ENV)
    assert os.getpriority(os.PRIO_PROCESS,0)>=10
    source=Path(reg['source_root']);root=source.parent;pwa=source/'pwa'
    out={'status':'running','pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),
         'python':sys.executable,'source_sha256':sha(__file__),
         'registration_sha256':sha(HERE/'registration.json')}
    def progress(stage,**counts):
        out['last_complete_stage']=stage;out.update(counts)
        save('checkpoint.json',dict(out,seconds=time.monotonic()-START))
        print(json.dumps({'stage':stage,'pid':os.getpid(),**counts}),flush=True)
    try:
        for item in reg['readonly_inputs']:
            assert sha(item['path'])==item['sha256'],item['path']
            if 'bytes' in item:assert Path(item['path']).stat().st_size==item['bytes']
        ul,uh=root_box(2);vl,vh=root_box(3)
        for i,j in product(range(5),repeat=2):MB[i,j]=(ul**i*vl**j,uh**i*vh**j)
        a=add(ONE,scale(power(U,3),Q(1,2)));s=scale(add(ONE,U),Q(1,2));threea=scale(a,3)
        geom=load(root/'lower_round025/tail_geometry.json');families={};expected_raw=[]
        assert len(geom['raw_inventory'])==1024
        oldfam={r['family_id']:r for r in geom['families']}
        for fid,lengths in enumerate(product(range(1,5),repeat=5)):
            n=prod(lengths);r=radius(n);tail=[]
            for ell in lengths[1:]:
                j=-((ell-1)//2)
                tail.extend((pmax(ZERO,sub(scale(s,j+ell-1),r)),pmin(scale(s,Q(1,2)),add(scale(s,j),r))))
            positive=all(cmp(tail[k+1],tail[k])>0 for k in range(0,8,2))
            row=geom['raw_inventory'][fid]
            assert row['family_id']==fid and row['lengths']==list(lengths) and row['N']==n
            assert row['status']==('positive_tail_box' if positive else 'tail_empty_or_fixed_zero_width')
            if positive:
                old=oldfam[fid];assert old['lengths']==list(lengths) and old['N']==n
                assert dec(old['radius'])==r and [dec(x) for x in old['tail_box']]==tail
                assert scale(power(r,5),32)==mono(0,0,n)
                assert cmp(r,mono(0,0,2))<=0 and cmp(r,scale(a,lengths[0]))<0
                if lengths[0]==4:assert cmp(threea,mono(0,0,4))>0
                families[fid]=(lengths,n,r)
            for j in range(-1,5):
                expected_raw.append({'family_id':fid,'start':j,'first_phase':j%3,'last_phase':(j+lengths[0]-1)%3,
                                     'status':'fixed_tail_empty_or_zero' if not positive else
                                              'ell1_4_span3a_global_empty' if lengths[0]==4 else 'retained_for_entire_window'})
        assert len(families)==len(oldfam)==213
        cap=load(source/'capture_forms.json');assert cap['all6144_raw_inventory']==expected_raw
        assert cap['global_ell1_4_excluded_families']==[fid for fid,(ls,n,r) in families.items() if ls[0]==4]
        assert cap['counts']==dict(raw=6144,fixed_tail_empty_or_zero=4866,
                                   ell1_4_span3a_global_empty=540,retained_for_entire_window=738)
        forms=[(dec(x['constant']),x['slope']) for x in cap['forms']]
        expected_forms={(ZERO,0),(threea,0)};z=(1,0,-1)
        for fid,(lengths,n,r) in families.items():
            ell=lengths[0]
            if ell==4:continue
            for j in range(-1,5):
                expected_forms.add((sub(scale(a,j+ell-1),r),z[(j+ell-1)%3]))
                expected_forms.add((add(scale(a,j),r),z[j%3]))
        assert len(forms)==len(set(forms))==279 and set(forms)==expected_forms
        eventsfile=load(source/'events.json');events=eventsfile['events'];E=[dec(x['value']) for x in events]
        assert len(E)==len(set(E))==965 and E[0]==ZERO and E[-1]==a
        for i,row in enumerate(events):
            assert row['event_index']==i;inside_saved(E[i],row['isolating_interval'])
            if i:assert cmp(E[i],E[i-1])>0
            if 'witness_form_indices' in row:
                f,g=[forms[k] for k in row['witness_form_indices']]
                assert scale(sub(g[0],f[0]),Q(1,f[1]-g[1]))==E[i]
        seen={ZERO,a};pairs=0
        for (c,m),(d,n) in combinations(forms,2):
            pairs+=1
            if m!=n:
                e=scale(sub(d,c),Q(1,m-n))
                if cmp(e)>=0 and cmp(e,a)<=0:seen.add(e)
        assert pairs==38781 and seen==set(E)
        progress('full_geometry_raw_forms_events',tail_shapes=213,raw_descriptors=6144,form_pairs=pairs,closed_events=965)

        # Read existing mask weights, independently compress by N groups (not set-bit walking).
        groups=[];allowed=sum(1<<fid for fid in families)
        for ell in (1,2,3):
            byN={}
            for fid,(lengths,n,r) in families.items():
                if lengths[0]==ell:byN[n]=byN.get(n,0)|(1<<fid)
            groups.append(sorted(byN.items(),reverse=True))
        aggregated={};maskseen=set();maskrows=0
        for key,coeff in stream_histogram(source/'full_tail_histogram.json.gz'):
            mask=int(key);assert 0<=mask and not(mask&~allowed) and mask not in maskseen
            maskseen.add(mask);desc=tuple(next((n for n,bits in group if mask&bits),0) for group in groups)
            target=aggregated.setdefault(desc,{})
            for ij,c in dec(coeff):target[ij]=target.get(ij,Q(0))+c
            maskrows+=1;assert maskrows<=139843
            if maskrows%40000==0:print(json.dumps({'mask_rows_read':maskrows,'pid':os.getpid()}),flush=True)
        assert maskrows==139843
        del maskseen
        saved=load(source/'descriptors.json.gz')['descriptors'];interface=load(pwa/'interface.json.gz')
        descs=interface['descriptors'];assert len(saved)==len(descs)==len(aggregated)==91
        assert {tuple(row['max_N_by_ell1']) for row in saved}==set(aggregated)
        fmap={row:i for i,row in enumerate(forms)};weights=[]
        for row,iface in zip(saved,descs):
            key=tuple(row['max_N_by_ell1']);weight=canonical(aggregated[key])
            assert weight==dec(row['physical_tail_weight'])==dec(iface['physical_tail_weight']) and cmp(weight)>0
            assert iface['max_N_by_ell1']==list(key)
            expected=[]
            for ell,n in enumerate(key,1):
                if not n:continue
                r=radius(n)
                for j in range(-1,5):expected.append([fmap[(sub(scale(a,j+ell-1),r),z[(j+ell-1)%3])],
                                                    fmap[(add(scale(a,j),r),z[j%3])]])
            assert iface['raw_interval_form_pairs']==expected
            weights.append(weight)
        assert add(*weights)==power(scale(s,Q(1,2)),4)
        assert interface['zero_form']==fmap[ZERO,0] and interface['domain_form']==fmap[threea,0]
        progress('static_archive_reaggregated',static_masks_read=maskrows,descriptor_weights=91,static_cells_rebuilt=0)

        C=dec(load(root/'lower_round003/a2_equals_a3_results.json')['coefficients'])
        oldD=dec(load(root/'lower_round025/plus_derivative.json')['directional_derivative_coefficients'])
        first=load(source/'first_cell_certificate.json');inherited=load(pwa/'cell0000_inherited.json')
        assert cmp(oldD)<0 and dec(first['D_coefficients'])==oldD
        assert dec(inherited['ratio_constant'])==C and dec(inherited['ratio_slope'])==oldD
        assert inherited['source_certificate_sha256']==sha(source/'first_cell_certificate.json')
        assert dec(first['right_epsilon_coefficients'])==E[1]
        assert mul(oldD,E[1])==dec(first['right_contact_minus_C_coefficients'])
        formbounds=[bounds(c) for c,m in forms];closed=[];nextcell=0;weightedrows=0;order_checks=0
        previous=None;continuity=0;residual_digest=hashlib.sha256()
        preg=load(pwa/'registration.json');producer_sha=preg['source_sha256'];interface_sha=sha(pwa/'interface.json.gz')
        def cellcheck(cell):
            nonlocal nextcell,weightedrows,order_checks,previous,continuity
            i=cell['cell_index'];assert i==nextcell and cell['left_event']==i and cell['right_event']==i+1
            A=dec(cell['ratio_constant']);B=dec(cell['ratio_slope'])
            assert cell['source_sha256']==producer_sha
            if i:
                q=Q(cell['ordering_only_interior_dyadic'])
                assert bounds(E[i])[1]<q<bounds(E[i+1])[0]
                lows=[lo+m*q for (lo,hi),(c,m) in zip(formbounds,forms)]
                highs=[hi+m*q for (lo,hi),(c,m) in zip(formbounds,forms)]
                order=sorted(range(279),key=lambda j:lows[j])
                assert all(highs[j]<lows[k] for j,k in zip(order,order[1:]))
                assert hashlib.sha256(json.dumps(order,separators=(',',':')).encode()).hexdigest()==cell['ordering_sha256']
                order_checks+=279
                rows=cell['descriptor_responses'];assert len(rows)==cell['descriptor_dual_unions']==91
                AA=ZERO;BB=ZERO
                for j,row in enumerate(rows):
                    assert row['descriptor']==j and isinstance(row['length_slope'],int)
                    AA=add(AA,mul(weights[j],dec(row['length_constant'])))
                    BB=add(BB,scale(weights[j],row['length_slope']));weightedrows+=1
                assert scale(AA,Q(16,3))==A and scale(BB,Q(16,3))==B
                assert cell['interface_sha256']==interface_sha
            residuals=[add(sub(A,C),mul(B,E[e])) for e in (i,i+1)]
            for value,key in zip(residuals,('left','right')):
                assert value==dec(cell[key+'_residual']);inside_saved(value,cell[key+'_residual_interval'])
            if previous is None:closed.append(residuals[0])
            else:assert previous==residuals[0];continuity+=1
            previous=residuals[1];closed.append(previous)
            choose=1 if cmp(residuals[1],residuals[0])>0 else 0
            assert residuals[choose]==dec(cell['cell_max_residual']) and cell['cell_max_event']==i+choose
            inside_saved(residuals[choose],cell['cell_max_residual_interval'])
            nextcell+=1
        cellcheck(inherited)
        previous_receipt=None;receipt_pids=[];batch_cells=0
        receipts=sorted((pwa/'batches').glob('batch_*.receipt.json'));assert len(receipts)==16
        for receipt_path in receipts:
            receipt=load(receipt_path);brpath=receipt_path.with_name(receipt_path.name.replace('.receipt.json','.registration.json'))
            br=load(brpath);assert sha(brpath)==receipt['batch_registration_sha256']
            assert br['source_sha256']==receipt['source_sha256']==producer_sha
            assert br['interface_sha256']==interface_sha
            assert br['first_cell']==receipt['first_cell']==nextcell
            assert br['last_cell']==receipt['last_requested_cell'] and br['last_cell']-br['first_cell']+1<=64
            if previous_receipt:assert br['previous_receipt_sha256']==previous_receipt
            artifact=pwa/'batches'/receipt['cell_artifact'];assert sha(artifact)==receipt['cell_artifact_sha256']
            cells=load(artifact)['cells'];ids=list(range(br['first_cell'],br['last_cell']+1))
            assert [c['cell_index'] for c in cells]==receipt['completed_cell_indices']==ids
            for cell in cells:
                assert cell['batch_registration_sha256']==sha(brpath);cellcheck(cell)
            assert receipt['new_cells_computed_this_process']==len(ids)
            assert receipt['new_descriptor_dual_union_objects']==91*len(ids)
            assert receipt['covered_prefix_cells']==nextcell and receipt['remaining_cells']==964-nextcell
            receipt_pids.append(receipt['runtime']['pid']);batch_cells+=len(ids);previous_receipt=sha(receipt_path)
            progress('cell_archive_prefix',completed_cells=nextcell,weighted_response_rows=weightedrows)
        assert nextcell==964 and len(closed)==965 and continuity==963 and weightedrows==87633
        assert order_checks==963*279 and batch_cells==963 and len(set(receipt_pids))==16
        zero=[];negative=[];positive=[]
        for i,r in enumerate(closed):
            sign=cmp(r)
            (zero if sign==0 else negative if sign<0 else positive).append(i)
            residual_digest.update((json.dumps(enc(r),sort_keys=True,separators=(',',':'))+'\n').encode())
        assert zero==[0] and len(negative)==964 and not positive
        best=1
        for i in range(2,965):
            if cmp(closed[i],closed[best])>0:best=i
        window=load(pwa/'window_result.json');state=load(pwa/'state.json')
        assert window['closed_events']==965 and window['open_cells']==964
        assert window['zero_residual_event_indices']==zero and window['strict_negative_closed_events']==964
        assert window['strict_positive_closed_events']==[] and window['global_maximum_event']==0
        assert dec(window['global_maximum_residual_coefficients'])==ZERO
        assert dec(window['global_maximum_ratio_coefficients'])==C
        assert window['best_noninitial_contact']['event']==best==1
        assert dec(window['best_noninitial_contact']['residual_coefficients'])==closed[best]
        inside_saved(add(C,closed[-1]),window['collision_limit_ratio_interval'])
        assert state['next_cell']==964 and state['new_cells_computed']==963 and state['inherited_cells']==1
        assert dec(state['best_global_residual'])==ZERO
        assert state['interface_sha256']==interface_sha
        for item in reg['readonly_inputs']:assert sha(item['path'])==item['sha256'],item['path']
        out.update(status='passed_independent_algebra_readback_reusing_static_and_union_responses',
                   closed_events=965,closed_negative=964,closed_zero=[0],open_cells=964,
                   adjacent_exact_continuities=963,weighted_response_rows=weightedrows,
                   endpoint_positions_order_checked=order_checks,descriptor_weights_positive=91,
                   producer_batch_pids=receipt_pids,closed_residual_digest=residual_digest.hexdigest(),
                   root_denominator=str(1<<256),root_lower_numerators=[str(int(ul*(1<<256))),str(int(vl*(1<<256)))],
                   root_lower_fractions=[str(ul),str(vl)],old_D_upper=str(bounds(oldD)[1]),
                   best_noninitial_event=best,best_noninitial_residual_interval=list(map(str,bounds(closed[best]))),
                   collision_limit_ratio_interval=list(map(str,bounds(add(C,closed[-1])))),
                   true_union_calls=0,static_cells_rebuilt=0,old_baseline_calls=0,
                   static_weight_geometry_independently_rebuilt=False,
                   descriptor_response_geometry_independently_rebuilt=False,
                   producer_math_modules_imported=False,uniform_positive_epsilon_gap_claim=False,
                   C5_upper_bound_claim=False,formal_bound_changed=False)
    except Exception as error:
        out.update(status='stopped_failure_no_rerun',error=repr(error))
        raise
    finally:
        out['seconds_before_final_write']=time.monotonic()-START
        out['peak_rss_platform_units']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        save('results.json',out)
    print(json.dumps({k:out[k] for k in ('status','pid','closed_events','closed_negative','seconds_before_final_write','peak_rss_platform_units')}),flush=True)


if __name__=='__main__':main()
