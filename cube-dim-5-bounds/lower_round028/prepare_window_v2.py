#!/usr/bin/env python3
"""Whole-window raw capture events and a swept-interval upper relaxation only."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,combinations
from collections import defaultdict,Counter
import sys,os,time,sqlite3,json,gzip,hashlib
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent;OLD=ROOT/'lower_round025'
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from full_static import ENV,load,save,sha,rt,decode,encode,mul,addto,ZERO

def main(stage):
    start=time.monotonic();reg=load(P/'window_registration_v2.json');assert sha(Path(__file__))==reg['source_sha256']
    for row in reg['inputs']:assert sha(Path(row['path']))==row['sha256']
    F,u,v,bounds,_,_=field();zero=F(0);a=1+u**3/2;s=(1+u)/2
    def dec(z):return F(a=decode(z))
    geom=load(OLD/'tail_geometry.json');fams={r['family_id']:r for r in geom['families']};radii={fid:dec(r['radius']) for fid,r in fams.items()}
    z=(1,0,-1);db=sqlite3.connect('file:'+reg['scratch_database']+'?mode=ro',uri=True)
    meta=dict(db.execute('SELECT key,value FROM meta'));assert meta['completed_slabs']=='26'
    if stage=='plan':
        assert not (P/'window_plan.json').exists(),'windowinputplanalreadycomplete'
        forms={(zero,0),(3*a,0)};raw=[];counts=Counter();global4=[]
        for row in geom['raw_inventory']:
            fid=row['family_id'];ell=row['lengths'][0]
            if fid in fams:
                r=radii[fid];assert r<=F(2) and r<ell*a
                if ell==4:assert 2*r<=F(4)<3*a;global4.append(fid)
            for j in range(-1,5):
                counts['raw']+=1;status='fixed_tail_empty_or_zero' if fid not in fams else 'ell1_4_span3a_global_empty' if ell==4 else 'retained_for_entire_window'
                counts[status]+=1
                raw.append({'family_id':fid,'start':j,'first_phase':j%3,'last_phase':(j+ell-1)%3,'status':status})
                if status=='retained_for_entire_window':
                    r=radii[fid];forms.add(((j+ell-1)*a-r,z[(j+ell-1)%3]));forms.add((j*a+r,z[j%3]))
        assert counts['raw']==6144
        forms=sorted(forms);formrows=[{'constant':encode(c.a),'slope':m} for c,m in forms]
        save(P/'capture_forms.json',{'status':'COMPLETE_wholewindow_raw_capture_forms','forms':formrows,'all6144_raw_inventory':raw,'global_ell1_4_excluded_families':global4,'counts':dict(counts),'safe_start_proof':'j<=-2:p_j<=-2a andr<=2<2a;j>=5:lastindex>=ell1+4,p_last>(ell1+3)a,r<ell1*a,solower>3a','closed_upper_only_at_collision':'epsilon=a allowedonlyascontinuousupperendpoint,notdistinctunitatomwitness','source_sha256':reg['source_sha256']})
        # Exact descriptor compression: for each first-axis length, the largest
        # active radius dominates all smaller radii for everyepsilon andstart.
        descriptors={};streamed=0
        for maskstr,coeff in db.execute('SELECT mask,coefficients FROM hist'):
            bits=int(maskstr);best=[0,0,0]
            while bits:
                bit=bits&-bits;fid=bit.bit_length()-1;bits-=bit
                row=fams[fid];ell=row['lengths'][0]
                if ell<4:best[ell-1]=max(best[ell-1],row['N'])
            key=tuple(best)
            if key not in descriptors:descriptors[key]=[ZERO]*25
            addto(descriptors[key],decode(json.loads(coeff)));streamed+=1
        dp=P/'descriptors.json.gz'
        with dp.open('wb') as rawout:
            with gzip.GzipFile(fileobj=rawout,mode='wb',mtime=0) as out:out.write((json.dumps({'descriptors':[{'max_N_by_ell1':list(key),'physical_tail_weight':encode(w)} for key,w in sorted(descriptors.items())]},separators=(',',':'))+'\n').encode())
        # Inventory the actual independent2D-pair verification object, no4D AND.
        axes=[[dec(c) for c in axis] for axis in geom['axes']];boxes={fid:[dec(c) for c in r['tail_box']] for fid,r in fams.items()};segs=[]
        for k,axis in enumerate(axes):segs.append([((hi-lo).a,sum(1<<fid for fid,b in boxes.items() if b[2*k]<=lo and hi<=b[2*k+1])) for lo,hi in zip(axis,axis[1:])])
        pairhist=[]
        for offset in (0,2):
            h={}
            for (dx,mx),(dy,my) in product(segs[offset],segs[offset+1]):
                key=mx&my;value=mul(dx,dy)
                if key not in h:h[key]=value
                else:addto(h[key],value)
            pairhist.append(h)
        save(P/'pair_plan.json',{'status':'COMPLETE_new213_pair_input_inventory_no4Dverificationyet','twoD_cells_each':676,'signature_counts':list(map(len,pairhist)),'full_independent_AND_pairs':len(pairhist[0])*len(pairhist[1]),'new_independent_AND_pairs_executed':0,'histograms':[{str(m):encode(w) for m,w in h.items()} for h in pairhist],'source_sha256':reg['source_sha256']})
        plan={'status':'COMPLETE_wholewindow_static_descriptor_and_capture_plan','all_tail_families':213,'static_signatures_streamed':streamed,'raw_boxes':6144,'counts':dict(counts),'distinct_capture_affine_forms':len(forms),'all_form_pairs':len(forms)*(len(forms)-1)//2,'different_slope_pairs':sum(m!=n for (_,m),(_,n) in combinations(forms,2)),'descriptor_count':len(descriptors),'capture_forms_sha256':sha(P/'capture_forms.json'),'descriptors_sha256':sha(dp),'pair_plan_sha256':sha(P/'pair_plan.json'),'new_epsilon_union_evaluations':0,'runtime':rt(start),'source_sha256':reg['source_sha256']}
        save(P/'window_plan.json',plan);print(plan,flush=True)
    elif stage=='events':
        er=load(P/'event_registration.json');assert er['source_sha256']==reg['source_sha256'] and er['capture_forms_sha256']==sha(P/'capture_forms.json')
        if (P/'events.json').exists():raise ValueError('completeeventinventoryalreadyexists')
        forms=[(dec(r['constant']),r['slope']) for r in load(P/'capture_forms.json')['forms']]
        found={zero.a:{'value':encode(zero.a),'witness':'window_left'},a.a:{'value':encode(a.a),'witness':'window_collision_limit'}};done=0;cp=P/'event_checkpoint.json'
        if cp.exists():
            ck=load(cp);assert ck['source_sha256']==reg['source_sha256']
            done=ck['completed_pairs'];found={tuple(decode(r['value'])):r for r in ck['events']}
        initial=done;cursor=0
        def checkpoint(status):save(cp,{'status':status,'completed_pairs':done,'total_pairs':len(forms)*(len(forms)-1)//2,'events':list(found.values()),'source_sha256':reg['source_sha256'],'runtime':rt(start)})
        print({'stage':'finite_capture_event_pairs','pid':os.getpid(),'forms':len(forms),'all_pairs':len(forms)*(len(forms)-1)//2,'initial_prefix':initial},flush=True)
        for i,(c,m) in enumerate(forms):
            for j in range(i+1,len(forms)):
                if cursor<initial:cursor+=1;continue
                if done-initial>=er['max_new_pairs_per_process']:checkpoint('FINITE_PAIR_LIMIT');return
                d,n=forms[j]
                if m!=n:
                    eps=(d-c)/(m-n)
                    if zero<=eps<=a:found.setdefault(eps.a,{'value':encode(eps.a),'witness_form_indices':[i,j]})
                done+=1;cursor+=1
                if done%10000==0:checkpoint('EXACT_pair_prefix_saved')
        ordered=sorted((F(a=k),row) for k,row in found.items())
        rows=[]
        for index,(value,row) in enumerate(ordered):rows.append(dict(row,event_index=index,isolating_interval=list(map(str,bounds(value.a)))))
        result={'status':'COMPLETE_all_wholewindow_capture_topology_events','affine_forms':len(forms),'pairs_checked':done,'closed_event_count':len(rows),'open_cells':len(rows)-1,'events':rows,'capture_forms_sha256':sha(P/'capture_forms.json'),'source_sha256':reg['source_sha256'],'new_epsilon_union_evaluations':0,'runtime':rt(start),'scope':'safecompletearrangementpossiblyoverrefined;capture/clipping/zero-width/touchchangesincluded;collisionendpointupperonly'}
        save(P/'events.json',result);checkpoint('COMPLETE_event_inventory');print({k:v for k,v in result.items() if k!='events'},flush=True)
    elif stage=='relax':
        rr=load(P/'relaxation_registration.json');assert rr['source_sha256']==reg['source_sha256'] and rr['descriptors_sha256']==sha(P/'descriptors.json.gz')
        if (P/'relaxation_result.json').exists():raise ValueError('completerelaxationalreadyexists')
        with gzip.open(P/'descriptors.json.gz','rt') as f:rows=json.load(f)['descriptors']
        assert len(rows)==rr['descriptor_count']<=rr['max_descriptors']
        Nroots={r['N']:radii[fid] for fid,r in fams.items()};total=[ZERO]*25;prefix=0;cp=P/'relaxation_checkpoint.json';outrows=[]
        if cp.exists():
            ck=load(cp);assert ck['source_sha256']==reg['source_sha256'];prefix=ck['completed_descriptors'];total=decode(ck['partial_upper_volume']);outrows=ck['rows']
        def checkpoint(status):save(cp,{'status':status,'completed_descriptors':prefix,'total_descriptors':len(rows),'partial_upper_volume':encode(total),'rows':outrows,'source_sha256':reg['source_sha256'],'runtime':rt(start)})
        for index in range(prefix,len(rows)):
            row=rows[index];intervals=[]
            for ell,N in enumerate(row['max_N_by_ell1'],1):
                if not N:continue
                r=Nroots[N]
                for j in range(-1,5):
                    c=(j+ell-1)*a-r;v=z[(j+ell-1)%3];d=j*a+r;w=z[j%3]
                    lo=max(zero,min(c,c+v*a));hi=min(3*a,max(d,d+w*a))
                    if lo<hi:intervals.append((lo,hi))
            length=F(0);left=right=None
            for lo,hi in sorted(intervals):
                if left is None:left,right=lo,hi
                elif right<lo:length+=right-left;left,right=lo,hi
                else:right=max(right,hi)
            if left is not None:length+=right-left
            assert zero<=length<=3*a
            addto(total,mul(decode(row['physical_tail_weight']),length.a));outrows.append({'max_N_by_ell1':row['max_N_by_ell1'],'swept_interval_union_length_upper':encode(length.a)})
            prefix=index+1
            if prefix%128==0:checkpoint('EXACT_descriptor_upper_prefix_saved')
        volume=F(a=total);ratio=volume*Q(16,3);C=load(ROOT/'lower_round003/a2_equals_a3_results.json');residual=ratio-dec(C['coefficients']);lo,hi=bounds(residual.a)
        result={'status':'COMPLETE_safe_sectionwise_window_upper_relaxation','upper_ratio_coefficients':encode(ratio.a),'upper_ratio_interval':list(map(str,bounds(ratio.a))),'upper_minus_complete_C_coefficients':encode(residual.a),'upper_minus_C_interval':[str(lo),str(hi)],'upper_ratio_display_only':[float(x) for x in bounds(ratio.a)],'strictly_below_C_entire_window':hi<0,'relaxation_too_coarse_to_exclude_window':lo>0,'descriptor_count':len(rows),'rule':'for eachrawintervaltakeindependentminlower/maxupperoverclosed[0,a],uniontheseprojectedintervalswithinsection,thenweightsum;validuppernotrealizablewitnessratio','closed_collision_endpoint_used_only_as_upper':True,'new_actual_epsilon_union_evaluations':0,'no_unique_contact_selected':True,'source_sha256':reg['source_sha256'],'runtime':rt(start)}
        save(P/'relaxation_result.json',result);checkpoint('COMPLETE_window_upper_relaxation');print(result,flush=True)
    else:raise ValueError('stage plan|events|relax')
    db.close()
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    if len(sys.argv)!=2:raise SystemExit('usage: prepare_window.py plan|events|relax')
    main(sys.argv[1])
