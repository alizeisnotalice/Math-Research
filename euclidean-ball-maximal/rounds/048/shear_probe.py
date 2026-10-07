#!/usr/bin/env python3
"""Exact n1 coarse-source shear of the full common-record outer measure.
All source atoms retain their original weights. Rational polygon clipping and
affine pushforward cross-sections replace sampling grids.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from bisect import bisect_left,bisect_right
import argparse,importlib.util,json,time,hashlib,sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
LABELS={'minimal_four_atoms','four_atom_perturbation_s471016','geometric_contact_clouds_s491101','geometric_contact_clouds_s491108','chain_L4_alpha3/2','chain_L4_alpha2','chain_L8_alpha2','chain_L16_alpha2','minimal_cloud_with_two_coarse_atoms'}

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def clip(poly,A,B,C):
    """Closed half-plane Ax+Bz<=C; boundary has zero area."""
    if not poly:return []
    out=[]
    for p,q in zip(poly,poly[1:]+poly[:1]):
        vp=A*p[0]+B*p[1]-C;vq=A*q[0]+B*q[1]-C
        if vp<=0:out.append(p)
        if (vp<0<vq) or (vq<0<vp):
            t=vp/(vp-vq);out.append((p[0]+t*(q[0]-p[0]),p[1]+t*(q[1]-p[1])))
    clean=[]
    for p in out:
        if not clean or p!=clean[-1]:clean.append(p)
    if len(clean)>1 and clean[0]==clean[-1]:clean.pop()
    return clean

def area(poly):
    if len(poly)<3:return F(0)
    return abs(sum((p[0]*q[1]-p[1]*q[0] for p,q in zip(poly,poly[1:]+poly[:1])),F(0)))/2

def slab(poly,y,l,r):
    return clip(clip(poly,F(-1),F(-1),-l-y),F(1),F(1),r+y)

def x_slab(poly,y,R):
    return clip(clip(poly,F(-1),F(0),-y+R),F(1),F(0),y+R)

def cross_section(poly,s):
    vals=[]
    for p,q in zip(poly,poly[1:]+poly[:1]):
        a=p[0]+p[1];b=q[0]+q[1]
        if a==b:
            if s==a:vals.extend((p[1],q[1]))
        elif min(a,b)<=s<=max(a,b):
            vals.append(p[1]+(q[1]-p[1])*(s-a)/(b-a))
    return max(vals)-min(vals) if vals else F(0)

def atomic_maximal(atoms,x):
    distances=sorted({abs(x-y) for y,w in atoms})
    if distances[0]==0:return None
    return max(sum((w for y,w in atoms if abs(x-y)<=d),F(0))/(2*d) for d in distances)

def escape_certificate(atoms,alpha,radii,witness):
    poly=[(F(x),F(z)) for x,z in witness['polygon']]
    x=sum((p[0] for p in poly),F(0))/len(poly);z=sum((p[1] for p in poly),F(0))/len(poly)
    y=F(witness['source_y']);q=x+z-y;beta=(F(witness['beta_lo'])+F(witness['beta_hi']))/2
    def trace(t):return [sum((w for s,w in atoms if abs(t-s)<R),F(0))/(2*R) for R in radii[1:]]
    tx=trace(x);tz=trace(z);tq=trace(q)
    kx=next(k for k,g in enumerate(tx,1) if g>beta);kz=next(k for k,g in enumerate(tz,1) if g>beta)
    assert kx==witness['j'] and kz==witness['k']
    assert abs(x-z)>radii[kx] and abs(x-y)<radii[kx] and abs(z-y)<radii[kz]
    maximal=atomic_maximal(atoms,q);kind=witness['target_kind']
    if kind=='MP_gt_2alpha':assert maximal is None or maximal>2*alpha
    elif kind=='MP_le_alpha':assert maximal is not None and maximal<=alpha
    elif kind=='band_J1':assert maximal is not None and alpha<maximal<=2*alpha and tq[0]>alpha/8
    witness['representative_point']=dict(x=str(x),z=str(z),source_y=str(y),q=str(q),beta=str(beta),K_beta_x=kx,K_beta_z=kz,target_maximal=str(maximal) if maximal is not None else 'infinity',target_maximal_over_alpha=str(maximal/alpha) if maximal is not None else 'infinity',target_g1=str(tq[0]),target_J1_threshold=str(alpha/8),original_gj_x=str(tx[kx-1]),original_gk_z=str(tz[kz-1]),status='passed')

def density_events(poly,y,coefficient,events):
    sums=sorted({x+z for x,z in poly});values=[cross_section(poly,s) for s in sums]
    assert values[0]==values[-1]==0
    assert sum(((b-a)*(u+v)/2 for a,b,u,v in zip(sums,sums[1:],values,values[1:])),F(0))==area(poly)
    for a,b,u,v in zip(sums,sums[1:],values,values[1:]):
        slope=coefficient*(v-u)/(b-a);events[a-y]+=slope;events[b-y]-=slope

def target_cells(old,atoms,alpha,radii,saved):
    full_E=old.observer.__globals__['full_E'];high=full_E(atoms,2*alpha);superlevel=full_E(atoms,alpha)
    eligible=[(F(c['lo']),F(c['hi'])) for c in saved['observer_cells']]
    cuts={y for y,w in atoms}
    for intervals in (high,superlevel,eligible):
        for a,b in intervals:cuts.update((a,b))
    for y,w in atoms:
        for R in radii[1:]:cuts.update((y-R,y+R))
    cuts=sorted(cuts);cells=[]
    def inside(intervals,x):return any(a<x<b for a,b in intervals)
    for a,b in zip(cuts,cuts[1:]):
        mid=(a+b)/2
        if inside(high,mid):kind='MP_gt_2alpha'
        elif not inside(superlevel,mid):kind='MP_le_alpha'
        elif inside(eligible,mid):kind='eligible_E'
        else:
            # These finite test cases cover the full band within the chosen depth.
            # Assert this instead of silently classifying a late exit as J=1.
            g1=sum((w for y,w in atoms if abs(y-mid)<radii[1]),F(0))/(2*radii[1])
            assert g1>alpha/8, 'unclassified eligibility loss: extend target partition'
            kind='band_J1'
        trace=None;records=None
        if kind=='eligible_E':
            trace=[sum((w for y,w in atoms if abs(y-mid)<R),F(0))/(2*R) for R in radii[1:]]
            assert trace[0]<=alpha/8 and max(trace)>alpha/2
        cells.append((a,b,kind,trace))
    return cuts,cells

def evaluate(flow,record,old,batch,label,atoms,alpha,params):
    began=time.monotonic();saved,radii,cells,_=flow.context(old,atoms,alpha);lo=alpha/4;hi=alpha/2;I=hi-lo
    targets,target=target_cells(old,atoms,alpha,radii,saved)
    prepared=[]
    for a,b,kind,trace in target:
        prepared.append((a,b,kind,record.intervals(trace,lo,hi) if trace is not None else None))
    data=[]
    for xa,xb,J,trace in cells:
        rec=record.intervals(trace,lo,hi);mid=(xa+xb)/2
        members={k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k]) for k in rec}
        data.append((xa,xb,rec,members))
    O=compatible=same=innerfine=overlapfine=F(0);classes={k:F(0) for k in ('eligible_E','MP_le_alpha','MP_gt_2alpha','band_J1')};relations={k:F(0) for k in ('Kq_lt_j','Kq_eq_j','j_lt_Kq_lt_k','Kq_eq_k','Kq_gt_k')};events=defaultdict(F)
    polygons=classified_pieces=beta_target_pairs=0;escape_witness=None;escape_by_kind={};eligible_bad_record_witness=None
    for xa,xb,ix,mx in data:
        for za,zb,iz,mz in data:
            rectangle=[(xa,za),(xb,za),(xb,zb),(xa,zb)];polygon_cache={}
            for j,(la,lb) in ix.items():
                R=radii[j]
                polygon_cache[j]=[p for p in (clip(rectangle,F(-1),F(1),-R),clip(rectangle,F(1),F(-1),-R)) if area(p)>0]
                square=lambda s:max(F(0),s)**2
                inside_area=sum((sign*(square(zb+shift)-square(za+shift))/2 for shift,sign in ((R-xa,1),(R-xb,-1),(-R-xa,-1),(-R-xb,1))),F(0))
                assert sum((area(p) for p in polygon_cache[j]),F(0))==(xb-xa)*(zb-za)-inside_area
                for k,(ra,rb) in iz.items():
                    if j>=k:continue
                    beta_l=max(la,ra);beta_r=min(lb,rb)
                    if beta_l>=beta_r:continue
                    r=radii[k];weight_beta=(beta_r-beta_l)/I
                    for source in mx[j]&mz[k]:
                        y,w=atoms[source];density=w/(4*R*r)
                        for poly in polygon_cache[j]:
                            polygons+=1;A=area(poly);cost=A*density*weight_beta;O+=cost
                            density_events(poly,y,density*weight_beta,events)
                            finepoly=x_slab(poly,y,r);finearea=area(finepoly);innerfine+=finearea*density*weight_beta
                            # In n1 an outer coarse edge has no inner-fine source.
                            assert finearea==0
                            tmin=min(x+z-y for x,z in poly);tmax=max(x+z-y for x,z in poly)
                            start=max(0,bisect_right(targets,tmin)-1);stop=bisect_left(targets,tmax)
                            subtotal=F(0)
                            for ta,tb,kind,iq in prepared[start:stop]:
                                if tb<=tmin or ta>=tmax:continue
                                sub=slab(poly,y,ta,tb);subarea=area(sub)
                                if not subarea:continue
                                classified_pieces+=1;subcost=subarea*density*weight_beta;classes[kind]+=subcost;subtotal+=subcost
                                if kind!='eligible_E':
                                    if kind not in escape_by_kind:
                                        cert=dict(source_y=str(y),source_mass=str(w),j=j,k=k,beta_lo=str(beta_l),beta_hi=str(beta_r),x_interval=[str(xa),str(xb)],z_interval=[str(za),str(zb)],target_interval=[str(ta),str(tb)],target_kind=kind,positive_cost=str(subcost),polygon=[[str(x),str(z)] for x,z in sub])
                                        escape_by_kind[kind]=cert
                                        if escape_witness is None:escape_witness=cert
                                    continue
                                targetsubtotal=F(0)
                                for l,(ql,qr) in iq.items():
                                    a=max(beta_l,ql);b=min(beta_r,qr)
                                    if a>=b:continue
                                    beta_target_pairs+=1;factor=density*(b-a)/I;value=subarea*factor;targetsubtotal+=value
                                    relation='Kq_lt_j' if l<j else 'Kq_eq_j' if l==j else 'j_lt_Kq_lt_k' if l<k else 'Kq_eq_k' if l==k else 'Kq_gt_k'
                                    relations[relation]+=value
                                    if j<=l<k:
                                        compoly=x_slab(sub,y,radii[l]);comparea=area(compoly);compatible+=comparea*factor
                                        if l==j:same+=comparea*factor
                                        overlapfine+=area(x_slab(compoly,y,r))*factor
                                    elif eligible_bad_record_witness is None:
                                        eligible_bad_record_witness=dict(source_y=str(y),j=j,k=k,target_K=l,beta_lo=str(a),beta_hi=str(b),target_interval=[str(ta),str(tb)],positive_cost=str(value),polygon=[[str(x),str(z)] for x,z in sub])
                                assert targetsubtotal==subcost
                            assert subtotal==cost
    assert sum(classes.values(),F(0))==O and sum(relations.values(),F(0))==classes['eligible_E']
    knots=sorted(t for t,v in events.items() if v);slope=value=mass=F(0);maximum=F(0);maximum_t=None;previous=None
    for t in knots:
        if previous is not None:
            next_value=value+slope*(t-previous);mass+=(t-previous)*(value+next_value)/2;value=next_value
        assert value>=0
        if value>maximum:maximum=value;maximum_t=t
        slope+=events[t];previous=t
    assert slope==value==0 and mass==O
    X=alpha*F(saved['eligible_volume'])
    ref=flow.flow_probe(old,record,batch,label,atoms,alpha,params);B=F(ref['mean']['B']);M=F(ref['mean']['M'])
    assert compatible<=2*B and same<=B and innerfine<=M/2 and innerfine<=F(3,8)*X
    union=compatible+innerfine-overlapfine;remainder=O-union;assert 0<=union<=O and remainder>=0
    for cert in escape_by_kind.values():escape_certificate(atoms,alpha,radii,cert)
    return dict(batch=batch,label=label,params=params,N=len(atoms),alpha=str(alpha),X=str(X),O=str(O),O_over_X=str(O/X),classes={k:str(v) for k,v in classes.items()},class_fraction={k:str(v/O) if O else None for k,v in classes.items()},eligible_record_relation={k:str(v) for k,v in relations.items()},compatible=str(compatible),same_label=str(same),inner_fine=str(innerfine),compatible_inner_overlap=str(overlapfine),paid_union=str(union),remainder=str(remainder),compatible_over_O=str(compatible/O) if O else None,remainder_over_O=str(remainder/O) if O else None,Bbar=str(B),Mbar=str(M),compatible_bound=str(2*B),same_label_bound=str(B),inner_fine_bound=str(M/2),pushforward_density_max=str(maximum),pushforward_density_max_over_alpha=str(maximum/alpha),pushforward_density_max_location=str(maximum_t) if maximum_t is not None else None,pushforward_density_knot_count=len(knots),pushforward_density_integral=str(mass),source_weight_preserved=True,pushforward_density_method='exact affine cross-sections of rational outer polygons; supremum at combined vertex-projection knots',positive_outer_source_polygons=polygons,target_classification_pieces=classified_pieces,beta_target_record_pairs=beta_target_pairs,escape_witness=escape_witness,escape_witness_by_kind=escape_by_kind,eligible_bad_record_witness=eligible_bad_record_witness,atoms=[dict(location=str(y),mass=str(w)) for y,w in atoms],exact_error='0',elapsed_seconds=time.monotonic()-began)

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);p.add_argument('--labels',nargs='*');args=p.parse_args()
    flow=module(args.repo_root/'rounds/047/flow_probe.py','flow48');record=module(args.repo_root/'rounds/046/threshold_probe.py','record48');prefix=module(args.repo_root/'rounds/046/prefix_probe.py','prefix48');old=record.load(args.repo_root,args.output)
    receipt_path=args.repo_root/'rounds/047/verification.json';reference={}
    if receipt_path.exists():
        receipt=json.loads(receipt_path.read_text());reference={c['label']:c['kernel']['positive_over_X'] for c in receipt['stress']['full_spatial_checks']}
    started=time.monotonic();out=dict(status='running',dimension=1,arithmetic='Fraction exact',scope='full original eligible observer, common beta, source-preserving shear; exact polygons and pushforward marginal',cases=[],finite_tests_prove_dimension_independent_density_cap=False)
    chosen=LABELS if args.labels is None else set(args.labels)
    inputs=list(flow.cases(old,prefix))
    inputs.append((0,'minimal_cloud_with_two_coarse_atoms',[(F(-507,256),F(3,10)),(F(0),F(13,256)),(F(1,16),F(3,64)),(F(3,32),F(9,128)),(F(515,256),F(3,10)),(F(10),F(297,1280))],F(1),dict(family='two_coarse_atoms_for_J1_target')))
    for batch,label,atoms,alpha,params in inputs:
        if label not in chosen:continue
        c=evaluate(flow,record,old,batch,label,atoms,alpha,params)
        if label in reference:assert F(c['O_over_X'])==F(reference[label]);c['archived_outer_kernel_crosscheck']='passed'
        out['cases'].append(c);args.output.write_text(json.dumps(out)+'\n')
        print(batch,label,'O/X',float(F(c['O_over_X'])),'classes',c['class_fraction'],'compatible/O',c['compatible_over_O'],'densitymax/alpha',float(F(c['pushforward_density_max_over_alpha'])),'seconds',round(c['elapsed_seconds'],2),flush=True)
    out.update(status='passed',case_count=len(out['cases']),batch_counts=[sum(c['batch']==b for c in out['cases']) for b in range(3)],elapsed_seconds=time.monotonic()-started,sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),exact_error='0');args.output.write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['case_count'],'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
