#!/usr/bin/env python3
"""Pure-rational finite-catalog guard for source-once greedy witness cover.

Finite catalog/receiver checks do not certify coverage of all translated boxes.
The general infinite-catalog statement belongs to the separately audited theory.
"""
from pathlib import Path
from fractions import Fraction as F
import datetime,hashlib,json,time

HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'greedy_witness_cover_guard_20261007_results.json'


def linf(x,y):return max(abs(a-b) for a,b in zip(x,y))
def vector(n,x=F(0),y=F(0)):
    z=[F(0)]*n;z[0]=F(x);z[1]=F(y);return tuple(z)
def in_box(z,center,side):return linf(z,center)<=side/2
def intersect(c,d,h):return linf(c,d)<=h
def touch(c,d,h):return linf(c,d)==h


def construct(n):
    a=F(n);h=F(1)
    atoms=[vector(n,F(-4,5)),vector(n,F(-2,5))]
    weights=[F(3,20),F(3,10)]
    if n==4:
        atoms.append(vector(n,F(2,5)));weights.append(F(7,20))
    else:
        atoms.extend([vector(n,F(2,5),F(1,4)),vector(n,F(2,5),F(-1,4))])
        weights.extend([F(1,5),F(3,20)])
    atoms.extend([vector(n,F(6,5)),vector(n,3*a+F(1,5))])
    weights.extend([F(1,20),F(3,20)])
    assert sum(weights)==1 and all(w>0 for w in weights)
    # Positive catalogs are first so the half-max policy can choose a smaller
    # representative. Exact max instead first chooses the negative .45 box.
    centers=[]
    for x in [F(3,10),F(1,2),F(9,10),F(-1,2),F(-3,10),F(-9,10),F(8,5),3*a+F(1,5),5*a]:
        for y in [F(0),F(1,2),F(-1,2)]:centers.append(vector(n,x,y))
    centers.extend([vector(n,F(-1,2),F(3,2)),vector(n,F(1,2),F(-3,2))])
    boxes=[]
    for i,c in enumerate(centers):
        members=[j for j,z in enumerate(atoms) if in_box(z,c,h)]
        boundary=[j for j in members if linf(atoms[j],c)==h/2]
        boxes.append(dict(id=i,center=c,members=members,mass=sum((weights[j] for j in members),F(0)),boundary_members=boundary))
    receivers=[]
    for q in [F(1,10),F(1,5),F(2,5)]:receivers.append(dict(name=f'positive_face_offset_{q}',x=vector(n,a/2+q),expected_hardband=True))
    for y in [F(1,4),F(-1,4)]:receivers.append(dict(name=f'positive_shift_y_{y}',x=vector(n,a/2+F(1,5),y),expected_hardband=True))
    for y in [F(0),F(1,4)]:receivers.append(dict(name=f'negative_arrival_y_{y}',x=vector(n,-a/2-F(3,5),y),expected_hardband=True))
    for coordinate in [1,n-1]:
        x=[F(0)]*n;x[coordinate]=a/2+F(1,2)
        receivers.append(dict(name=f'other_face_{coordinate}',x=tuple(x),expected_hardband=True))
    receivers.extend([dict(name='origin_above_band',x=vector(n),expected_hardband=False),
                      dict(name='outlier_below_band',x=vector(n,3*a+F(1,5)),expected_hardband=False),
                      dict(name='negative_at_side_a_above_band',x=vector(n,-a/2-F(2,5)),expected_hardband=False)])
    # One common lambda per complete source, fixed before receiver evaluation.
    lam=F(1)/(10*a**n)
    return a,h,atoms,weights,boxes,receivers,lam


def true_window_winner(x,atoms,weights,n,a):
    arrivals=[2*linf(x,z) for z in atoms]
    candidates=sorted({a,2*a}|{R for R in arrivals if a<=R<=2*a})
    rows=[]
    for R in candidates:
        members=[j for j,d in enumerate(arrivals) if d<=R]
        mass=sum((weights[j] for j in members),F(0))
        rows.append(dict(R=R,mass=mass,u=mass/R**n,members=members,
                         boundary_members=[j for j in members if arrivals[j]==R]))
    winner=max(rows,key=lambda r:(r['u'],-r['R']))
    # Independent accumulation by COMPLETE exact equal-arrival groups.
    groups={}
    for j,d in enumerate(arrivals):groups[d]=groups.get(d,F(0))+weights[j]
    cumulative=sum((mass for d,mass in groups.items() if d<=a),F(0))
    alternate=[(cumulative/a**n,a)]
    for d in sorted(groups):
        if a<d<=2*a:
            cumulative+=groups[d];alternate.append((cumulative/d**n,d))
    mass_b=sum((mass for d,mass in groups.items() if d<=2*a),F(0))
    alternate.append((mass_b/(2*a)**n,2*a))
    assert max(alternate,key=lambda r:(r[0],-r[1]))==(winner['u'],winner['R'])
    return winner,rows,arrivals


def greedy(boxes,h,v,mode):
    remaining=[b['id'] for b in boxes if b['mass']>v]
    selected=[];assignment={};steps=[]
    while remaining:
        maxmass=max(boxes[i]['mass'] for i in remaining)
        if mode=='exact_max':representative=next(i for i in remaining if boxes[i]['mass']==maxmass)
        else:representative=next(i for i in remaining if boxes[i]['mass']>=maxmass/2)
        rep=boxes[representative];removed=[i for i in remaining if intersect(boxes[i]['center'],rep['center'],h)]
        assert rep['mass']>=maxmass/2
        for i in removed:
            assert boxes[i]['mass']<=2*rep['mass']
            assignment[i]=representative
        selected.append(representative)
        steps.append(dict(representative=representative,representative_mass=rep['mass'],current_max_mass=maxmass,
                          removed=removed,touch_only_removed=[i for i in removed if i!=representative and touch(boxes[i]['center'],rep['center'],h)]))
        remaining=[i for i in remaining if i not in set(removed)]
    assert set(assignment)=={b['id'] for b in boxes if b['mass']>v}
    for i,x in enumerate(selected):
        for y in selected[i+1:]:assert not intersect(boxes[x]['center'],boxes[y]['center'],h)
    # Closed-touch boxes are removed, so no boundary atom gets charged twice.
    members=[j for i in selected for j in boxes[i]['members']]
    assert len(members)==len(set(members))
    return selected,assignment,steps


def check_instance(n,eta,mode):
    a,h,atoms,weights,boxes,receivers,lam=construct(n)
    W=sum(weights);v=2*eta*lam*a**n
    selected,assignment,steps=greedy(boxes,h,v,mode)
    selected_mass=sum((boxes[i]['mass'] for i in selected),F(0));assert selected_mass<=W
    factor=(1+3/F(n))**n
    cover_fee=factor*selected_mass/eta;ceiling=factor*W/eta;assert cover_fee<=ceiling
    received=[];witness_count=0;unseen_rep_count=0;unseen_rep_mass_count=0;R_above_a=0;closed_winner_count=0
    for receiver in receivers:
        x=receiver['x'];winner,candidates,arrivals=true_window_winner(x,atoms,weights,n,a)
        R,m,u=winner['R'],winner['mass'],winner['u']
        E=2*lam<u<=4*lam
        assert E==receiver['expected_hardband']
        closed_winner_count+=bool(winner['boundary_members'])
        R_above_a+=bool(E and R>a)
        rows=[]
        for B in boxes:
            captured=[j for j in B['members'] if j in winner['members']]
            mass_captured=sum((weights[j] for j in captured),F(0))
            qualified=E and mass_captured>eta*m
            row=dict(box=B['id'],captured_mass=mass_captured,captured_members=captured,qualified_witness=qualified)
            if qualified:
                witness_count+=1;assert B['mass']>v
                rep_id=assignment[B['id']];rep=boxes[rep_id];M=rep['mass']
                assert intersect(B['center'],rep['center'],h) and B['mass']<=2*M
                assert mass_captured>eta*m and m>2*lam*R**n
                assert R**n<M/(eta*lam)
                distance=linf(x,rep['center'])
                assert distance<=R/2+3*h/2
                assert R/2+3*h/2<=(1+3/F(n))*R/2
                # Radical-free membership in the fixed enlarged output cover.
                assert (2*distance/(1+3/F(n)))**n<M/(eta*lam)
                rep_in_Q=in_box(rep['center'],x,R)
                rep_captured_mass=sum((weights[j] for j in rep['members'] if j in winner['members']),F(0))
                unseen_rep_count+=not rep_in_Q;unseen_rep_mass_count+=rep_captured_mass==0
                row.update(representative=rep_id,representative_mass=M,witness_mass=B['mass'],
                    witness_to_representative_mass_ratio=B['mass']/M,
                    R_power=R**n,R_power_upper=M/(eta*lam),
                    receiver_to_rep_center_distance=distance,
                    geometric_distance_upper=R/2+3*h/2,
                    scale_dilation_distance_upper=(1+3/F(n))*R/2,
                    radical_free_output_cover_lhs=(2*distance/(1+3/F(n)))**n,
                    radical_free_output_cover_rhs=M/(eta*lam),
                    representative_center_in_original_Q=rep_in_Q,representative_mass_captured_by_original_Q=rep_captured_mass,
                    all_chain_checks=True)
            if B['mass']<=v:assert not qualified
            rows.append(row)
        received.append(dict(name=receiver['name'],x=x,hardband=E,winner=winner,candidates=candidates,
                         atom_arrival_sides=arrivals,witness_checks=rows))
    assert witness_count>0 and unseen_rep_count>0 and unseen_rep_mass_count>0 and R_above_a>0
    stats=dict(catalog_boxes=len(boxes),eligible_boxes=sum(b['mass']>v for b in boxes),
       small_or_equal_mass_boxes=sum(b['mass']<=v for b in boxes),mass_exactly_v_boxes=sum(b['mass']==v for b in boxes),
       selected_boxes=len(selected),source_mass_once=selected_mass,
       boundary_atom_box_incidences=sum(len(b['boundary_members']) for b in boxes),
       catalog_closed_touch_pairs=sum(touch(boxes[i]['center'],boxes[j]['center'],h) for i in range(len(boxes)) for j in range(i+1,len(boxes))),
       greedy_touch_removals=sum(len(s['touch_only_removed']) for s in steps),
       receivers=len(receivers),hardband_receivers=sum(r['hardband'] for r in received),
       hardband_winners_above_a=R_above_a,receivers_with_atoms_on_winner_boundary=closed_winner_count,
       qualified_receiver_witness_pairs=witness_count,assigned_rep_centers_outside_original_Q=unseen_rep_count,
       assigned_rep_source_mass_not_captured_by_original_Q=unseen_rep_mass_count)
    core=dict(n=n,a=a,b=2*a,h=h,W=W,eta=eta,common_lambda=lam,v=v,greedy_mode=mode,
        atoms=atoms,weights=weights,boxes=boxes,selected=selected,assignment=assignment,greedy_steps=steps,
        fixed_output_cover_volume_fee=cover_fee,general_source_once_ceiling=ceiling,
        exact_factor=(1+3/F(n))**n,receiver_records=received,counts=stats)
    core['input_sha256']=hashlib.sha256(json.dumps(dict(n=n,a=a,h=h,atoms=atoms,weights=weights,lambda_=lam,eta=eta,box_centers=[b['center'] for b in boxes]),default=str,sort_keys=True).encode()).hexdigest()
    return core


def main():
    start=time.perf_counter();results=[]
    for round_,n in enumerate([4,16,64],1):
        for eta in [F(1,4),F(1,2)]:
            for mode in ['exact_max','half_max_first']:
                result=check_instance(n,eta,mode);result['round']=round_;results.append(result)
                print('ROUND',round_,'n',n,'eta',eta,mode,result['counts'],flush=True)
    data=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),records=results,
       arithmetic='Fractions only for source weights, coordinates, box/cube membership, maximal responses, thresholds, greedy decisions, power bounds and cover volumes. No Monte Carlo, radius grid, numerical root or floating geometric comparison.',
       complete_source='All positive atomic masses retained with W=1; fixed common lambda for each complete n-dimensional source and all its receivers. Modes/etas never modify the source or lambda.',
       winner='Exact centered full-side continuous [a,2a] maximal winner: endpoints a,2a and every atom arrival side in the window; complete equal-distance groups; smaller side selected for response tie. Direct inclusion and independent grouped accumulation agree exactly.',
       greedy='Eligible catalog mu(B)>v=2eta lambda a^n. Select exact max or first box >=1/2 of remaining max; remove all CLOSED-intersecting boxes, including boundary-touching. Every eligible listed witness is assigned once; selected support boxes are pairwise disjoint and sum selected masses<=W.',
       cover='For a qualified receiver/witness, captured mass>eta*m and hardband u>2lambda give mu(B)>2eta lambda R^n. Since mu(B)<=2Mrep, R^n<Mrep/(eta lambda). Intersecting h-boxes and one captured atom give distance(x,crep)<=R/2+3h/2<=(1+3/n)R/2. Radical-free nth-power membership checks certify the fixed enlarged output cover. Sum lambda-volume <=(1+3/n)^n W/eta.',
       edge_cases='Exact mass=v catalog boxes are excluded, closed atom faces included, closed box contact removed, actual winners above a included, and selected representatives can have zero source mass captured by the original Qwinner.',
       finite_scope='Finite box catalog and finitely many exact receiver points certify only this implementation chain. Coverage of ALL arbitrary translated h-boxes, approximate maximum selection existence and the measurable infinite cover require the separately audited general theory; this catalog does not certify them.',
       construction_scope='Positive atomic components and overlapping catalog positions are mechanism-level B38/B41 inspiration only; no grid-density LP/MIP optimization, complete tensor lattice, A13/A14 or original geom/FIRST qualification is claimed.',
       skill_workflow='L02 finite feasible-input reconstruction and exact independent recomputation; no external optimization theorem invoked.',
       runtime_seconds=time.perf_counter()-start,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with OUTPUT.open('x') as f:json.dump(data,f,default=str,ensure_ascii=False,indent=2)
    print('OUTPUT',OUTPUT,'RUNTIME',data['runtime_seconds'],flush=True)


if __name__=='__main__':main()
