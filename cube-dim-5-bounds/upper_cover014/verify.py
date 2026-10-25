#!/usr/bin/env python3
"""Three registered exact instances; no maximal receiver sampling."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import argparse
import hashlib
import json
import time

SIGNS=tuple(product((-1,1),repeat=5))
ETA=F(1,65536)


def volume(box):
    answer=F(1)
    for lo,hi in box: answer*=max(F(0),hi-lo)
    return answer
def box(center,radius): return tuple((c-radius,c+radius) for c in center)
def intersect(*boxes): return tuple((max(b[j][0] for b in boxes),min(b[j][1] for b in boxes)) for j in range(5))
def inside(point,B): return all(lo<x<hi for x,(lo,hi) in zip(point,B))
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clique_geometry(cubes,common):
    assert len(cubes)==32
    core=intersect(*(Q['box'] for Q in cubes))
    assert volume(core)>0
    by_sign={tuple(1 if c>z else -1 if c<z else 0 for c,z in zip(Q['center'],common)):Q for Q in cubes}
    assert set(by_sign)==set(SIGNS)
    adjacent=0
    for sig,Q in by_sign.items():
        for j in range(5):
            if sig[j]==1: continue
            other=list(sig);other[j]=1;P=by_sign[tuple(other)]
            diff=abs(Q['center'][j]-P['center'][j])
            assert diff>max(Q['r'],P['r'])
            length=min(Q['box'][j][1],P['box'][j][1])-max(Q['box'][j][0],P['box'][j][0])
            assert 0<length<min(Q['r'],P['r'])
            adjacent+=1
    R=max(Q['r'] for Q in cubes)
    widths=[hi-lo for lo,hi in core]
    assert adjacent==80 and all(w<R for w in widths)
    return core,{'adjacent_orthant_pairs_checked':adjacent,'R':str(R),
                 'r_min':str(min(Q['r'] for Q in cubes)),
                 'core_box':[[str(x) for x in interval] for interval in core],
                 'core_widths':list(map(str,widths))}


def family_check(cubes):
    ratios=[max(abs(x-y) for x,y in zip(Q['center'],P['center']))/max(Q['r'],P['r'])
            for Q,P in combinations(cubes,2)]
    assert min(ratios)>1
    return str(min(ratios))


def input_account(cubes,source_boxes,masses,union_volume):
    m=sum(masses,F())
    integrals=[sum((mu*volume(intersect(Q['box'],S))/volume(S) for S,mu in zip(source_boxes,masses)),F()) for Q in cubes]
    assert all(I>volume(Q['box']) for I,Q in zip(integrals,cubes))
    # Each source box is wholly in exactly32 of the selected cubes.
    for S in source_boxes:
        hits=[volume(intersect(Q['box'],S))/volume(S) for Q in cubes]
        assert set(hits)<={F(0),F(1)} and sum(hits)==32
    summed=sum(integrals,F());sumvol=sum((volume(Q['box']) for Q in cubes),F())
    d_mass=32*m-summed;d_avg=summed-sumvol;d_geometry=sumvol-union_volume
    assert d_mass==0 and d_avg>0 and d_geometry>=0
    assert 32*m-union_volume==d_mass+d_avg+d_geometry
    return m,{'input_mass':str(m),'minimum_witness_average':str(min(I/volume(Q['box']) for I,Q in zip(integrals,cubes))),
        'sum_witness_integrals':str(summed),'sum_cube_volumes':str(sumvol),
        'source_outside_multiplicity32_mass':'0','multiplicity_defect_over_mass':str(d_mass/m),
        'average_slack_over_mass':str(d_avg/m),'union_volume_defect_over_mass':str(d_geometry/m),
        'selected_union_volume_over_mass':str(union_volume/m),
        'identity_warning':'Receiver in this finite identity is the selected union, not asserted to be a maximal strict level set.'}


def unequal():
    q=F(1,16);a=F(33,34);cubes=[]
    for sig in SIGNS:
        r=F(1) if sig==(1,)*5 else q
        c=tuple(a*r*s for s in sig)
        cubes.append({'r':r,'center':c,'box':box(c,r)})
    core,clique=clique_geometry(cubes,(F(0),)*5)
    assert core==box((F(0),)*5,F(1,544))
    endpoints=sorted({x for Q in cubes for interval in Q['box'] for x in interval})
    intervals=list(zip(endpoints,endpoints[1:]));assert len(intervals)==5
    U=F();N_integral=F();N32_volume=F()
    for cell in product(intervals,repeat=5):
        mid=tuple((lo+hi)/2 for lo,hi in cell);N=sum(inside(mid,Q['box']) for Q in cubes);v=volume(cell)
        if N: U+=v
        N_integral+=N*v
        if N==32: N32_volume+=v
    H=F(67,544)
    assert U==32+(2*H)**5-F(83,544)**5
    assert N_integral==sum((volume(Q['box']) for Q in cubes),F()) and N32_volume==volume(core)
    m,account=input_account(cubes,[core],[(1+ETA)*32],U)
    halo=(F(1,544)+1+ETA)**5/(1+ETA)
    assert halo<F(101,100)
    return {'id':'unequal_single32','cube_count':32,'minimum_normalized_pair_distance':family_check(cubes),
        'cliques':[clique],'arrangement_cells_checked':5**5,'account':account,
        'true_full_receiver_over_mass_upper':str(halo),
        'receiver_bound_method':'Core support halo, using (1+ETA)^(1/5)<1+ETA; no receiver sampling.'}


def multi(K,concentrated):
    a=F(4095,4096);b=F(5,2);L=16
    packets=[];cubes=[];cores=[];Ds=[];union_volume=F();base_masses=[];cliques=[]
    for j in range(K):
        r=F(1,L**j);common=(b*r,)+(F(0),)*4
        packet=[]
        for sig in SIGNS:
            c=tuple(z+a*r*s for z,s in zip(common,sig))
            packet.append({'r':r,'center':c,'box':box(c,r)})
        core,clique=clique_geometry(packet,common)
        U=box(common,(1+a)*r);D=box(common,(2+a)*r)
        packets.append((packet,U));cores.append(core);Ds.append(D);cliques.append(clique);cubes+=packet
        union_volume+=volume(U);base_masses.append((1+ETA)*(2*r)**5)
    for (_,U),(_,T) in combinations(packets,2): assert volume(intersect(U,T))==0
    assert all(volume(intersect(Ds[j],Ds[j+1]))==volume(Ds[j+1]) for j in range(K-1))
    assert all(inside((F(0),)*5,D) for D in Ds)
    masses=list(base_masses)
    if concentrated: masses[-1]+=100*sum(base_masses,F())
    m,account=input_account(cubes,cores,masses,union_volume)
    incidence=[[volume(intersect(D,S))/volume(S) for S in cores] for D in Ds]
    assert incidence==[[F(int(p>=j)) for p in range(K)] for j in range(K)]
    dilation_charge=sum((sum(mu*v for mu,v in zip(masses,row)) for row in incidence),F())
    row={'id':'three_packets_concentrated' if concentrated else 'two_packets_balanced',
        'cube_count':len(cubes),'packet_count':K,'minimum_normalized_pair_distance':family_check(cubes),
        'a':str(a),'b':str(b),'radius_ratio':L,'cliques':cliques,'account':account,
        'packet_unions_pairwise_disjoint':True,'twice_cube_packet_unions_nested':True,
        'dilated_packet_multiplicity_at_origin':K,'dilated_packet_source_incidence':[[int(x) for x in r] for r in incidence],
        'dilated_packet_source_charge_over_mass':str(dilation_charge/m),
        'receiver_warning':'The family is a valid finite Besicovitch witness family; it is not claimed to be a forced greedy cover of the entire maximal level set.'}
    if concentrated:
        assert dilation_charge/m>F(29,10)
    else:
        assert union_volume/m>F(3197,100)
        tau=F(1,100);Rbar=F(1003,1000);h=1-a
        assert Rbar**5>(1+ETA)/(1-tau)
        # Largest core halo + known weak32 for the small-radius source tail.
        halo=(h+Rbar)**5/(1+ETA)+F(32,1)/tau/F(L**5-1)
        assert halo<F(51,50)
        row['true_full_receiver_over_mass_upper']=str(halo)
        row['receiver_bound_method']='Split off largest core at threshold99/100; bound all remaining mass by geometric series and weak32 at threshold1/100.'
        row['weighted_dilation_charge_less_than128_over127']=dilation_charge/m<F(128,127)
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();root=Path(__file__).resolve().parent;started=time.monotonic()
    assert F(6,5)**5>2 and F(17,10)**5+16<F(127,4)<32
    assert 32-F(1,4)==F(127,4)
    rows=[unequal(),multi(2,False),multi(3,True)]
    delta_neighborhood=ETA*F(1,4096)/100
    assert (1+ETA)*(1-ETA/F(40))>1
    result={'status':'unequal_clique_lemma_and_source_localization_obstruction_checked_not_global_C5_improvement',
        'C5_improved':False,'local_unequal_radius_32_clique_bound':'127/4',
        'universal_global_defect_identity':'32m-lambda|E| = [32m-sum integral_Q f] + [sum integral_Q f-lambda sum |Q|] + lambda[sum |Q|-|E|]',
        'local_core_bound':'(17/10)^5','outside_core_mass_split':'1/4','rows':rows,
        'nonzero_receiver_neighborhood':{'relative_radius_for_multiscale_examples':str(delta_neighborhood),
            'average_lower_factor':'(1+ETA)*(1-ETA/40)>1','purpose':'Shows Lebesgue-positive receiver neighborhoods, not a full maximal level-set cover.'},
        'runtime':{'elapsed_seconds':time.monotonic()-started,'arithmetic':'Python Fraction exact rational geometry',
            'new_instances':3,'arrangement_cells':3125,'receiver_grid_samples':0,'MATLAB_used':False},
        'source_sha256':{p.name:sha(p) for p in [Path(__file__),root/'registration.json']},
        'remaining_gap':'No global selection/charging theorem controls the source cost of all32-cliques while retaining the real receivers. Whole-mass local bounds cannot be added.'}
    args.output_dir.mkdir(parents=True,exist_ok=True);out=args.output_dir/'results.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'output':str(out),'sha256':sha(out),'runtime':result['runtime'],
        'instances':[{'id':r['id'],'union_volume_over_mass':r['account']['selected_union_volume_over_mass'],
            'receiver_upper':r.get('true_full_receiver_over_mass_upper'),
            'dilation_charge_over_mass':r.get('dilated_packet_source_charge_over_mass'),
            'average_slack_over_mass':r['account']['average_slack_over_mass']} for r in rows]},indent=2))


if __name__=='__main__': main()
