#!/usr/bin/env python3
"""One127-label pass; inherited base ranks only, no union/exposure algorithm."""
from fractions import Fraction as Q
from pathlib import Path
from math import prod
import hashlib,json,os,time

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')

def velocity(mass,l1,side):
    if l1%2==0:
        return '0'
    m=mass;i=j=0
    while m%2==0:i+=1;m//=2
    while m%3==0:j+=1;m//=3
    assert m==1
    sigma=1 if l1==1 else -1
    c=Q(sigma*(1 if side else -1),12*l1)*2**(i//6)*3**(j//6)
    return f'{c}*u^{i%6}*v^{j%6}'

def main():
    started=time.process_time()
    ip=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
    ep=MAIN/'audit/lower008/parametric_union/candidate_results.json'
    packet=json.loads(ip.read_text()); accepted=json.loads(ep.read_text())
    assert len(packet['box_lengths'])==len(packet['boxes'])==127
    assert accepted['input_boxes']==127 and accepted['checked_adjacent_order_gaps']==120
    ties=accepted['boundary_ties']
    assert len(ties)==1 and ties[0]['right_equal'] and not ties[0]['left_equal']
    occurrences={};cell_faces=0
    for axis,groups in enumerate(packet['order_certificate']):
        ranks={key:rank for rank,group in enumerate(groups) for key in group}
        # Exactly the accepted t=1 boundary tie, without a new field comparison.
        if axis==ties[0]['axis']:
            left,right=ties[0]['adjacent_groups']
            ranks={key:(left if rank==right else rank) for key,rank in ranks.items()}
        last=len(groups)-1
        for index,(ell,box) in enumerate(zip(packet['box_lengths'],packet['boxes'])):
            mass=prod(ell)
            for side,key in enumerate(box[axis]):
                rank=ranks[key]
                if rank in (0,last):
                    cell_faces+=1
                    continue
                entry={'label_index':index,'lengths':ell,'mass':mass,'l1':ell[0],
                       'endpoint_id':key,'velocity':velocity(mass,ell[0],side)}
                occurrences.setdefault((axis,rank,side),[]).append(entry)
    split=[]
    for (axis,rank,side),items in occurrences.items():
        velocities={item['velocity'] for item in items}
        if len(velocities)>1:
            split.append({'axis_one_based':axis+1,'inherited_rank':rank,
                          'orientation':'upper' if side else 'lower',
                          'velocity_count':len(velocities),'occurrences':items})
    out={'status':'inventory_complete_exposure_not_decided','pid':os.getpid(),
         'raw_label_count':127,'label_passes':1,'base_order_source':'accepted E008 exact rank/gap certificate',
         'base_boundary_ties_reused':ties,'same_orientation_noncell_velocity_split_groups':split,
         'split_group_count':len(split),'cell_face_occurrences_excluded':cell_faces,
         'new_union_calls':0,'new_endpoint_field_sign_tests':0,'new_epsilon_candidates':0,
         'CPU_seconds':time.process_time()-started,
         'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ip,ep)},
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'actual_process_status':'finished; further exposure work requires main authorization'}
    (HERE/'facet_inventory_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='same_orientation_noncell_velocity_split_groups'},indent=2))
    for group in split:
        print(json.dumps({'axis':group['axis_one_based'],'orientation':group['orientation'],
                          'rank':group['inherited_rank'],'n':len(group['occurrences']),
                          'velocities':sorted({x['velocity'] for x in group['occurrences']}),
                          'l1_classes':sorted({x['l1'] for x in group['occurrences']})}))

if __name__=='__main__':main()
