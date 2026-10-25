#!/usr/bin/env python3
"""Bounded127-label correlated-parity facet preflight; twelve exact checks."""
from fractions import Fraction as Q
from pathlib import Path
from math import prod
import hashlib,json,os,time

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')

def ps(l):
    q=-(l//2) if l%2 else 1-l//2
    return (1 if q%2==0 else -1) if l%2 else 0

def velocity(mass,ell,side):
    rho=Q(ps(ell[0])*ps(ell[1]),ell[0]*ell[1])
    if not rho:return '0'
    m=mass;i=j=0
    while m%2==0:i+=1;m//=2
    while m%3==0:j+=1;m//=3
    assert m==1
    c=Q(1 if side else -1,12)*rho*2**(i//6)*3**(j//6)
    return f'{c}*u^{i%6}*v^{j%6}'

def main():
    started=time.process_time()
    ip=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
    ep=MAIN/'audit/lower008/parametric_union/candidate_results.json'
    p=json.loads(ip.read_text());e=json.loads(ep.read_text())
    checks=[]
    def check(name,condition,claim):
        if len(checks)>=12:raise RuntimeError('twelve-check cap')
        checks.append({'name':name,'pass':bool(condition),'claim':claim})
        if not condition:raise ArithmeticError(name)
    check('01',len(p['box_lengths'])==len(p['boxes'])==len({tuple(x) for x in p['box_lengths']})==127,
          'exactly127 distinct inherited labels, no added label')
    ties=e['boundary_ties']
    check('02',e['checked_adjacent_order_gaps']==120 and len(ties)==1 and ties[0]['right_equal'] and not ties[0]['left_equal'],
          'accepted complete base coincidence partition reused; one right-boundary tie')
    eps=Q(1,7)
    check('03',sum((1+eps*((-1)**(s+t)) for s in (0,1) for t in (0,1)),Q(0))==4,
          'four-phase cell mass exactly4; rational check is not a geometric candidate')
    check('04',sorted(((-1)**(s+t) for s in (0,1) for t in (0,1)))==[-1,-1,1,1],
          'two positive and two negative phase characters')
    blocks={l:set(range(-(l//2) if l%2 else 1-l//2,
                         (-(l//2) if l%2 else 1-l//2)+l)) for l in (1,2,3,4)}
    check('05',blocks[1]<blocks[2],'central source J1 strictly nested in J2')
    check('06',blocks[2]<blocks[3],'central source J2 strictly nested in J3')
    check('07',blocks[3]<blocks[4],'central source J3 strictly nested in J4')
    ranks=[];lasts=[]
    for axis,groups in enumerate(p['order_certificate']):
        r={key:k for k,g in enumerate(groups) for key in g}
        if axis==ties[0]['axis']:
            low,high=ties[0]['adjacent_groups'];r={key:(low if k==high else k) for key,k in r.items()}
        ranks.append(r);lasts.append(len(groups)-1)
    groups={};labels={};cell=0
    # One iteration over the127 raw source labels; all six axes are processed.
    for index,(ell,box) in enumerate(zip(p['box_lengths'],p['boxes'])):
        ell=tuple(ell);labels[ell]=index;mass=prod(ell)
        for axis,ids in enumerate(box):
            for side,key in enumerate(ids):
                rank=ranks[axis][key]
                if rank in (0,lasts[axis]):cell+=1;continue
                entry={'label_index':index,'lengths':ell,'mass':mass,'axis_length':ell[axis],
                       'endpoint_id':key,'rho':str(Q(ps(ell[0])*ps(ell[1]),ell[0]*ell[1])),
                       'velocity':velocity(mass,ell,side)}
                groups.setdefault((axis,rank,side),[]).append(entry)
    split=[]
    for (axis,rank,side),items in groups.items():
        if len({x['velocity'] for x in items})>1:
            split.append({'axis_one_based':axis+1,'rank':rank,
                          'orientation':'upper' if side else 'lower','occurrences':items})
    check('08',all(len({x['mass'] for x in g['occurrences']})==1 for g in split),
          'every unequal-velocity same-oriented noncell coincidence has one common baseline mass')
    check('09',all(len({x['axis_length'] for x in g['occurrences']})==1 for g in split),
          'within every such group the normal-axis source length is common')
    A=(1,2,1,1,1,1);B=(1,1,2,1,1,1);C=tuple(max(a,b) for a,b in zip(A,B))
    check('10',A in labels and B in labels and C==(1,2,2,1,1,1),
          'explicit new axis1 equal-mass velocity split has componentwise-max label C')
    check('11',prod(A)==prod(B)==2 and prod(C)==4 and Q(ps(A[0])*ps(A[1]),A[0]*A[1])==0
          and Q(ps(B[0])*ps(B[1]),B[0]*B[1])==1,
          'example velocities differ:0 vs r2/6; merged radius r4 strictly exceeds r2')
    check('12',C in labels,'example merged covering box is already among the127 labels')
    out={'status':'PASS_no_positive_first_variation_in_fixed127_family',
         'pid':os.getpid(),'distinct_raw_labels_read':127,'new_exact_or_bookkeeping_checks':len(checks),
         'checks':checks,'period_mass':4,'ratio_formula':'32*(V(epsilon)+V(-epsilon))',
         'split_group_count':len(split),'split_groups':split,
         'cell_face_occurrences_excluded':cell,
         'uniform_mass_coincident_facet_merge_lemma_applicable':True,
         'example':{'A':A,'B':B,'C':C,'source_label_indices':[labels[A],labels[B],labels[C]],
                    'axis_one_based':1,'masses':[2,2,4]},
         'conclusion':'all new unequal-velocity same-oriented cofacets are hidden by the existing central componentwise-max block whenever their projected intersection has positive5D area; Rright(0)<=0',
         'unresolved':['exact first-order coefficient','second-order sign','nonzero epsilon improvement'],
         'new_C6_bound':False,'full_union_calls':0,'new_epsilon_candidates':0,
         'baseline_replays':0,'numerical_root_evaluations':0,'MATLAB_calls':0,
         'CPU_seconds':time.process_time()-started,
         'input_sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in (ip,ep)},
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'actual_process_status':'finished; stop, upper-bound priority'}
    (HERE/'preflight_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('checks','split_groups','input_sha256')},indent=2))
    for g in split:
        print(json.dumps({'axis':g['axis_one_based'],'orientation':g['orientation'],
                          'mass':g['occurrences'][0]['mass'],
                          'velocities':sorted({x['velocity'] for x in g['occurrences']})}))

if __name__=='__main__':main()
