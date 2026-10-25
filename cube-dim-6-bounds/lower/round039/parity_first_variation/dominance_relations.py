#!/usr/bin/env python3
"""Twelve bounded exact relation certificates; no union or numerical roots."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,os,time

HERE=Path(__file__).resolve().parent
IP=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/lower/round008/spacing_chamber/chamber_interface.json')

def main():
    started=time.process_time()
    inventory=json.loads((HERE/'facet_inventory_results.json').read_text())
    labels={tuple(x) for x in json.loads(IP.read_text())['box_lengths']}
    assert len(labels)==127
    records=[]
    def check(name,condition,meaning):
        assert condition,name
        records.append({'name':name,'status':'PASS','meaning':meaning})
    check('01',(Q(6,5))**6>2,'u<6/5, hence a_tail<11/10')
    check('02',(Q(4,3))**3>2,'x=u^2<4/3')
    check('03',(Q(3,2))**3>3,'y=v^2<3/2, hence a2<7/6')
    check('04',Q(6,5)**2<Q(3,2),'r216=sqrt(6)/2>6/5')
    check('05',Q(7,6)<Q(6,5),'r216>a2')
    check('06',Q(11,10)<Q(6,5),'r216>a_tail')
    check('07',216>108,'r216>r108')
    check('08',324>216,'r324>r216')
    check('09',216>162,'r324>r162 by transitivity')
    fixed={}
    for group in inventory['same_orientation_noncell_velocity_split_groups']:
        assert group['orientation']=='upper'
        for item in group['occurrences']:
            if item['l1']==2:
                fixed[item['label_index']]=item
    f108=[x for x in fixed.values() if x['mass']==108]
    f162=[x for x in fixed.values() if x['mass']==162]
    check('10',len(f108)==20 and all(sorted(x['lengths'][1:])==[1,2,3,3,3] for x in f108),
          'all20 fixed108 labels have other lengths perm[1,2,3,3,3]')
    check('11',len(f162)==5 and len(fixed)==25 and
          all(sorted(x['lengths'][1:])==[1,3,3,3,3] for x in f162),
          'all5 fixed162 labels have other lengths perm[1,3,3,3,3]; no other fixed mixed labels')
    dominators={}
    pairs=[]
    for item in fixed.values():
        ell=tuple(item['lengths'])
        c=list(ell);k=c.index(1);assert k>0;c[k]=2;c=tuple(c)
        dominators[c]=True
        pairs.append({'dominated_label_index':item['label_index'],'B':ell,'C':c,
                      'changed_axis_one_based':k+1,'mass_B':item['mass'],
                      'mass_C':2*item['mass']})
    check('12',len(dominators)==15 and set(dominators)<=labels,
          'all15 exact dominating C labels already belong to inherited127 family')
    assert len(records)==12
    out={'status':'PASS','pid':os.getpid(),'new_exact_relation_count':12,
         'relations':records,'dominated_fixed_labels':25,'dominating_existing_labels':15,
         'dominance_pairs':sorted(pairs,key=lambda x:x['dominated_label_index']),
         'conclusion':'all noncell same-orientation velocity-split facets have fixed-phase contributors hidden by existing fixed boxes; no positive |epsilon| first variation',
         'full_union_calls':0,'new_epsilon_candidates':0,'numerical_root_evaluations':0,
         'CPU_seconds':time.process_time()-started,
         'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (IP,HERE/'facet_inventory_results.json')},
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'actual_process_status':'finished; no more calculations'}
    (HERE/'dominance_relations_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('relations','dominance_pairs')},indent=2))

if __name__=='__main__':main()
