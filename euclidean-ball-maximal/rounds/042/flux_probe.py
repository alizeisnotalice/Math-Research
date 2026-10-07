#!/usr/bin/env python3
"""Exact stage42 row/commutator ledger on retained actual stage41 certificates."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
sys.dont_write_bytecode=True
import argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'runtime/work/general_input_20261003'))
from stage30_transition_weight_probe import build_h
from stage24_clock_replacement_probe import observer
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--row-input',type=Path)
args=parser.parse_args()
data={'retained':{}}
for k in (10,12,14):
    r=F(2)**(2-k);N=2**(k-8)
    atoms=sorted([(-8*r*(i+F(1,2)),3*r/2) for i in range(N)]+[(F(9,32),F(27,64)),(F(10),F(71,128))])
    saved,_,radii,_=observer(atoms,F(1),'band',4)
    assert F(saved['entrance_a'])==F(1,8)
    cells=[dict(c,radius=str(radii[c['K']])) for c in saved['observer_cells']]
    data['retained']['stage41_micro_k'+str(k)]=dict(source=[dict(location=str(y),mass=str(w)) for y,w in atoms],alpha='1',observer_cells=cells)
if args.row_input:
    row_data=json.loads(args.row_input.read_text())
    cert=row_data['retained']['maximum_row_case'];summary=cert['summary']
    data['retained'][summary['label']]=dict(source=summary['source'],alpha=summary['alpha'],observer_cells=cert['observer_cells'])
out=[]
for label,cert in data['retained'].items():
    atoms=[(F(v['location']),F(v['mass'])) for v in cert['source']]
    alpha=F(cert['alpha']);cells=cert['observer_cells']
    hs=build_h(cells,alpha/8,'eligible');ks=sorted(hs);rows=[]
    for j in ks:
        hj=hs[j];tail=[k for k in ks if k>j];radj=F(4,2**j)/alpha;vj=2*radj
        xj=alpha*sum(F(c['hi'])-F(c['lo']) for c in cells if c['K']==j)
        row=sum(w*hj.value(y)*sum(hs[k].value(y) for k in tail) for y,w in atoms)
        bj=F(0);pos=F(0);neg=F(0)
        for cell in cells:
            k=cell['K']
            if k<=j:continue
            l=F(cell['lo']);r=F(cell['hi']);rad=F(cell['radius']);vk=2*rad
            cuts=sorted({l,r}|{z for z in hj.knots if l<z<r}|{z for y,w in atoms for z in (y-rad,y+rad) if l<z<r})
            for lo,hi in zip(cuts,cuts[1:]):
                mid=(lo+hi)/2
                mass=sum(w for y,w in atoms if abs(y-mid)<rad)
                bj+=(hi-lo)*hj.value(mid)*mass/vk
            for y,w in atoms:
                lo=max(l,y-rad);hi=min(r,y+rad)
                if lo>=hi:continue
                cuts=sorted({lo,hi}|{z for z in hj.knots if lo<z<hi})
                val=hj.value(y)
                for left,right in zip(cuts,cuts[1:]):
                    dl=val-hj.value(left);dr=val-hj.value(right)
                    nodes=[left,right]
                    if dl*dr<0:nodes.insert(1,left+(right-left)*dl/(dl-dr))
                    for u,v in zip(nodes,nodes[1:]):
                        area=(v-u)*(val-hj.value((u+v)/2))*w/vk
                        if area>=0:pos+=area
                        else:neg-=area
        assert row==bj+pos-neg
        assert 0<=bj<=xj
        occupation=F(0);conditional=F(0)
        for cell in cells:
            if cell['K']!=j:continue
            l=F(cell['lo']);r=F(cell['hi'])
            cuts=sorted({l,r}|{z for y,w in atoms for z in (y-radj,y+radj) if l<z<r})
            for lo,hi in zip(cuts,cuts[1:]):
                x=(lo+hi)/2
                within=[(y,w) for y,w in atoms if abs(y-x)<radj]
                q=sum(w for y,w in within)
                occ=sum(w*sum(hs[k].value(y) for k in tail) for y,w in within)
                assert alpha*vj/2<q<=alpha*vj
                occupation=max(occupation,occ/(alpha*vj));conditional=max(conditional,occ/q)
        rows.append(dict(j=j,Xj=str(xj),row=str(row),base=str(bj),positive_flow=str(pos),negative_flow=str(neg),
            row_over_Xj=str(row/xj),net_flow_over_Xj=str((pos-neg)/xj),positive_flow_over_Xj=str(pos/xj),
            negative_flow_over_Xj=str(neg/xj),max_actual_ball_occupation=str(occupation),max_conditional_occupation=str(conditional)))
    out.append(dict(label=label,rows=rows,max_row_over_Xj=str(max(F(r['row_over_Xj']) for r in rows)),
                    max_occupation=str(max(F(r['max_actual_ball_occupation']) for r in rows))))
result=dict(status='passed',arithmetic='Fraction exact',scope='regenerated micro family and supplied current-row maximum; original band and J/K',
            identity='row = base + positive_flow - negative_flow',cases=out,dimension_free_bound_proved=False)
path=args.output
path.write_text(json.dumps(result,indent=2)+'\n')
for c in out:
    print(c['label'],'max row',float(F(c['max_row_over_Xj'])),'max occupation',float(F(c['max_occupation'])))
print(path)
