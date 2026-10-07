#!/usr/bin/env python3
"""Exact common-threshold record partitions, fixed radii and full eligible band."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import importlib.util,json,sys,argparse,time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]

def load(root,out):
    spec=importlib.util.spec_from_file_location('record44',root/'rounds/044/tail_probe.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.ROOT=root
    return m.archived43(out)

def intervals(trace,lo,hi):
    record=F(0);result={}
    for k,g in enumerate(trace,1):
        a=max(lo,record);b=min(hi,g)
        if a<b:result[k]=(a,b)
        record=max(record,g)
    assert sum(b-a for a,b in result.values())==hi-lo
    return result

def pair_audit(x,z,lo,hi):
    ix=intervals(x,lo,hi);iz=intervals(z,lo,hi)
    joint={};ordered=F(0)
    for j,(a,b) in ix.items():
        for k,(c,d) in iz.items():
            width=max(F(0),min(b,d)-max(a,c))
            joint[j,k]=width
            if j<k:
                assert width<=max(F(0),x[j-1]-z[j-1]);ordered+=width
    assert ordered<=hi-lo and sum(joint.values())==hi-lo
    # Independent first-exit reconstruction on every scalar threshold interval.
    cuts=sorted({lo,hi}|{g for g in x+z if lo<g<hi});direct=defaultdict(F)
    for a,b in zip(cuts,cuts[1:]):
        beta=(a+b)/2;j=next(k for k,g in enumerate(x,1) if g>beta);k=next(k for k,g in enumerate(z,1) if g>beta)
        direct[j,k]+=b-a
    assert all(direct[key]==width for key,width in joint.items()) and all(key in joint for key in direct)
    return len(joint)

def cases(old):
    base={'separated_equal_4','separated_binary','two_near_half'}
    for b,label,atoms,alpha,p in old.cases():
        if label in base:yield 0,label,atoms,alpha
    for k in (10,12):
        r=F(2)**(2-k);N=2**(k-8)
        yield 1,f'micro{k}',sorted([(-8*r*(i+F(1,2)),3*r/2) for i in range(N)]+[(F(9,32),F(27,64)),(F(10),F(71,128))]),F(1)
    for L in (6,12):yield 1,f'chain{L}',old.old.old.chain(L,den=4,position=F(9,8),weight=F(27,16)),F(1)
    yield 1,'enhanced12',old.old.old.chain(12,den=16,position=F(9,8),weight=F(7,4)),F(1)
    yield 1,'minimal_positive_prefix',[(F(0),F(13,256)),(F(1,16),F(3,64)),(F(3,32),F(9,128)),(F(10),F(213,256))],F(1)
    e=F(1,2048);p=F(63,2048)
    yield 1,'old_density_gap',[(F(-1,8),2*e),(F(1,10),F(1,4)-p-e),(F(511,2048),p),(F(10),F(3,4)-e)],F(1)
    labels={'random_two_clouds_s420701_L6','random_two_clouds_s420720_L10','random_two_clouds_s420739_L14','asymmetric_tree_d3_s420924','asymmetric_tree_d4_s420947','asymmetric_tree_d5_s421154','big_contact_L6_b8_s-1_side1','band_edge_sliver_L6_b10_side1'}
    for b,label,atoms,alpha,p in old.cases():
        if label in labels:yield 2,label,atoms,alpha
    yield 2,'deep_chain24',old.old.old.chain(24,den=4,position=F(9,8),weight=F(27,16)),F(2)

def evaluate(old,batch,label,atoms,alpha):
    saved,_,radii,_=old.observer(atoms,alpha,'band',4);D=saved['D'];lo=alpha/4;hi=alpha/2
    X=alpha*F(saved['eligible_volume']);assert X>0
    edges={y+sign*radii[k] for y,w in atoms for k in range(1,D+1) for sign in (-1,1)}
    cells=[];cuts={lo,hi}
    for c in saved['observer_cells']:
        a=F(c['lo']);b=F(c['hi']);split=sorted({a,b}|{x for x in edges if a<x<b})
        for l,r in zip(split,split[1:]):
            z=(l+r)/2
            trace=[sum((w for y,w in atoms if abs(y-z)<radii[k]),F(0))/(2*radii[k]) for k in range(1,D+1)]
            assert trace[0]<=alpha/8 and max(trace)>hi
            assert all(trace[k]<=2*trace[k-1] for k in range(1,D))
            cells.append((l,r,c['J'],trace));cuts.update(g for g in trace if lo<g<hi)
    cuts=sorted(cuts);averageQ=averageM=averageF=averageB=averageD=F(0);steps=[]
    def integral(fun,l,r):
        cs=sorted({l,r}|{z for z in fun.knots if l<z<r})
        return sum(((b-a)*fun.value((a+b)/2) for a,b in zip(cs,cs[1:])),F(0))
    for a,b in zip(cuts,cuts[1:]):
        beta=(a+b)/2;assigned=[]
        for l,r,J,trace in cells:
            K=next(k for k,g in enumerate(trace,1) if g>beta);g=trace[K-1]
            assert J<=K and beta<g<=2*beta<=alpha
            assigned.append(dict(lo=str(l),hi=str(r),J=J,K=K,radius=str(radii[K])))
        hs=old.build_h(assigned,alpha/8,'eligible');A=old.add(hs.values())
        M=sum((w*A.value(y) for y,w in atoms),F(0));Q=sum((w*A.value(y)**2 for y,w in atoms),F(0))
        diagonal=sum((w*sum((h.value(y)**2 for h in hs.values()),F(0)) for y,w in atoms),F(0))
        prefixes={k:old.add(hs[j] for j in hs if j<k) for k in hs}
        B=F(0)
        for row,cell in zip(assigned,cells):
            l,r,J,trace=cell;k=row['K'];B+=trace[k-1]*integral(prefixes[k],l,r)
        flow=(Q-diagonal)/2-B
        assert beta*X/alpha<=M<=2*beta*X/alpha and M*M<=Q
        assert diagonal<=M and B<=2*beta*X/alpha
        averageM+=(b-a)*M/(hi-lo);averageQ+=(b-a)*Q/(hi-lo)
        averageF+=(b-a)*flow/(hi-lo);averageB+=(b-a)*B/(hi-lo);averageD+=(b-a)*diagonal/(hi-lo)
        steps.append(dict(lo=str(a),hi=str(b),M=str(M),Q=str(Q),B=str(B),F=str(flow)))
    assert F(7,48)*X*X<=averageQ
    assert averageQ==averageD+2*averageB+2*averageF and averageQ<=F(9,4)*X+2*averageF
    selected=[cells[i][3] for i in sorted({i*(len(cells)-1)//11 for i in range(12)})]
    pairs=sum(pair_audit(x,z,lo,hi) for x in selected for z in selected)
    pair_integral=None
    if label in ('two_near_half','minimal_positive_prefix','big_contact_L6_b8_s-1_side1'):
        raw=F(0)
        for xa,xb,J,tx in cells:
            ix=intervals(tx,lo,hi)
            for za,zb,J,tz in cells:
                iz=intervals(tz,lo,hi)
                for j,(la,lb) in ix.items():
                    for k,(ra,rb) in iz.items():
                        width=max(F(0),min(lb,rb)-max(la,ra))
                        if j>=k or not width:continue
                        R=radii[j];r=radii[k]
                        length=lambda a,b,y,s:max(F(0),min(b,y+s)-max(a,y-s))
                        first=sum((w*length(xa,xb,y,R)*length(za,zb,y,r)/(4*R*r) for y,w in atoms),F(0))
                        # Rectangle integral of the coarse ball kernel, by exact hinge primitive.
                        square=lambda x:max(F(0),x)**2
                        area=F(0)
                        for shift,sign in ((R-xa,1),(R-xb,-1),(-R-xa,-1),(-R-xb,1)):
                            area+=sign*(square(zb+shift)-square(za+shift))/2
                        assert area>=0
                        second=tz[k-1]*area/(2*R)
                        raw+=width*(first-second)
        pair_integral=raw/(hi-lo);assert pair_integral==averageF
    return dict(batch=batch,label=label,N=len(atoms),D=D,alpha=str(alpha),X=str(X),observer_cell_count=len(cells),threshold_interval_count=len(steps),mean_M=str(averageM),mean_Q=str(averageQ),mean_Q_over_X=str(averageQ/X),mean_F=str(averageF),mean_B=str(averageB),mean_diagonal=str(averageD),direct_pair_F=str(pair_integral) if pair_integral is not None else None,pair_entries_checked=pairs,steps=steps)

def gap_test():
    rows=[]
    for power in (10,12,16,24,32,48):
        e=F(1,2**power);p=F(63,2048)
        atoms=[(F(-1,8),2*e),(F(1,10),F(1,4)-p-e),(F(511,2048),p),(F(10),F(3,4)-e)]
        def trace(z):return [sum((w for y,w in atoms if abs(y-z)<F(2)**(2-k)),F(0))/(2*F(2)**(2-k)) for k in range(1,9)]
        tx=trace(F(0));tz=trace(F(271,1024));ix=intervals(tx,F(1,4),F(1,2));iz=intervals(tz,F(1,4),F(1,2))
        a,b=ix[4];c,d=iz[8];width=max(F(0),min(b,d)-max(a,c))
        assert width==2*e
        rows.append(dict(epsilon=str(e),width_4_8=str(width),probability_4_8=str(4*width),density_gap=str(tx[3]-tz[3])))
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);a=p.parse_args();old=load(a.repo_root,a.output)
    result=dict(status='running',scope='n1 actual fixed observer, fixed radii, one common uniform threshold',cases=[])
    for batch,label,atoms,alpha in cases(old):
        r=evaluate(old,batch,label,atoms,alpha);result['cases'].append(r);a.output.write_text(json.dumps(result));print(label,r['threshold_interval_count'],flush=True)
    result.update(status='passed',gap_tests=gap_test(),batch_counts=[sum(r['batch']==b for r in result['cases']) for b in range(3)])
    a.output.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
