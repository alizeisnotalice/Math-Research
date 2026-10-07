#!/usr/bin/env python3
"""Fraction-exact original-band global source/observer ledgers, stage43.
Uses archived observer and piecewise linear kernels; writes only the requested output.
"""
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import sys, json, time, argparse, hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--enhanced-only',action='store_true');args=parser.parse_args();OUT=args.output
sys.path.insert(0,str(ROOT/'rounds/042'))
# row_probe imports historical constructors with their shared output parser only.
sys.argv=[sys.argv[0],'--output',str(OUT)]
import row_probe as old
from stage24_clock_replacement_probe import observer
from stage30_transition_weight_probe import build_h,Linear

def add(fs):
    ev=defaultdict(F)
    for f in fs:
        for z,d in f.events.items():ev[z]+=d
    return Linear(ev)

def integrate_flow(fun,y,l,r,w,v):
    """Exact positive/negative linear drift on one atom's ball and one cell."""
    p=n=F(0);hy=fun.value(y)
    nodes=sorted({l,r}|{z for z in fun.knots if l<z<r})
    for lo,hi in zip(nodes,nodes[1:]):
        dl=hy-fun.value(lo);dr=hy-fun.value(hi);cs=[lo,hi]
        if dl*dr<0:cs.insert(1,lo+(hi-lo)*dl/(dl-dr))
        for u,t in zip(cs,cs[1:]):
            a=(t-u)*(hy-fun.value((u+t)/2))*w/v
            if a>=0:p+=a
            else:n-=a
    return p,n

def evaluate(batch,label,atoms,alpha,params):
    started=time.monotonic();saved,_,radii,_=observer(atoms,alpha,'band',4)
    assert F(saved['entrance_a'])==alpha/8 and F(saved['terminal_b'])==alpha/2
    cells=[dict(c,radius=str(radii[c['K']])) for c in saved['observer_cells']]
    hs=build_h(cells,alpha/8,'eligible');ks=sorted(hs);A=add(hs.values())
    all_h_knots={z for h in hs.values() for z in h.knots}
    prefixes={k:add(hs[j] for j in ks if j<k) for k in ks}
    X=alpha*F(saved['eligible_volume']);assert X>0
    vectors=[(y,w,{k:h.value(y) for k,h in hs.items()}) for y,w in atoms]
    M=sum((w*sum(v.values()) for y,w,v in vectors),F(0))
    Q=sum((w*sum(v.values())**2 for y,w,v in vectors),F(0))
    diag=sum((w*sum(t*t for t in v.values()) for y,w,v in vectors),F(0))
    rows={j:sum((w*v[j]*sum(v[k] for k in ks if k>j) for y,w,v in vectors),F(0)) for j in ks}
    row=sum(rows.values(),F(0));B=T=Dobs=C=SV=AV=F(0);pos=neg=fullpos=fullneg=F(0)
    observer_square=F(0)
    lebesgue_square=sum(((hi-lo)*(A.value(lo)**2+A.value(lo)*A.value(hi)+A.value(hi)**2)/3 for lo,hi in zip(A.knots,A.knots[1:])),F(0))
    bases=defaultdict(F);xby=defaultdict(F);Smax=F(0);Aobsmax=F(0);witness=None
    for c in cells:
        k=c['K'];l=F(c['lo']);r=F(c['hi']);rad=radii[k];vk=2*rad;sk=prefixes[k];hk=hs[k]
        xby[k]+=alpha*(r-l)
        cuts=sorted({l,r}|{z for z in all_h_knots if l<z<r}|{z for y,w in atoms for z in (y-rad,y+rad) if l<z<r})
        for lo,hi in zip(cuts,cuts[1:]):
            z=(lo+hi)/2;g=sum((w for y,w in atoms if abs(y-z)<rad),F(0))/vk
            assert alpha/2<g<=alpha
            sv=sk.value(z);av=A.value(z);hv=hk.value(z);width=hi-lo
            observer_square+=width*(A.value(lo)**2+A.value(lo)*A.value(hi)+A.value(hi)**2)*g/3
            B+=width*sv*g;T+=width*(av-sv-hv)*g;Dobs+=width*hv*g;C+=width*av*g;SV+=width*sv;AV+=width*av
            for j in ks:
                if j<k:bases[j]+=width*hs[j].value(z)*g
            for z0 in (lo,hi):
                sv0=sk.value(z0);av0=A.value(z0)
                if sv0>Smax:Smax=sv0;witness=dict(z=str(z0),K=k,S=str(sv0),cell=c)
                Aobsmax=max(Aobsmax,av0)
        for y,w in atoms:
            lo=max(l,y-rad);hi=min(r,y+rad)
            if lo>=hi:continue
            p,n=integrate_flow(sk,y,lo,hi,w,vk);pos+=p;neg+=n
            p,n=integrate_flow(A,y,lo,hi,w,vk);fullpos+=p;fullneg+=n
    net=row-B;full=Q-C
    assert net==pos-neg and full==fullpos-fullneg
    assert Q==diag+2*B+2*net and C==B+T+Dobs
    assert sum(bases.values())==B and sum(xby.values())==X
    assert 0<=B<=alpha*SV<=X and 0<=C<=alpha*AV<=X
    assert X/2<=M<=X and diag<=M
    rownets=[rows[j]-bases[j] for j in ks]
    rowpositive=sum((max(F(0),v) for v in rownets),F(0));rownegative=sum((max(F(0),-v) for v in rownets),F(0))
    assert net==rowpositive-rownegative
    sourceAmax=max((sum(v.values()) for y,w,v in vectors),default=F(0))
    # Independent direct-overlap kernel reconstruction at every source.
    for y,w,v in vectors:
        for k in ks:
            direct=sum((max(F(0),min(F(c['hi']),y+radii[k])-max(F(c['lo']),y-radii[k]))/(2*radii[k]) for c in cells if c['K']==k),F(0))
            assert direct==v[k]
    assert observer_square<=alpha*lebesgue_square
    vals=dict(lebesgue_square=alpha*lebesgue_square,observer_square=observer_square,X=X,M=M,Q=Q,diag=diag,row=row,B=B,F=net,coarse_positive=pos,coarse_negative=neg,
              aggregate_cancellation=min(pos,neg),row_positive_net=rowpositive,row_negative_net=rownegative,
              fine_observer=T,fine_drift=row-T,diag_observer=Dobs,diag_drift=diag-Dobs,
              full_observer=C,full_drift=full,full_positive=fullpos,full_negative=fullneg,
              S_integral=SV,A_observer_integral=AV,max_observer_S=Smax,max_observer_A=Aobsmax,max_source_A=sourceAmax)
    ratios={key+'_over_X':str(v/X) for key,v in vals.items() if key not in ('X','max_observer_S','max_observer_A','max_source_A','S_integral','A_observer_integral')}
    ratios.update(S_average_over_E=str(SV/(X/alpha)),A_average_over_E=str(AV/(X/alpha)),coarse_cancellation_fraction=str(2*min(pos,neg)/(pos+neg)) if pos+neg else '0')
    return dict(batch=batch,label=label,parameters=params,N=len(atoms),D=saved['D'],K_count=len(ks),cell_count=len(cells),alpha=str(alpha),
                **{key:str(v) for key,v in vals.items()},**ratios,
                source=[dict(location=str(y),mass=str(w),A=str(sum(v.values()))) for y,w,v in vectors],
                rows=[dict(j=j,Xj=str(xby[j]),row=str(rows[j]),B=str(bases[j]),F=str(rows[j]-bases[j])) for j in ks],
                S_supremum_witness=witness,observer_cells=cells,original_intervals=saved['original_intervals'],
                J1_volume=saved['J1_volume'],original_volume=saved['original_volume'],
                assertions='original band/J/K; source direct overlaps; Q=diag+2B+2F; signed flow; B<=alpha integral S<=X; complete observer<=X',
                exact_error='0',elapsed_seconds=time.monotonic()-started)

def cases():
    # Baseline: separated probabilities and translated versions, no coarse cross.
    for n in (4,8):
        yield 0,f'separated_equal_{n}',[(F(5*i),F(1,n)) for i in range(n)],F(1),{}
    yield 0,'separated_binary',[(F(0),F(1,4)),(F(5),F(1,8)),(F(10),F(1,8)),(F(15),F(1,2))],F(1),{}
    yield 0,'two_near_half',[(F(-1,8),F(3,8)),(F(9,32),F(7,16)),(F(10),F(3,16))],F(1),{}
    # Retained old positive row witnesses, full global recomputation.
    for k in (10,12,14):
        r=F(2)**(2-k);N=2**(k-8)
        atoms=sorted([(-8*r*(i+F(1,2)),3*r/2) for i in range(N)]+[(F(9,32),F(27,64)),(F(10),F(71,128))])
        yield 1,f'old41_micro_k{k}',atoms,F(1),dict(k=k)
    for L in (12,24,48):
        yield 1,f'old41_chain_L{L}',old.old.chain(L,den=4,position=F(9,8),weight=F(27,16)),F(1),dict(L=L)
    yield 1,'old41_nested_d5',old.old.finish(old.old.tree(5)),F(1),dict(depth=5)
    # More complicated original cells: contact, unequal clouds, and asymmetric trees.
    wanted={'big_contact_L6_b8_s-1_side1','band_edge_sliver_L6_b10_side1',
        'random_two_clouds_s420701_L6','random_two_clouds_s420720_L10','random_two_clouds_s420739_L14',
        'asymmetric_tree_d2_s420901','asymmetric_tree_d3_s420924','asymmetric_tree_d4_s420947','asymmetric_tree_d5_s421154'}
    for _,label,atoms,alpha,params in old.cases():
        if label in wanted:yield 2,label,atoms,alpha,params

def finalize(out):
    for c in out['cases']:
        c['source_occupancy_tails']={}
        for level in (1,2,4,8):
            sources=[v for v in c['source'] if F(v['A'])>=level]
            mass=sum((F(v['mass']) for v in sources),F(0))
            q=sum((F(v['mass'])*F(v['A'])**2 for v in sources),F(0))
            c['source_occupancy_tails'][str(level)]=dict(mass=str(mass),Q_contribution=str(q),Q_contribution_over_X=str(q/F(c['X'])))
    out['independent_band_clock_verifications']=[]
    for label in ('old41_micro_k12','asymmetric_tree_d5_s421154'):
        c=next(c for c in out['cases'] if c['label']==label)
        summary=dict(c,rows=[dict(r,Rj=r['row']) for r in c['rows']],max_row=None)
        check=old.independent(dict(summary=summary))
        out['independent_band_clock_verifications'].append(dict(label=label,**check))
    out['sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return out

def enhanced_windows():
    results=[]
    for L in (4,12,24):
        atoms=old.old.chain(L,den=16,position=F(9,8),weight=F(7,4))
        cert=evaluate(3,'enhanced_chain_L'+str(L),atoms,F(1),dict(L=L))
        cells=cert['observer_cells'];hs=build_h(cells,F(1,8),'eligible');ks=sorted(hs)
        prefixes={k:add(hs[j] for j in ks if j<k) for k in ks}
        windows=[]
        for i in range(L):
            ri=F(1,4*8**i);lo=ri/4;hi=11*ri/16;k=4+3*i
            parts=[(max(lo,F(c['lo'])),min(hi,F(c['hi'])),c) for c in cells if max(lo,F(c['lo']))<min(hi,F(c['hi']))]
            assert sum((b-a for a,b,c in parts),F(0))==hi-lo
            assert all(c['J']==2 and c['K']==k for a,b,c in parts)
            windows.append(dict(i=i,K=k,lo=str(lo),hi=str(hi),full_window_actual_J2_K_verified=True))
        rl=F(1,4*8**(L-1));lo=rl/4;hi=11*rl/16;k=4+3*(L-1);sk=prefixes[k];bound=F(7*(L-1),32)
        cuts=sorted({lo,hi}|{z for z in sk.knots if lo<z<hi})
        minimum=min(sk.value(z) for z in cuts);assert minimum>=bound
        # Exact full observer heavy tail: threshold cuts are rational roots of S.
        tailvolume=tailbase=F(0)
        for c in cells:
            l=F(c['lo']);r=F(c['hi']);rad=F(c['radius']);s=prefixes[c['K']]
            cuts=sorted({l,r}|{z for z in s.knots if l<z<r})
            g=sum((w for y,w in atoms if abs(y-(l+r)/2)<rad),F(0))/(2*rad)
            for a,b in zip(cuts,cuts[1:]):
                dl=s.value(a)-bound;dr=s.value(b)-bound;nodes=[a,b]
                if dl*dr<0:nodes.insert(1,a+(b-a)*dl/(dl-dr))
                for u,v in zip(nodes,nodes[1:]):
                    mid=(u+v)/2
                    if s.value(mid)>=bound:tailvolume+=v-u;tailbase+=(v-u)*s.value(mid)*g
        X=F(cert['X'])
        results.append(dict(L=L,source=cert['source'],windows=windows,observer_S_lower_bound=str(bound),
          min_S_on_last_window=str(minimum),last_window_volume=str(hi-lo),X=cert['X'],B_over_X=cert['B_over_X'],F_over_X=cert['F_over_X'],Q_over_X=cert['Q_over_X'],
          max_observer_S=cert['max_observer_S'],heavy_tail_threshold=str(bound),heavy_tail_volume=str(tailvolume),heavy_tail_volume_over_E=str(tailvolume/X),heavy_tail_B_over_X=str(tailbase/X)))
        print('ENHANCED',L,'S_bound',bound,'minimum',minimum,'tailvol/E',float(tailvolume/X),'tailB/X',float(tailbase/X),flush=True)
    return dict(status='passed',arithmetic='Fraction exact',exact_error='0',scope='enhanced chain exact windows plus full original observer heavy tail',cases=results)

def main():
    if args.enhanced_only:
        OUT.write_text(json.dumps(enhanced_windows(),indent=2)+'\n')
        return
    start=time.monotonic();out=dict(status='running',arithmetic='Fraction exact',dimension=1,scope='full original band; a=alpha/8,b=alpha/2; eligible J>=2; all true K<=D',cases=[])
    for batch,label,atoms,alpha,params in cases():
        result=evaluate(batch,label,atoms,alpha,params);out['cases'].append(result)
        OUT.write_text(json.dumps(out)+'\n')
        print(batch,label,'N',result['N'],'B/X',float(F(result['B_over_X'])),'F/X',float(F(result['F_over_X'])),'Q/X',float(F(result['Q_over_X'])),'Smax',float(F(result['max_observer_S'])),'A_source_max',float(F(result['max_source_A'])),'secs',round(result['elapsed_seconds'],2),flush=True)
    out['batches']=[]
    for batch in sorted({c['batch'] for c in out['cases']}):
        group=[c for c in out['cases'] if c['batch']==batch]
        maximum={}
        for key in ('B_over_X','F_over_X','Q_over_X','max_observer_S','max_source_A','full_drift_over_X','coarse_positive_over_X','coarse_negative_over_X'):
            best=max(group,key=lambda c:F(c[key]));maximum[key]=dict(label=best['label'],exact=best[key],display=float(F(best[key])))
        out['batches'].append(dict(batch=batch,count=len(group),maximum=maximum,positive_F_count=sum(F(c['F'])>0 for c in group)))
    out.update(status='passed',case_count=len(out['cases']),elapsed_seconds=time.monotonic()-start,exact_error='0',finite_tests_prove_general_bound=False,
               sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    finalize(out)
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('PASSED',len(out['cases']),'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
