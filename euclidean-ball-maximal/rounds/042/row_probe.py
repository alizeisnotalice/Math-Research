#!/usr/bin/env python3
"""Exact one-dimensional actual-band terminal-row testing, ratio4 throughout."""
from fractions import Fraction as F
from collections import defaultdict
from bisect import bisect_left
from pathlib import Path
import sys,json,time,random
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'041'))
import square_probe as old
from stage24_clock_replacement_probe import observer
from stage30_transition_weight_probe import build_h
OUT=old.OUT  # shared --output argument; no prior-round main is executed

def normalized(raw,alpha=F(1)):
    atoms=defaultdict(F)
    for x,w in raw:atoms[x/alpha]+=w
    assert all(w>0 for w in atoms.values()) and sum(atoms.values())<1
    atoms[F(12)/alpha]+=1-sum(atoms.values())
    return sorted(atoms.items())

def microcloud(L,r=F(1,32),q=8,sign=1,seed=None):
    rng=random.Random(seed);raw=[]
    for i in range(L):
        ri=r/F(q**i);w=F(27,16)*ri
        if seed is not None:w*=F(rng.randint(40,88),64)
        raw.append((sign*F(9,8)*ri,w))
    raw.append((F(0),r/F(q**(L-1))/16))
    return raw

def cases():
    for label,L in [('old41_chain_L48',48)]:
        yield 0,label,old.chain(L,den=4,position=F(9,8),weight=F(27,16)),F(1),{}
    yield 0,'old41_nested_d5',old.finish(old.tree(5)),F(1),{}
    # Threshold mass and double-radius contact of a big source with a multilevel cloud.
    for i,(L,bits,shift,side) in enumerate([(6,8,-1,1),(12,12,-1,1),(18,20,-1,1),(24,28,-1,1),
        (12,12,1,1),(18,20,1,1),(12,12,-1,-1),(18,20,-1,-1),
        (12,16,-2,1),(18,24,2,-1),(24,32,0,1),(24,32,0,-1)]):
        r=F(1,4);eps=F(1,2**bits);raw=microcloud(L,r/16,sign=side)
        raw.append((side*(2*r+shift*eps*r),r-eps*r))
        yield 1,f'big_contact_L{L}_b{bits}_s{shift}_side{side}',normalized(raw),F(1),dict(L=L,epsilon=str(eps),shift=shift,side=side)
    for L,bits,sign in [(6,10,1),(12,20,1),(24,40,1),(12,20,-1),(24,40,-1),(18,28,1)]:
        r=F(1,4);eps=F(1,2**bits);raw=microcloud(L,r/16,sign=sign)
        # For this one-sided chain the extreme alpha superlevel edge belongs to its biggest atom.
        x,w=raw[0];edge=x+sign*w/2
        raw.append((sign*r+edge-sign*eps*r,r-eps*r))
        yield 1,f'band_edge_sliver_L{L}_b{bits}_side{sign}',normalized(raw),F(1),dict(L=L,epsilon=str(eps),side=sign,target_coarse_j=4,edge=str(edge))
    # Unequal random weights, two offset clouds, almost dyadic big mass; fixed fresh seeds.
    for i in range(14):
        seed=420701+i*19;rng=random.Random(seed);L=[6,10,14,18][i%4];r=F(1,2**(3+i%3));eps=F(1,2**(8+i))
        raw=microcloud(L,r/12,sign=(-1 if i%2 else 1),seed=seed)
        raw += [(F(3,8)*r+x,F(3,5)*w) for x,w in microcloud(max(3,L//2),r/40,sign=-1,seed=seed+1)]
        raw.append(((F(rng.randint(8,24),8)+eps)*r,(1+(-1 if i%3 else 1)*eps)*r))
        alpha=[F(1),F(3,4),F(5,4)][i%3]
        yield 2,f'random_two_clouds_s{seed}_L{L}',normalized(raw,alpha),alpha,dict(seed=seed,L=L,epsilon=str(eps))
    # Asymmetric nested clouds; shifted masses around independent contact scales.
    for i in range(12):
        depth=[2,3,4,5][i%4];seed=420901+23*i;r=F(1,8)
        raw=old.tree(depth,q=[6,8,12][i%3],position=F(9,8),weight=F(27,16),jitter=True,seed=seed,leaf=F(1,8),skew=True)
        eps=F(1,2**(10+i));raw.append(((F(3,2)+(-1)**i*eps)*r,(F(1,2)+(-1)**(i//2)*eps)*r))
        raw.append((-(F(5,2)-eps)*r,F(1,3)*r))
        alpha=[F(1),F(7,8),F(9,8)][i%3]
        yield 3,f'asymmetric_tree_d{depth}_s{seed}',normalized(raw,alpha),alpha,dict(depth=depth,seed=seed,epsilon=str(eps))

def evaluate(atoms,alpha,label):
    saved,en,radii,H=observer(atoms,alpha,'band',4)
    a=F(saved['entrance_a']);b=F(saved['terminal_b']);V=F(saved['eligible_volume']);J1=F(saved['J1_volume']);volume=F(saved['original_volume'])
    assert a==alpha/8 and b==alpha/2 and saved['b_over_a']==4
    assert volume==J1+V and all(v['J']>=2 for v in saved['observer_cells'])
    cells=[dict(v,radius=str(radii[v['K']])) for v in saved['observer_cells']]
    hs=build_h(cells,a,'eligible');ks=sorted(hs)
    rows=[];I={};vectors=[(y,w,{k:h.value(y) for k,h in hs.items()}) for y,w in atoms]
    Xby=defaultdict(F)
    for c in cells:Xby[c['K']]+=alpha*(F(c['hi'])-F(c['lo']))
    for j in ks:
        Mj=sum((w*vs[j] for y,w,vs in vectors),F(0));Rj=F(0);pair=[]
        assert Xby[j]/2<=Mj<=Xby[j]
        for k in ks:
            if k<=j:continue
            v=sum((w*vs[j]*vs[k] for y,w,vs in vectors),F(0));Rj+=v;I[(j,k)]=v
            if v:pair.append(dict(k=k,gap=k-j,I=str(v),I_over_Xj=str(v/Xby[j])))
        row=dict(j=j,radius=str(radii[j]),Ej_volume=str(Xby[j]/alpha),Ej_volume_over_radius=str(Xby[j]/alpha/radii[j]),Xj=str(Xby[j]),Mj=str(Mj),Rj=str(Rj),
            Rj_over_Xj=str(Rj/Xby[j]),pairs=pair,display=float(Rj/Xby[j]))
        rows.append(row)
    Q=sum((w*sum(vs.values())**2 for y,w,vs in vectors),F(0));M=sum((w*sum(vs.values()) for y,w,vs in vectors),F(0))
    diag=sum((w*sum(v*v for v in vs.values()) for y,w,vs in vectors),F(0))
    assert Q==diag+2*sum(I.values()) and sum(Xby.values())==alpha*V
    best=max(rows,key=lambda r:F(r['Rj_over_Xj'])) if rows else None
    summary=dict(label=label,N=len(atoms),D=saved['D'],alpha=str(alpha),a=str(a),b=str(b),b_over_a=4,
        original_volume=str(volume),J1_volume=str(J1),eligible_volume=str(V),J1_separate_volume_assertion=True,
        eligible_J_min=min((c['J'] for c in cells),default=None),cell_count=len(cells),row_count=len(rows),M=str(M),Q=str(Q),X=str(alpha*V),
        Q_over_X=str(Q/(alpha*V)) if V else None,max_row=best,rows=rows,
        max_row_display=best['display'] if best else 0,source=[dict(location=str(y),mass=str(w)) for y,w in atoms])
    cert=dict(summary=summary,observer_cells=cells,original_intervals=saved['original_intervals'])
    return summary,cert

def independent(cert):
    """Independent distance maximal function, original band, all J/K, and row overlap sums."""
    s=cert['summary'];atoms=[(F(v['location']),F(v['mass'])) for v in s['source']];alpha=F(s['alpha']);D=s['D']+1
    radii={j:F(4,2**j)/alpha for j in range(1,D+1)};cuts={y for y,w in atoms}
    for y,w in atoms:
        for r in radii.values():cuts.update([y-r,y+r])
    for i in range(len(atoms)):
        m=F(0)
        for j in range(i,len(atoms)):
            m+=atoms[j][1]
            for beta in [alpha,2*alpha]:cuts.update([atoms[j][0]-m/(2*beta),atoms[i][0]+m/(2*beta)])
    cs=sorted(cuts);cells=defaultdict(list);volume=J1=F(0)
    for l,r in zip(cs,cs[1:]):
        x=(l+r)/2;dist=defaultdict(F)
        for y,w in atoms:dist[abs(x-y)]+=w
        points=sorted(dist);prefix=[F(0)];mx=F(0)
        for z in points:prefix.append(prefix[-1]+dist[z]);mx=max(mx,prefix[-1]/(2*z))
        if not alpha<mx<=2*alpha:continue
        volume+=r-l;g={j:prefix[bisect_left(points,rad)]/(2*rad) for j,rad in radii.items()}
        J=next(j for j in g if g[j]>alpha/8);K=next(j for j in g if g[j]>alpha/2)
        if J==1:J1+=r-l
        else:cells[K].append((l,r))
    def h(x,k):return sum((max(F(0),min(r,x+radii[k])-max(l,x-radii[k]))/(2*radii[k]) for l,r in cells[k]),F(0))
    vectors=[(y,w,{k:h(y,k) for k in cells}) for y,w in atoms]
    for row in s['rows']:
        j=row['j'];Xj=alpha*sum(r-l for l,r in cells[j]);Rj=sum(w*vs[j]*sum(v for k,v in vs.items() if k>j) for y,w,vs in vectors)
        assert Xj==F(row['Xj']) and Rj==F(row['Rj'])
    assert volume==F(s['original_volume']) and J1==F(s['J1_volume'])
    return dict(status='passed',all_rows_Xj_Rj_exactly_recomputed=True,original_band_and_J_K_reconstructed=True,
        J1_and_eligible_independently_separated=True,audit_D=D,row_count=len(s['rows']),maximum=s['max_row'])

def main():
    start=time.monotonic();out=dict(status='running',dimension=1,arithmetic='Fraction exact',I_jk_definition='integral h_j h_k dP, no factor2',
        candidate='max_j (sum_k>j I_jk)/(alpha |E_j|), eligible J>=2',cases=[],retained={});best=None
    for batch,label,atoms,alpha,params in cases():
        s,c=evaluate(atoms,alpha,label);s.update(batch=batch,parameters=params);out['cases'].append(s)
        if best is None or F(s['max_row']['Rj_over_Xj'])>F(best['summary']['max_row']['Rj_over_Xj']):best=c
        out['retained']['maximum_row_case']=best
        OUT.write_text(json.dumps(out)+'\n');print(batch,label,'max row',s['max_row_display'],'j',s['max_row']['j'] if s['max_row'] else None,'N',s['N'],flush=True)
    out['independent_maximum_case_verification']=independent(best)
    out.update(status='passed',elapsed_seconds=time.monotonic()-start,new_batches=3,new_case_count=sum(s['batch']>0 for s in out['cases']),
        maximum_row=best['summary']['max_row'],maximum_row_label=best['summary']['label'],uniform_row_bound_proved=False,
        uniform_row_counterexample_found=False)
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('saved',OUT,'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
