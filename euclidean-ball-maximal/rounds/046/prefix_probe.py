#!/usr/bin/env python3
"""Bounded exact stress search for positive source-prefix/full-observer deficit.
Uses actual original-band first-J/K constructors; finite n=1 samples only.
"""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left,bisect_right
from collections import defaultdict
import argparse,importlib.util,json,random,time,sys,hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]

def constructors(root,output):
    spec=importlib.util.spec_from_file_location('prefix46_archived44',root/'rounds/044/tail_probe.py')
    old44=importlib.util.module_from_spec(spec);spec.loader.exec_module(old44);old44.ROOT=root
    return old44.archived43(output)

def affine_positive(x0,x1,width):
    if min(x0,x1)>=0:return width*(x0+x1)/2
    if max(x0,x1)<=0:return F(0)
    root=x0/(x0-x1)
    if x0>0:return width*root*x0/2
    return width*(1-root)*x1/2

def cases(seed):
    for batch,count in enumerate((10,20,30)):
        for index in range(count):
            s=seed+10000*batch+index;rng=random.Random(s)
            n=rng.randint(3,8) if batch==0 else rng.randint(5,20)
            raw=[]
            if batch==0:
                # Independent masses, variable separation, some almost-contact atoms.
                centers=[F(rng.randint(-48,48),2**rng.randint(2,7)) for _ in range(n-1)]
                for i,x in enumerate(centers):raw.append((x,F(rng.randint(1,32),2**rng.randint(0,7))))
                family='unequal_small_configuration'
            elif batch==1:
                # 2--4 very narrow clouds, cluster weights chosen independently.
                nc=rng.randint(2,4);centers=[F(rng.randint(-24,24),2**rng.randint(1,5)) for _ in range(nc)]
                scales=[F(1,2**rng.randint(4,12)) for _ in centers]
                masses=[F(rng.randint(1,16),2**rng.randint(0,5)) for _ in centers]
                for i in range(n-1):
                    cluster=i%nc;x=centers[cluster]+rng.randint(-8,8)*scales[cluster]
                    raw.append((x,masses[cluster]*F(rng.randint(1,24),2**rng.randint(0,5))))
                family='narrow_multiple_clouds'
            else:
                # Correlated geometric mass/distance chains in both directions,
                # plus independent heavy atoms near contact with coarse clouds.
                depth=rng.randint(3,7);q=rng.choice((2,3,4,6,8));scale=F(1,2**rng.randint(1,5))
                offset=F(rng.randint(-8,8),16);position=F(rng.randint(4,20),8)
                for i in range(n-3):
                    level=i%depth;r=scale/F(q**level);side=rng.choice((-1,1))
                    x=offset+side*position*r+F(rng.randint(-4,4),32)*r
                    raw.append((x,F(rng.randint(4,64),16)*r))
                for i in range(2):raw.append((offset+F(rng.randint(-48,48),16)*scale,F(rng.randint(2,32),16)*scale))
                family='geometric_contact_clouds'
            merged=defaultdict(F)
            for x,w in raw:merged[x]+=w
            alpha=rng.choice((F(1,4),F(1,2),F(3,4),F(1),F(3,2),F(2),F(4)))
            # Main cloud below half total mass avoids the all-J1 degeneracy.
            # Distant completion stays outside its coarse balls; original J/K
            # is still recomputed, including the completion atom's own band.
            total=sum(merged.values(),F(0));rho=F(rng.randint(4,15),32)
            atoms=sorted([(x,rho*w/total) for x,w in merged.items()]+[(F(100)/alpha,1-rho)])
            assert sum(w for y,w in atoms)==1 and all(w>0 for y,w in atoms)
            yield batch,f'{family}_s{s}',atoms,alpha,dict(seed=s,family=family,raw_N=n,main_cloud_mass=str(rho),probability_completion_location=str(F(100)/alpha))

def evaluate(old,batch,label,atoms,alpha,params):
    began=time.monotonic();saved,_,radii,_=old.observer(atoms,alpha,'band',4)
    assert F(saved['entrance_a'])==alpha/8 and F(saved['terminal_b'])==alpha/2
    cells=[dict(c,radius=str(radii[c['K']])) for c in saved['observer_cells']]
    assert all(c['J']>=2 for c in cells)
    hs=old.build_h(cells,alpha/8,'eligible');ks=sorted(hs);A=old.add(hs.values())
    vectors={};av={};source_overlap_count=0
    by_k=defaultdict(list)
    for c in cells:by_k[c['K']].append(c)
    for y,w in atoms:
        v={};S=F(0)
        for k in ks:
            h=hs[k].value(y);assert 0<=h<=1
            direct=sum((max(F(0),min(F(c['hi']),y+radii[k])-max(F(c['lo']),y-radii[k]))/(2*radii[k]) for c in by_k[k]),F(0))
            assert direct==h;source_overlap_count+=1;v[k]=(S,h);S+=h
        vectors[y]=v;av[y]=S;assert A.value(y)==S
    X=alpha*F(saved['eligible_volume']);M=sum((w*av[y] for y,w in atoms),F(0))
    assert X/2<=M<=X
    knots=sorted({z for h in hs.values() for z in h.knots});ys=[y for y,w in atoms]
    cache={};T=mass=Cobs=F(0);raw=defaultdict(F);edge_count=positive_count=0
    max_excess=None;witness=None
    def value(z):
        if z not in cache:cache[z]=A.value(z)
        return cache[z]
    for c in cells:
        l=F(c['lo']);r=F(c['hi']);k=c['K'];rad=radii[k];v=2*rad
        cuts=sorted({l,r}|set(knots[bisect_right(knots,l):bisect_left(knots,r)])|{z for y,w in atoms for z in (y-rad,y+rad) if l<z<r})
        for a,b in zip(cuts,cuts[1:]):
            mid=(a+b)/2;width=b-a;a0=value(a);a1=value(b)
            present=atoms[bisect_right(ys,mid-rad):bisect_left(ys,mid+rad)]
            g=sum((w for y,w in present),F(0))/v;assert alpha/2<g<=alpha
            Cobs+=width*g*(a0+a1)/2
            for y,w in present:
                edge_count+=1;S,h=vectors[y][k];assert h>0
                mass+=width*w/v;raw[y]+=width/v
                excess=max(S-a0,S-a1)
                if max_excess is None or excess>max_excess:max_excess=excess
                integral=affine_positive(S-a0,S-a1,width)
                independent=width*max(F(0),S-a0) if a0==a1 else width*(max(F(0),S-a1)**2-max(F(0),S-a0)**2)/(2*(a0-a1))
                assert integral==independent and integral>=0
                contribution=integral*w/v;T+=contribution
                if contribution:
                    positive_count+=1
                    if witness is None:
                        witness=dict(y=str(y),source_mass=str(w),alpha=str(alpha),lo=str(a),hi=str(b),J=c['J'],K=k,radius=str(rad),S=str(S),h=str(h),A_lo=str(a0),A_hi=str(a1),positive_prefix_excess_max=str(excess),edge_density=str(w/v),T_contribution=str(contribution),parent_observer_cell=c)
    assert mass==M and all(raw[y]==av[y] for y,w in atoms)
    return dict(batch=batch,label=label,parameters=params,N=len(atoms),alpha=str(alpha),D=saved['D'],cell_count=len(cells),K_count=len(ks),X=str(X),M=str(M),Cobs=str(Cobs),original_volume=saved['original_volume'],J1_volume=saved['J1_volume'],T=str(T),T_over_X=str(T/X) if X else None,max_source_A=str(max(av.values())),max_prefix_minus_observer_A=str(max_excess) if max_excess is not None else None,positive_edge_count=positive_count,edge_count=edge_count,source_overlap_checks=source_overlap_count,witness=witness,atoms=[dict(location=str(y),mass=str(w)) for y,w in atoms],exact_error='0',elapsed_seconds=time.monotonic()-began)

def independent_witness(c):
    """Rebuild the centered maximal band and every J/K without observer()."""
    atoms=[(F(a['location']),F(a['mass'])) for a in c['atoms']];alpha=F(c['alpha'])
    D=c['D']+1;radii={j:F(4,2**j)/alpha for j in range(1,D+1)}
    cuts={y for y,w in atoms}
    for y,w in atoms:
        for r in radii.values():cuts.update((y-r,y+r))
    for i in range(len(atoms)):
        mass=F(0)
        for j in range(i,len(atoms)):
            mass+=atoms[j][1]
            for beta in (alpha,2*alpha):cuts.update((atoms[j][0]-mass/(2*beta),atoms[i][0]+mass/(2*beta)))
    cells=[];volume=J1=F(0)
    cuts=sorted(cuts)
    for lo,hi in zip(cuts,cuts[1:]):
        mid=(lo+hi)/2;dist=defaultdict(F)
        for y,w in atoms:dist[abs(mid-y)]+=w
        mass=mx=F(0)
        for d in sorted(dist):mass+=dist[d];mx=max(mx,mass/(2*d))
        if not alpha<mx<=2*alpha:continue
        volume+=hi-lo
        g={j:sum((w for y,w in atoms if abs(mid-y)<rad),F(0))/(2*rad) for j,rad in radii.items()}
        J=next(j for j in g if g[j]>alpha/8);K=next(j for j in g if g[j]>alpha/2)
        if J==1:J1+=hi-lo
        else:cells.append((lo,hi,J,K))
    assert volume==F(c['original_volume']) and J1==F(c['J1_volume'])
    assert alpha*sum((hi-lo for lo,hi,J,K in cells),F(0))==F(c['X'])
    def h(x,k):
        rad=radii[k]
        return sum((max(F(0),min(hi,x+rad)-max(lo,x-rad))/(2*rad) for lo,hi,J,K in cells if K==k),F(0))
    ks=sorted({K for lo,hi,J,K in cells})
    M=sum((w*sum((h(y,k) for k in ks),F(0)) for y,w in atoms),F(0));assert M==F(c['M'])
    w=c['witness'];y=F(w['y']);lo=F(w['lo']);hi=F(w['hi']);k=w['K']
    matching=[(max(lo,a),min(hi,b)) for a,b,J,K in cells if J==w['J'] and K==k and max(lo,a)<min(hi,b)]
    assert sum((b-a for a,b in matching),F(0))==hi-lo
    S=sum((h(y,j) for j in ks if j<k),F(0));a0=sum((h(lo,j) for j in ks),F(0));a1=sum((h(hi,j) for j in ks),F(0))
    assert S==F(w['S']) and h(y,k)==F(w['h']) and a0==F(w['A_lo']) and a1==F(w['A_hi'])
    contribution=affine_positive(S-a0,S-a1,hi-lo)*F(w['source_mass'])/(2*radii[k])
    assert contribution==F(w['T_contribution']) and contribution>0
    return dict(status='passed',original_band_and_all_J_K_reconstructed=True,original_and_J1_volumes_exact=True,source_marginal_M_exact=True,witness_direct_overlap_S_A_exact=True,positive_contribution=str(contribution),audit_D=D,independent_cell_count=len(cells),exact_error='0')

def interval_certificate(old,c):
    """Compact full eligible E_k lists and constant witness h profiles."""
    atoms=[(F(v['location']),F(v['mass'])) for v in c['atoms']];alpha=F(c['alpha'])
    saved,_,radii,_=old.observer(atoms,alpha,'band',4)
    by_k=defaultdict(list)
    for row in saved['observer_cells']:by_k[row['K']].append((F(row['lo']),F(row['hi'])))
    merged={}
    for k,intervals in by_k.items():
        result=[]
        for lo,hi in sorted(intervals):
            if result and lo==result[-1][1]:result[-1]=(result[-1][0],hi)
            else:result.append((lo,hi))
        merged[k]=result
    def h(x,k):
        r=radii[k]
        return sum((max(F(0),min(hi,x+r)-max(lo,x-r))/(2*r) for lo,hi in merged[k]),F(0))
    w=c['witness'];y=F(w['y']);lo=F(w['lo']);hi=F(w['hi']);rows=[]
    for k,intervals in sorted(merged.items()):
        nodes=sorted({lo,hi}|{z for a,b in intervals for z in (a-radii[k],a+radii[k],b-radii[k],b+radii[k]) if lo<z<hi})
        values=[h(z,k) for z in nodes]
        assert len(set(values))==1
        rows.append(dict(K=k,radius=str(radii[k]),E_intervals=[[str(a),str(b)] for a,b in intervals],E_volume=str(sum((b-a for a,b in intervals),F(0))),source_h=str(h(y,k)),observer_h_constant=str(values[0]),observer_h_knots_checked=len(nodes)))
    assert sum((F(r['source_h']) for r in rows if r['K']<w['K']),F(0))==F(w['S'])
    assert sum((F(r['observer_h_constant']) for r in rows),F(0))==F(w['A_lo'])==F(w['A_hi'])
    assert F(c['T'])==F(w['T_contribution']) and c['positive_edge_count']==1
    return dict(status='passed',all_eligible_J_values=sorted({row['J'] for row in saved['observer_cells']}),E_by_K=rows,global_T=str(c['T']),positive_edge_segment_count=c['positive_edge_count'],all_other_edge_segment_deficits_zero=True,global_T_equals_witness_contribution=True,exact_error='0')

def summary(out):
    witnesses=[c for c in out['cases'] if F(c['T'])>0]
    out['positive_case_count']=len(witnesses);out['case_count']=len(out['cases']);out['batches']=[]
    for batch in sorted({c['batch'] for c in out['cases']}):
        group=[c for c in out['cases'] if c['batch']==batch];active=[c for c in group if c['T_over_X'] is not None]
        best=max(active,key=lambda c:F(c['T_over_X'])) if active else None
        out['batches'].append(dict(batch=batch,count=len(group),positive_count=sum(F(c['T'])>0 for c in group),eligible_empty_count=sum(F(c['X'])==0 for c in group),max_T_over_X=dict(label=best['label'],exact=best['T_over_X'],display=float(F(best['T_over_X']))) if best else None))
    out['exact_checks']=dict(error='0',source_direct_overlap_count=sum(c['source_overlap_checks'] for c in out['cases']),edge_segment_count=sum(c['edge_count'] for c in out['cases']),mu_source_mass_reconstruction=True,observer_cap=True,independent_positive_affine_integrals=True)
    out['sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,default=461007);p.add_argument('--time-budget',type=float,default=120);args=p.parse_args()
    old=constructors(args.repo_root,args.output);began=time.monotonic()
    out=dict(status='running',arithmetic='Fraction exact',dimension=1,scope='actual original band/J/K; prefix/full-observer positive deficit',seed=args.seed,planned_counts=[10,20,30],time_budget_seconds=args.time_budget,cases=[],finite_search_proves_general_zero=False)
    for batch,label,atoms,alpha,params in cases(args.seed):
        if time.monotonic()-began>args.time_budget:break
        c=evaluate(old,batch,label,atoms,alpha,params);out['cases'].append(c);args.output.write_text(json.dumps(out)+'\n')
        print(batch,label,'N',c['N'],'T/X',c['T_over_X'],'max_excess',c['max_prefix_minus_observer_A'],'sec',round(c['elapsed_seconds'],2),flush=True)
    out['status']='passed' if len(out['cases'])==60 else 'budget_completed';out['elapsed_seconds']=time.monotonic()-began;summary(out)
    positives=[c for c in out['cases'] if F(c['T'])>0]
    for c in positives:c['independent_audit']=independent_witness(c)
    minimal=evaluate(old,3,'minimal_four_atoms',[(F(0),F(13,256)),(F(1,16),F(3,64)),(F(3,32),F(9,128)),(F(10),F(213,256))],F(1),{})
    minimal['independent_audit']=independent_witness(minimal);minimal['interval_certificate']=interval_certificate(old,minimal)
    out['minimal_witness']=minimal
    out['total_elapsed_seconds']=time.monotonic()-began
    args.output.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['case_count'],'positive',out['positive_case_count'],'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
