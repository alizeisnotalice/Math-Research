#!/usr/bin/env python3
"""Exact continuous occupation-rank deficits on round44 actual n=1 cases.

No observation sets are frozen or reselected: archived constructors recompute
the original band and first J/K. Only the requested JSON/Markdown are written.
"""
from fractions import Fraction as F
from pathlib import Path
from bisect import bisect_left, bisect_right
from collections import defaultdict
import argparse, importlib.util, sys, json, time, hashlib
sys.dont_write_bytecode = True
SHIFTS = (F(0), F(1,4), F(1,2), F(1), F(3,2))
RANK_SHIFT = F(3,2)
DEFAULT_ROOT = Path(__file__).resolve().parents[2]

def archived44(root):
    spec = importlib.util.spec_from_file_location('archived44_rank', root/'rounds/044/tail_probe.py')
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.ROOT = root
    return mod

def positive_square_integral(x0, x1, width):
    """Independent cubic-antiderivative integral of (affine x)_+ squared."""
    if x0 == x1: return width*max(F(0),x0)**2
    return width*(max(F(0),x1)**3-max(F(0),x0)**3)/(3*(x1-x0))

def positive_affine_integral(x0, x1, width):
    if x0 == x1: return width*max(F(0),x0)
    return width*(max(F(0),x1)**2-max(F(0),x0)**2)/(2*(x1-x0))

def rank_integrals(S,h,C,a0,a1,width):
    """Integrate both rank signs; Simpson exact after threshold splitting.

    Returned values have not yet been multiplied by the true edge density w/V.
    The independent method uses positive-part cubic antiderivatives.
    """
    assert h > 0
    lo,hi = S-C, S+h-C
    if hi <= min(a0,a1):
        return F(0), width*((a0+a1)/2+C-S-h/2), 0
    if lo >= max(a0,a1):
        return width*(S+h/2-C-(a0+a1)/2), F(0), 0
    cuts = {F(0), F(1)}
    if a1 != a0:
        for t in (lo,hi):
            u = (t-a0)/(a1-a0)
            if 0<u<1: cuts.add(u)
    cuts = sorted(cuts)
    def values(u):
        A = a0+(a1-a0)*u
        d = (max(F(0),hi-A)**2-max(F(0),lo-A)**2)/(2*h)
        n = (max(F(0),A-lo)**2-max(F(0),A-hi)**2)/(2*h)
        return d,n
    d=n=F(0)
    for u,v in zip(cuts,cuts[1:]):
        d0,n0=values(u);dm,nm=values((u+v)/2);d1,n1=values(v)
        d += width*(v-u)*(d0+4*dm+d1)/6
        n += width*(v-u)*(n0+4*nm+n1)/6
    direct_d = (positive_square_integral(hi-a0,hi-a1,width)-positive_square_integral(lo-a0,lo-a1,width))/(2*h)
    direct_n = (positive_square_integral(a0-lo,a1-lo,width)-positive_square_integral(a0-hi,a1-hi,width))/(2*h)
    assert d==direct_d and n==direct_n and d>=0 and n>=0
    assert d-n == width*(S+h/2-C-(a0+a1)/2)
    return d,n,len(cuts)-1

def self_test():
    count=0
    for S in (F(0),F(1,3),F(2)):
        for h in (F(1,7),F(1,2),F(1)):
            for C in SHIFTS:
                for a0,a1 in ((F(0),F(0)),(F(0),F(3)),(F(3),F(0)),(F(1,4),F(3,4)),(F(3,4),F(1,4))):
                    d,n,_=rank_integrals(S,h,C,a0,a1,F(5,7))
                    assert d-n==F(5,7)*(S+h/2-C-(a0+a1)/2)
                    count+=1
    return count

def evaluate(old,batch,label,atoms,alpha,params):
    start=time.monotonic();saved,_,radii,_=old.observer(atoms,alpha,'band',4)
    assert F(saved['entrance_a'])==alpha/8 and F(saved['terminal_b'])==alpha/2
    cells=[dict(c,radius=str(radii[c['K']])) for c in saved['observer_cells']]
    assert all(c['J']>=2 for c in cells)
    hs=old.build_h(cells,alpha/8,'eligible');ks=sorted(hs);A=old.add(hs.values())
    ys=[y for y,w in atoms];assert ys==sorted(set(ys));vectors={};av={};source_shift={C:F(0) for C in SHIFTS}
    for y,w in atoms:
        running=F(0);vectors[y]={}
        for k in ks:
            h=hs[k].value(y);assert 0<=h<=1
            vectors[y][k]=(running,h);running+=h
            direct=sum((max(F(0),min(F(c['hi']),y+radii[k])-max(F(c['lo']),y-radii[k]))/(2*radii[k]) for c in cells if c['K']==k),F(0))
            assert direct==h
        av[y]=running;assert A.value(y)==running
        for C in SHIFTS:source_shift[C]+=w*max(F(0),running-C)**2
    X=alpha*F(saved['eligible_volume']);M=sum((w*av[y] for y,w in atoms),F(0));Q=sum((w*av[y]**2 for y,w in atoms),F(0))
    assert X>0 and X/2<=M<=X
    knots=sorted({z for h in hs.values() for z in h.knots});cache={}
    def a_value(z):
        if z not in cache:cache[z]=A.value(z)
        return cache[z]
    D={C:F(0) for C in SHIFTS};N={C:F(0) for C in SHIFTS}
    Cobs=T=edgeM=edgeS=F(0);max_rank_excess=F(0);positive_witness=None
    edge_count=segment_count=quadratic_pieces=0;positive_edges=0
    raw=defaultdict(F)
    for c in cells:
        k=c['K'];l=F(c['lo']);r=F(c['hi']);rad=radii[k];vk=2*rad
        cuts=sorted({l,r}|set(knots[bisect_right(knots,l):bisect_left(knots,r)])|{z for y,w in atoms for z in (y-rad,y+rad) if l<z<r})
        for a,b in zip(cuts,cuts[1:]):
            segment_count+=1;mid=(a+b)/2;width=b-a;a0=a_value(a);a1=a_value(b)
            present=atoms[bisect_right(ys,mid-rad):bisect_left(ys,mid+rad)]
            g=sum((w for y,w in present),F(0))/vk;assert alpha/2<g<=alpha
            Cobs+=width*(a0+a1)*g/2
            for y,w in present:
                edge_count+=1;S,h=vectors[y][k];assert h>0
                mass=width*w/vk;edgeM+=mass;edgeS+=mass*(S+h/2);raw[y]+=width/vk
                excess=S+h-min(a0,a1);assert excess<=RANK_SHIFT
                max_rank_excess=max(max_rank_excess,excess)
                T+=positive_affine_integral(S-a0,S-a1,width)*w/vk
                for C in SHIFTS:
                    d,n,pieces=rank_integrals(S,h,C,a0,a1,width);quadratic_pieces+=pieces
                    D[C]+=d*w/vk;N[C]+=n*w/vk
                    if C==0 and d:
                        positive_edges+=1
                        if positive_witness is None:
                            positive_witness=dict(y=str(y),lo=str(a),hi=str(b),K=k,S=str(S),h=str(h),A_lo=str(a0),A_hi=str(a1),edge_density=str(w/vk),D0_contribution=str(d*w/vk))
    assert edgeM==M and edgeS==Q/2 and all(raw[y]==av[y] for y,w in atoms)
    assert Cobs<=X and D[F(3,2)]==0 and D[F(1)]<=T
    assert all(D[SHIFTS[i]]>=D[SHIFTS[i+1]] for i in range(len(SHIFTS)-1))
    shifted=[]
    for C in SHIFTS:
        signed=Q/2-C*M-Cobs
        ordinary_bound=source_shift[C]/2
        rank_bound=(source_shift[C]-source_shift[RANK_SHIFT])/2
        assert D[C]-N[C]==signed and N[C]<=C*M+Cobs
        assert D[C]<=ordinary_bound and D[C]<=rank_bound and D[C]<=(RANK_SHIFT-C)*M
        assert ordinary_bound<=Cobs+D[C] and Q/2<=C*M+Cobs+D[C]
        shifted.append(dict(C=str(C),D=str(D[C]),N=str(N[C]),signed=str(signed),D_over_X=str(D[C]/X),D_over_M=str(D[C]/M),N_over_X=str(N[C]/X),signed_over_X=str(signed/X),source_shifted_square=str(source_shift[C]),ordinary_moment_bound=str(ordinary_bound),rank_difference_moment_bound=str(rank_bound),rank_mass_bound=str((RANK_SHIFT-C)*M),negative_part_bound=str(C*M+Cobs),shifted_transfer_bound=str(2*(Cobs+D[C])),shifted_transfer_slack=str(2*(Cobs+D[C])-source_shift[C]),identity_error='0'))
    return dict(batch=batch,label=label,parameters=params,N=len(atoms),D_depth=saved['D'],K_count=len(ks),observer_cell_count=len(cells),max_source_A=str(max(av.values())),X=str(X),M=str(M),Q=str(Q),Q_over_X=str(Q/X),Cobs=str(Cobs),Cobs_over_X=str(Cobs/X),signed_Qhalf_minus_Cobs=str(Q/2-Cobs),signed_Qhalf_minus_Cobs_over_X=str((Q/2-Cobs)/X),prefix_deficit_T=str(T),T_over_X=str(T/X),D1_over_T=str(D[F(1)]/T) if T else None,shifted=shifted,positive_D0=bool(D[F(0)]),positive_D0_edge_count=positive_edges,positive_D0_witness=positive_witness,max_rank_excess=str(max_rank_excess),edge_count=edge_count,observer_affine_segments=segment_count,simpson_quadratic_subpieces=quadratic_pieces,exact_error='0',elapsed_seconds=time.monotonic()-start)

def finalize(out):
    out['batches']=[]
    for batch in sorted({c['batch'] for c in out['cases']}):
        group=[c for c in out['cases'] if c['batch']==batch]
        shifts=[]
        for C in SHIFTS:
            rows=[(c,next(r for r in c['shifted'] if F(r['C'])==C)) for c in group]
            maxima={}
            for key in ('D_over_X','D_over_M','N_over_X','signed_over_X'):
                c,r=max(rows,key=lambda v:F(v[1][key]));maxima[key]=dict(label=c['label'],exact=r[key],display=float(F(r[key])))
            shifts.append(dict(C=str(C),positive_count=sum(F(r['D'])>0 for c,r in rows),maxima=maxima))
        best=max(group,key=lambda c:F(c['T_over_X']))
        out['batches'].append(dict(batch=batch,count=len(group),shifted=shifts,T_over_X_maximum=dict(label=best['label'],exact=best['T_over_X'],display=float(F(best['T_over_X'])))))
    out['audit']=dict(exact_error='0',edge_count=sum(c['edge_count'] for c in out['cases']),simpson_quadratic_subpieces=sum(c['simpson_quadratic_subpieces'] for c in out['cases']),source_direct_overlap_reconstruction=True,rank_source_uniform_marginal_Qhalf=True,observer_marginal_Cobs=True,all_defect_signed_identities=True,all_negative_part_bounds=True,all_shifted_moment_bounds=True,all_shifted_moment_transfer_bounds=True,all_D1_le_T=True,all_D_three_halves_zero=True,independent_cubic_antiderivative_check_on_every_mixed_edge=True)
    out['positive_D0_cases']=[c['label'] for c in out['cases'] if c['positive_D0']]
    out['case_count']=len(out['cases']);out['status']='passed';out['sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,default=DEFAULT_ROOT);p.add_argument('--output',type=Path,required=True);p.add_argument('--self-test-only',action='store_true');p.add_argument('--labels',nargs='*');args=p.parse_args()
    tests=self_test()
    if args.self_test_only: print('PASSED local integration tests',tests);return
    started=time.monotonic();round44=archived44(args.repo_root);old=round44.archived43(args.output)
    out=dict(status='running',arithmetic='Fraction exact',dimension=1,scope='actual original band/J/K, continuous rank coupling44',shifts=[str(C) for C in SHIFTS],local_integration_tests=tests,cases=[],finite_tests_prove_dimension_independent_bound=False)
    for batch,label,atoms,alpha,params in round44.cases(old):
        if args.labels is not None and label not in args.labels:continue
        result=evaluate(old,batch,label,atoms,alpha,params);out['cases'].append(result)
        args.output.write_text(json.dumps(out)+'\n')
        print(batch,label,'D/X',[(r['C'],float(F(r['D_over_X']))) for r in result['shifted']],'T/X',float(F(result['T_over_X'])),'seconds',round(result['elapsed_seconds'],2),flush=True)
    finalize(out);out['elapsed_seconds']=time.monotonic()-started;args.output.write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['case_count'],'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
