#!/usr/bin/env python3
"""Saved static-section, two position derivatives and signed symmetry audit."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations,product
from math import prod
from collections import Counter
import sys,os,time,gzip,hashlib
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path[:0]=[str(ROOT/'lower'),str(P)]
from certify_rectangular import field
from position_sections import ENV,load,save,sha,rt,enc
def main():
    start=time.monotonic();F,u,v,bounds,_,_=field();a=1+u**3/2;s=(1+u)/2;zero=F(0)
    def dec(d):
        z=[Q(0)]*25
        for ij,c in d.items():i,j=map(int,ij.split(','));z[5*i+j]=Q(c)
        return F(a=z)
    reg=load(P/'registration.json');repair=load(P/'repair_registration.json')
    receipt=load(P/'prefix_archive_receipt.json')
    for source,rr in [('position_sections.py',reg),('resume_projected_sections.py',repair)]:
        assert sha(P/source)==rr['source_sha256']
        for x in rr['inputs']:
            if Path(x['path']).name=='tail_checkpoint.json' and not Path(x['path']).exists():assert x['sha256']==receipt['uncompressed_sha256']
            else:assert sha(Path(x['path']))==x['sha256']
    for oldname in ('lower_round023','lower_round024'):
        old=ROOT/oldname
        for x in load(old/'archive_whitelist.json')['files']:assert sha(old/x['path'])==x['sha256']
    packpath=P/receipt['compressed_path'];local_prefix_hash_checked=False
    if packpath.exists():
        packed=packpath.read_bytes();assert hashlib.sha256(packed).hexdigest()==receipt['compressed_sha256']
        assert hashlib.sha256(gzip.decompress(packed)).hexdigest()==receipt['uncompressed_sha256'];local_prefix_hash_checked=True
    geo=load(P/'tail_geometry.json');projection=load(P/'projection_certificate.json');hh=load(P/'tail_histogram.json')
    assert len(geo['raw_inventory'])==1024 and len({tuple(x['lengths']) for x in geo['raw_inventory']})==1024
    keep=projection['kept_first_span_families'];keepmask=sum(1<<i for i in keep);hist={int(m):dec(z) for m,z in hh['histogram'].items()}
    assert all(m&~keepmask==0 for m in hist) and sum(hist.values(),F(0))==(s/2)**4
    assert hh['old_cells_recomputed']==0 and hh['successfully_reused_fine_cells']==193336 and hh['new_remaining_coarse_cells']==36864
    assert a-F(Q(1,4))>zero
    families={x['family_id']:x for x in geo['families']};atomchecks=0;marginals=0
    for fid in keep:
        f=families[fid];ls=f['lengths'];n=prod(ls);r=dec(f['radius']);assert 32*r**5==F(n)
        box=list(map(dec,f['tail_box']));assert (ls[0]-1)*a<=2*r
        target=F(1)
        for k in range(4):target*=box[2*k+1]-box[2*k]
        acc=[Q(0)]*25
        for m,w in hist.items():
            if m&(1<<fid):acc=[x+y for x,y in zip(acc,w.a)]
        assert F(a=acc)==target;marginals+=1
        tailstarts=[-(ell-1)//2 if ell%2 else 1-ell//2 for ell in ls[1:]]
        for j in range(-1,5):
            starts=[j]+tailstarts;atoms=set(product(*(range(st,st+ell) for st,ell in zip(starts,ls))))
            assert len(atoms)==n;atomchecks+=n
    emax=Q(projection['valid_epsilon_upper']);assert 0<emax<=Q(1,8)
    for x in projection['excluded_strict_span_families']:
        f=families[x['family_id']];r=dec(f['radius']);gap=(f['lengths'][0]-1)*a-2*r
        assert 2*F(emax)<gap and dec(x['strict_negative_first_axis_width'])==-gap
    # Independent factorization of the entire final static histogram. Tail axes
    # are reversed, then integrated as two 2D signature histograms and one AND.
    # This small certificate replay removes any need to archive the31MB fine prefix.
    vr=load(P/'verification_registration.json');assert sha(Path(__file__))==vr['source_sha256']
    for row in vr['inputs']:assert sha(Path(row['path']))==row['sha256']
    boxes={fid:list(map(dec,families[fid]['tail_box'])) for fid in keep};segments=[]
    for k in (3,2,1,0):
        cuts=sorted({zero,s/2}|{b[2*k+j] for b in boxes.values() for j in (0,1)})
        segments.append([(hi-lo,sum(1<<fid for fid,b in boxes.items() if b[2*k]<=lo and hi<=b[2*k+1])) for lo,hi in zip(cuts,cuts[1:])])
    pairh=[]
    for offset in (0,2):
        hh2={}
        for (dx,mx),(dy,my) in product(segments[offset],segments[offset+1]):
            m=mx&my;hh2[m]=hh2.get(m,F(0))+dx*dy
        pairh.append(hh2)
    assert len(pairh[0])*len(pairh[1])<=vr['max_joint_pairs']
    independent={};jointpairs=0
    for m,w in pairh[0].items():
        for n,vv in pairh[1].items():
            key=m&n;independent[key]=independent.get(key,F(0))+w*vv;jointpairs+=1
    assert independent==hist
    derivatives=[];checks=[]
    for name in ('plus','minus'):
        result=load(P/f'{name}_derivative.json');ck=load(P/f'{name}_checkpoint.json')
        assert ck['status']=='COMPLETE_all_section_derivatives' and ck['completed_sections']==ck['total_sections']==len(hist)==2224
        D=F(0);speeds=Counter()
        for row in ck['sections']:
            m=int(row['family_mask']);speed=Q(row['oneD_union_length_derivative']);D+=hist[m]*speed;speeds[str(speed)]+=1
        assert D==dec(ck['partial_derivative_volume'])
        ratioD=D*Q(16,3);assert ratioD==dec(result['directional_derivative_coefficients']);derivatives.append(ratioD)
        assert list(map(str,bounds(ratioD.a)))==result['derivative_interval'] and Q(result['derivative_interval'][1])<0
        assert result['counts']['raw_candidates']==6144 and result['counts']['active_raw']==179
        assert result['both_oneD_union_engines_agree_all_sections'] and not result['new_zero_width_openings']
        checks.append({'profile':name,'section_velocity_counts':dict(speeds),'derivative_interval':result['derivative_interval'],'result_sha256':sha(P/f'{name}_derivative.json')})
    assert derivatives[0]==derivatives[1]
    # G={even permutation, minus odd permutation}; a genuinely signed D3 action.
    perms=list(permutations(range(3)))
    def parity(p):return -1 if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2 else 1
    def act(z,p):return tuple(parity(p)*z[p[i]] for i in range(3))
    plus=(Q(1),Q(0),Q(-1));minus=tuple(-x for x in plus)
    po={act(plus,p) for p in perms};mo={act(minus,p) for p in perms};assert len(po)==len(mo)==3 and po.isdisjoint(mo)
    tforms=((Q(1,2),Q(-1,2)),(Q(0),Q(1)),(Q(-1,2),Q(-1,2)))
    flipped=tuple(tuple(-x for x in tforms[i]) for i in (2,1,0))
    assert flipped==tuple((A,-B) for A,B in tforms)
    out={'status':'PASS_two_period3_position_directions_and_static_sections','profiles':checks,'same_exact_derivative':True,'static_section_count':len(hist),'complete1024inventory':True,'eligible_tail_marginals_verified':marginals,'explicit_distinct_atom_checks':atomchecks,'safe_ordering_at_epsilon_1_8':True,'signed_D3_plus_minus_orbits_disjoint':True,'next_fan_qualification_only':'two halfedges +/-((1-t)/2,t,-(1+t)/2),0<=t<=1/3 suffice under signedG;no events orfan computed','local_compressed_prefix_hash_checked':local_prefix_hash_checked,'entire_static_histogram_independent_pair_factorization_agrees':True,'independent_static_joint_pairs':jointpairs,'fine_prefix_not_required_for_final_certificate':True,'old023024_archive_files_unchanged':True,'new_5D_or_parameter_union_jobs_in_verifier':0,'runtime':rt(start),'source_sha256':sha(Path(__file__))}
    save(P/'verification.json',out);print({k:v for k,v in out.items() if k!='profiles'})
if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
