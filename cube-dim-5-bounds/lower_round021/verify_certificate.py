#!/usr/bin/env python3
"""Audit saved certificates without rerunning geometry, histograms, or AND."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import json,hashlib,os,time
from finite_union import G,A,inds,volume,ENV,runtime,dec
P=Path(__file__).resolve().parent
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    start=time.monotonic()
    repair=load(P/'repair_registration.json');reg=load(P/'registration.json')
    for row in repair['locked_inputs']:assert sha(P/row['path'])==row['sha256']
    for row in reg['inputs']:assert sha(P.parent/row['path'])==row['sha256']
    assert sha(P/'continue_and.py')==repair['source_sha256']
    geom=load(P/'candidate_geometry.json');r=load(P/'results.json');ck=load(P/'and_checkpoint.json')
    raw=load(P.parent/'lower_round017/raw_endpoint_catalog.json')['raw_cases']
    assert len(raw)==len(geom['raw_inventory'])==1024
    for a,b in zip(raw,geom['raw_inventory']):
        assert (a['shape'],a['K'],a['n'])==(b['shape'],b['K'],b['n'])
    sl,sh=map(Q,geom['s_interval']);guard_count=0
    for f in geom['families']:
        rl=Q(f['strict_common_radius']);assert (2*rl)**5<f['N']==f['K']*f['n']
        assert len(f['shapes'])==len(f['inner_boxes'])==len(f['outer_boxes'])
        for shape,ib,ob in zip(f['shapes'],f['inner_boxes'],f['outer_boxes']):
            assert all(ob[2*j]<=ib[2*j] and ib[2*j+1]<=ob[2*j+1] for j in range(3))
            for j,(ell,t,sc) in enumerate(zip(shape,(f['alpha_grid'],f['rho_grid'],f['rho_grid']),(A,sh,sh))):
                ii=inds(ell);assert len(ii)==ell
                if volume(ib):
                    assert Q(max(ii))-Q(ib[2*j],G)<=Q(t,G)
                    assert Q(ib[2*j+1],G)-Q(min(ii))<=Q(t,G)
                    assert Q(t,G)*sc<=rl;guard_count+=2
    atoms=load(P.parent/'lower_round019/atoms.json')
    assert len(atoms)==len({(x['m'],x['n']) for x in atoms})==len({(x['x_grid'],x['y_grid']) for x in atoms})==225
    assert ck['status']=='COMPLETE_saved_histogram_AND'
    assert ck['completed_pairs']==ck['total_pairs']==262*649==r['completed_AND_pairs']
    assert ck['next_pair_indices']==[262,0] and ck['joint_histogram']==r['joint_histogram']
    assert Q(sum(r['joint_histogram'].values()),G**5)==Q(r['normalized_inner_stripe_union_volume'])
    assert sum(ck['joint_histogram'].values())==ck['union_grid_volume']
    C=load(P.parent/'lower_round003/a2_equals_a3_results.json')
    D=int(C['root_isolation_denominator'])
    roots=[(Q(n,D),Q(n+1,D)) for n in C['root_lower_numerators']]
    for (lo,hi),n in zip(roots,[2,3]):assert lo**5<n<hi**5
    cl=ch=Q(0)
    for ij,c in C['coefficients'].items():
        i,j=map(int,ij.split(','));c=Q(c)
        lo=roots[0][0]**i*roots[1][0]**j;hi=roots[0][1]**i*roots[1][1]**j
        cl+=c*(lo if c>=0 else hi);ch+=c*(hi if c>=0 else lo)
    assert [str(cl),str(ch)]==C['ratio_interval']
    low,high=map(Q,r['nominal_finite_family_ratio_interval'])
    old=load(P.parent/'lower_round016/delivery007/result.json');ol,oh=map(Q,old['ratio_interval'])
    assert low>oh and high<cl
    assert [str(low-ch),str(high-cl)]==r['difference_from_C_cert_interval']
    assert [str(low-oh),str(high-ol)]==r['difference_from_old016_interval']
    out={'status':'PASS_saved_certificate_audit','all_locked_original021_and017020_inputs_unchanged':True,'raw_status_counts':dict(Counter(x['status'] for x in geom['raw_inventory'])),'common_radius_unrotated_endpoint_guards_checked':guard_count,'distinct_planar_atoms':225,'complete_AND_prefix':ck['completed_pairs'],'C_cert_roots_and25term_interval_reconstructed_exactly':True,'strict_low_minus_old016_high_positive':str(low-oh),'strict_C_cert_low_minus_nominal_high_positive':str(cl-high),'C_gap_decimal_display':[dec(cl-high),dec(ch-low)],'old016_improvement_decimal_display':[dec(low-oh),dec(high-ol)],'inner_union_ratio_decimal_display':[dec(Q(x)) for x in r['inner_union_ratio_interval']],'recomputed_geometry_histogram_AND_cells':0,'runtime':runtime(start),'source_sha256':sha(Path(__file__)),'result_sha256':sha(P/'results.json')}
    (P/'verification.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    assert all(os.environ.get(k)=='1' for k in ENV) and os.getpriority(os.PRIO_PROCESS,0)>=10
    main()
