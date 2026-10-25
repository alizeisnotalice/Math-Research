#!/usr/bin/env python3
"""Check saved exact certificates; no Arb evaluation or spatial kernel calls."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import json,hashlib,time,os
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bounds(r):
    a=Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent'])
    b=Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
    assert a<=b and b-a<=max(Q(1),abs((a+b)/2))*Q(2)**-450
    return a,b
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def check():
    tic=time.monotonic();reg=read(HERE/'registration.json');run=read(HERE/'run_receipt.json');res=read(HERE/'results.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for f,h in reg['code_sha256'].items():assert sha(HERE/f)==h,f
    assert run['status']=='COMPLETED' and run['exit_code']==0
    assert run['attempted_configurations']==run['completed_configurations']==6
    assert run['registration_sha256']==sha(HERE/'registration.json') and run['results_sha256']==sha(HERE/'results.json')
    assert res['status']=='COMPLETE_SIX_CONFIGURATIONS_TRUE_SHORT_CAP_UPPER'
    rows=res['configurations'];assert [r['id']for r in rows]==list(range(1,7))
    assert run['native_calls']==run['spatial_grid_points']==run['matlab_calls']==0
    geo=read(HERE/'geometry_certificate.json');assert sha(HERE/'geometry_certificate.json')==rows[0]['geometry_certificate_sha256']
    v0=(1,)*6;other=geo['body_without_ones'];assert len(other)==31 and all(tuple(v)!=v0 for v in other)
    assert len(set(map(tuple,other)))==31 and all(v[0]==1 and all(x in(-1,1)for x in v)for v in other)
    line=[];pair=[];keys=set()
    for row in geo['lines']:
        v=other[row['index']];d=Q(6)-Q(dot(v,v0)**2,6)
        assert list(v)==row['direction'] and d==Q(row['distance_squared']);line.append(str(d))
    for row in geo['pairs']:
        i,j=row['indices'];assert 0<=i<j<31 and (i,j)not in keys;keys.add((i,j))
        v,w=other[i],other[j];a,b,s=dot(v,v0),dot(w,v0),dot(v,w)
        determinant=36-s*s
        d=Q(6)-Q(6*a*a-2*s*a*b+6*b*b,determinant)
        assert determinant==row['pair_Gram'] and d*determinant==row['triple_Gram']
        assert d==Q(row['distance_squared']);pair.append(str(d))
    assert len(line)==31 and len(keys)==465 and min(map(Q,line))==Q(10,3) and min(map(Q,pair))==2
    assert dict(Counter(line))==rows[0]['line_distance_histogram'] and dict(Counter(pair))==rows[0]['pair_distance_histogram']
    peaks=read(HERE.parents[1]/'round033/exact_projection_peaks/results.json');pc={r['id']:r for r in peaks['cases']}
    total=Q(0)
    for name,n in rows[1]['unordered_counts'].items():
        old=pc['r3_'+name];assert n==Q(29,32)*old['triple_orbit_count']
        total+=n*Q(old['R_exact'])
    assert sum(rows[1]['unordered_counts'].values())==4495 and rows[1]['ordered_total']==6*4495
    assert total==Q(rows[1]['sum_nprime_R_exact'])
    cn=24576*total;assert cn==Q(rows[1]['Crem_times_pi6_exact'])==Q(44743079,2970)
    # Reuse only the saved pi enclosure, not any old failed budget evaluation.
    olddir=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round036/exact_peak_budget')
    manifest=read(olddir/'manifest.json')
    assert sha(olddir/'results.json')==manifest['whitelist_sha256']['results.json']
    oldrows=read(olddir/'results.json')['configurations']
    piL,piU=bounds(oldrows[0]['pi']);assert Q(3)<piL<piU<Q(22,7)
    g=rows[2];rho=Q(g['rho']);eps=Q(g['epsilon']);assert rho==Q(1,64) and eps==Q(1,1000)
    assert g['time_denominator_coefficients']==[32*31,32*30,32*29]
    assert g['projected_frame_lower_eigenvalue']==26 and g['operator_upper_eigenvalue']==32
    etaL,etaU=bounds(g['eta']);assert 0<etaL and etaL**2<=6*rho*rho<=etaU**2
    gaps=[bounds(g[f'c{i}'])for i in range(3)];independent_gap_lowers=[]
    for (a,b),square in zip(gaps,[Q(3,2),Q(5,6),Q(1,2)]):
        lower=a*(1-Q(2)**-400)
        assert 0<lower<=a<=b and (lower+etaU)**2<=square
        independent_gap_lowers.append(lower)
    sqrt6U=etaU/rho
    volumeU=4*sqrt6U*piU*piU/45*rho**5
    # Exact endpoint inequalities verify the published norm upper bound without
    # recomputing any Arb scalar. A2^-384 outward enlargement absorbs differing
    # evaluation trees of the already stored512-bit enclosures.
    enlarge=1+Q(2)**-384
    goodL,goodU=bounds(rows[3]['good_branch_short_norm_upper']);goodU*=enlarge
    R1=Q(pc['r1']['R_exact'])
    assert goodU**6 >= Q(6,7)**6*R1/(8*piL**11)*eps**7
    lowerden=1
    for a in independent_gap_lowers:lowerden*=a*a
    remainder_densityU=cn/(piL**6*lowerden)
    badL,badU=bounds(rows[4]['bad_branch_short_norm_upper']);badU*=enlarge
    assert badU**6>=remainder_densityU**6*volumeU**5*eps**6
    simple=Q(rows[5]['simple_rational_upper'])
    assert goodU+badU<simple<Q(4,48750)
    final=bounds(rows[5]['short_cap_norm_upper']);assert final[1]<simple
    assert rows[4]['discarded_parameter_events']==0 and not rows[4]['Q2_used']
    assert rows[4]['ones_density_integrated_then_restored'] and not rows[5]['time_tail_above_epsilon_paid']
    return {'status':'PASS_SAVED_GEOMETRY_AND_RATIONAL_NORM_UPPER','checker_pid':os.getpid(),'elapsed_seconds':time.monotonic()-tic,'new_Arb_configurations':0,'new_kernel_calls':0,'saved_line_certificates_checked':31,'saved_pair_certificates_checked':465,'triple_counts_verified_by_transitive_incidence':True,'triple_enumeration_repeated':False,'simple_rational_short_cap_norm_upper':str(simple),'scope':'rho1/64,0<t<=1/1000,one positive-face S1 cap only','pi_source':str(olddir/'results.json'),'pi_source_sha256':sha(olddir/'results.json'),'no_old_budget_recomputed':True,'new_I_certificate':False,'new_C6_bound':False,'checker_source_sha256':sha(__file__)}
if __name__=='__main__':
    assert sys.argv[1:]in([],['--write-receipt'])
    result=check()
    if sys.argv[1:]==['--write-receipt']:
        assert not(HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
