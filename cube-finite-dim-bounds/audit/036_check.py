"""Independent record acceptance, no author imports and no Green node evaluation.
Default --main is project root containing audit/ and upper/. Optional raw audit
checks lossless compaction against external immutable original records.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, time, platform
from flint import arb, ctx
HERE=Path(__file__).resolve().parent

def load(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def aq(x):
    x=F(x); return arb(x.numerator)/x.denominator

def unpack(v):
    assert len(v)==4 and int(v[2])==0
    return F(int(v[0]))*F(2)**int(v[1])
def endpoint(x, upper=True): return F(str((x.upper() if upper else x.lower()).fmpq()))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--main',type=Path,default=HERE.parent)
    ap.add_argument('--raw-archive',type=Path);ap.add_argument('--save',type=Path,default=HERE/'036_results.json')
    args=ap.parse_args();root=args.main.resolve();t=time.perf_counter();ctx.prec=512
    snapshot=load(HERE/'036_source_snapshot.json')
    for rel,h in snapshot.items(): assert sha(root/rel)==h, rel
    b=root/'upper/continuous_bridge';s=b/'cb11_integral_cheb_final_staging'
    manifest=load(b/'cb11_integral_cheb_final_manifest.json')
    assert len(manifest['payload'])==25
    # Source snapshot already locks each exact payload; no executable from it is run.
    expected_config={'schema':1,'cone_sha256':'23258e3d0236700ff442d602ccbb738a17c92728f3b40d95ddf9dd5d9ef924f7','epsilon':'9/5000','tau':'1/3','degree':32,'grid':128,'node_bits':128,'precision_bits':384,'rho':'5/2'}
    checker='4f89a377fbdf263061c2d6525d2101a773f7a6ba4127bb367bd3b7328ed57eb5'
    total=F(0);weight_sum=0;seen=set();raw_hash={};cached=[];fresh=0
    provenance=load(b/'cb11_integral_cheb_final_provenance.json')
    for part in range(12):
        name=f'cb11_integral_cheb_part_{part:02}.json';d=load(s/name)
        assert d['config']==expected_config and d['checker_sha256']==checker
        assert d['part']==part and d['parts']==12 and len(d['patches'])==688
        wanted={(i,j) for i in range(128) for j in range(i+1) if (i*(i+1)//2+j)%12==part}
        assert {tuple(r['patch']) for r in d['patches']}==wanted
        raw=None
        if args.raw_archive:
            p=args.raw_archive/name;h=sha(p);assert h==provenance['raw_parts'][name]['sha256'];raw_hash[name]=h
            raw=load(p);assert raw['config']==expected_config and raw['checker_sha256']==checker
            assert raw['part']==part and raw['parts']==12
            assert len(raw['patches'])==688
        for k,r in enumerate(d['patches']):
            ij=tuple(r['patch']);assert ij not in seen;seen.add(ij);i,j=ij
            L,U,m,c=map(unpack,(r['L'],r['U'],r['mean_H_upper'],r['IH_contribution_upper']))
            assert 0<=L<=m<=U and U>0 and c==F(r['IH_contribution_exact'])
            w=1 if i==j else 2;weight_sum+=w
            # Evaluate the secant at the rigorous mean using positive square roots,
            # independently at 512 bits, from the exact dyadic recorded endpoints.
            l,u=aq(L).sqrt(),aq(U).sqrt()
            chord=l**3+(u*u+u*l+l*l)/(u+l)*aq(m-L)
            own=chord*w/128**2
            assert endpoint(own)<=c, ('chord',ij)
            total+=c
            if raw is not None:
                rr=raw['patches'][k]
                assert all(rr[key]==value for key,value in r.items()),ij
                assert rr['nodes']==1089 and rr['weight']==w
                assert list(map(F,rr['whole_square']))==[F(i,128),F(i+1,128),F(j,128),F(j+1,128)]
                assert unpack(rr['chord_mean_upper'])*w/128**2==c
                if rr['cache_hit']: cached.append(list(ij))
                else: fresh+=rr['nodes']
    assert len(seen)==8256 and weight_sum==16384 and F(weight_sum,128**2)==1
    saved=load(s/'cb11_integral_cheb_all_result.json')
    assert saved['config']==expected_config and saved['patch_count']==8256
    assert F(saved['IH_upper_exact'])==total==F(provenance['IH_upper_exact'])
    # The L^(3/2) whole-plane short-time bound, then restriction to [0,1]^2.
    eps=aq(F(9,5000));pi=arb.pi()
    et=(aq(F(1,16)).root(3)*aq(F(3,4))*(aq(F(256,3))/pi**3).root(3)
        *(aq(F(16,3))/pi).root(3)**2*eps.root(3)**4)
    packed=unpack(saved['I_upper']);assert saved['I_upper']==manifest['I_upper_packed']
    combined=(aq(total).root(3)**2+et)
    full=combined*combined.sqrt()
    assert endpoint(full)<=packed
    I0=F(23129127,40000000000);g0=F(2621323989,500000000000);U0=F(1979,250)
    gcert=load(root/'audit/033_results.json')
    assert gcert['patch_count']==136 and len(gcert['patches'])==136
    glow=min(F(r['H_lower_exact']) for r in gcert['patches'])
    assert glow==F(gcert['minimum_H_lower_exact']) and 0<g0<=glow
    assert packed<=I0
    cost2=27*I0**2/g0**3;gap=U0**2-cost2;assert gap>0
    old=load(root/'results/round011_n3_full_certificate_cost.json')
    assert cost2==F(old['cost_squared_exact'])
    assert g0==F(old['g_outward_simple']) and I0==F(old['I_outward_simple'])
    # Exact independent Machin arithmetic; the angle branch proof is in report.
    def atan_sum(x,m): return sum(((-1)**k*x**(2*k+1)/F(2*k+1) for k in range(m)),F(0))
    mach=16*atan_sum(F(1,5),13)-4*atan_sum(F(1,239),4)
    assert mach<F(355,113)
    assert (F(120,119)-F(1,239))/(1+F(120,119*239))==1
    if args.raw_archive:
        assert sorted(cached)==sorted([[127,57],[127,127],[127,0],[0,0]])
        assert fresh==8986428 and fresh+len(cached)*1089==8990784
    result={'status':'PASS C3 < 1979/250 = 7.916, complete fixed candidate',
      'source_sha256':snapshot,'precision_bits':512,'python':platform.python_version(),
      'patches':len(seen),'parts':12,'patches_per_part':688,'weighted_area_exact':'1',
      'independent_chords_checked':len(seen),'IH_upper_exact':str(total),
      'Etime_interval':[str(endpoint(et,False)),str(endpoint(et))],
      'combined_I_interval':[str(endpoint(full,False)),str(endpoint(full))],
      'packed_I_exact':str(packed),'I_outward':str(I0),'g033_full_lower':str(glow),'g_outward':str(g0),
      'cost_squared_exact':str(cost2),'cost_interval':str(aq(cost2).sqrt()),
      'strict_square_gap':str(gap),'accepted_U3':str(U0),
      'machin_pi_upper':str(mach),'machin_gap_to_355_113':str(F(355,113)-mach),
      'raw_records_checked':8256 if args.raw_archive else 0,'raw_sha256':raw_hash,
      'cached_patch_indices':cached,'original_fresh_nodes_verified':fresh if args.raw_archive else None,
      'original_total_nodes_verified':8990784 if args.raw_archive else None,
      'new_Green_node_evaluations':0,'shared_arithmetic':'python-flint Arb; no author program imported or executed',
      'wall_seconds':time.perf_counter()-t}
    args.save.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','patches','weighted_area_exact','cost_interval','new_Green_node_evaluations','wall_seconds']},indent=2))
if __name__=='__main__':main()
