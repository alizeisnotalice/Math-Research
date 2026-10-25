#!/usr/bin/env python3
"""Independent exact reconstruction of saved B045 proposals and endpoints."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,importlib.util,itertools,json,os,time
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
PREP=HERE.parents[1]/'round015/spatial_minimum_search/prepare.py'
spec=importlib.util.spec_from_file_location('frozen_forms_only',PREP)
prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def iv(d):
    assert set(d)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
    lo,hi,e=int(d['lower_mantissa']),int(d['upper_mantissa']),int(d['common_binary_exponent'])
    assert lo<=hi and abs(e)<10000
    return F(lo)*F(2)**e,F(hi)*F(2)**e
def nearest(x):
    q,r=divmod(x.numerator,x.denominator)
    return q+(2*r>x.denominator or (2*r==x.denominator and q%2))
def key(n,d):return tuple(str(F(a,d))for a in n)

def check():
    start=time.monotonic();reg=read(HERE/'registration.json');run=read(HERE/'run_receipt.json');results=read(HERE/'results.json')
    assert reg['max_native']==4 and reg['g0']=='1/65000'
    assert reg['total_timeout_seconds']==600 and reg['per_native_timeout_seconds']==120
    assert sha(HERE/'run_batch.py')==reg['driver_sha256']
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    assert run['registration_sha256']==sha(HERE/'registration.json')
    assert run['results_sha256']==sha(HERE/'results.json')
    assert run['status']=='COMPLETED' and run['exit_code']==0 and run['all_child_handles_closed']
    assert run['elapsed_seconds']<600 and not run['seed_recomputed']
    assert len(run['native_calls'])==len(results['points'])<=4
    _,rays=prep.fan_data();forms=prep.all_zero_forms(rays)
    assert len(forms)==9308
    assert hashlib.sha256(json.dumps(forms,separators=(',',':')).encode()).hexdigest()==reg['zero_forms_sha256']
    seen=set(map(tuple,reg['old_point_registry']));assert len(seen)==reg['registry_final_count']==167
    den=10**14;gap=F(1,10**8);target=F(1,65000)
    nums=reg['seed_numerators'];fixed=nums[:4];oldpath=Path(reg['seed_native_path']);old=read(oldpath)
    assert key(nums,den)in seen and sha(oldpath)==reg['seed_native_sha256']
    rows=[];expected_status=None
    for i,(rec,saved)in enumerate(zip(run['native_calls'],results['points']),1):
        assert expected_status is None,'a call occurred after a stop condition'
        c=read(HERE/rec['candidate_path']);out=read(HERE/rec['stdout'])
        assert i==c['index']==rec['index']==saved['index']
        assert rec['status']=='COMPLETE_VALIDATED' and rec['exit_code']==0 and rec['elapsed_seconds']<120
        assert sha(HERE/rec['candidate_path'])==rec['candidate_sha256']
        assert sha(HERE/rec['stdout'])==rec['stdout_sha256']
        assert sha(HERE/rec['stderr'])==rec['stderr_sha256'] and (HERE/rec['stderr']).read_bytes()==b''
        assert c['previous_numerators']==nums and c['previous_native_path']==str(oldpath)
        assert c['previous_native_sha256']==sha(oldpath)
        grad=[sum(iv(v))/2 for v in old['gradient']]
        mag=max(abs(grad[4]),abs(grad[5]));assert mag>0
        direction=[F(0)]*4+[-grad[4]/mag,-grad[5]/mag]
        raw=[F(n,den)+d/20 for n,d in zip(nums,direction)]
        x5=min(F(fixed[3],den)-gap,max(2*gap,raw[4]))
        x6=min(x5-gap,max(gap,raw[5]))
        clipped=list(map(lambda n:F(n,den),fixed))+[x5,x6]
        quant=[nearest(x*den)for x in clipped]
        assert c['gradient_exact_midpoints']==list(map(str,grad))
        assert F(c['normalization_last_two'])==mag and c['direction']==list(map(str,direction))
        assert c['raw_proposal']==list(map(str,raw)) and c['clipped_proposal']==list(map(str,clipped))
        assert c['quantized_numerators']==quant
        assert list(map(F,c['quantization_errors']))==[F(n,den)-x for n,x in zip(quant,clipped)]
        pf=prep.preflight(quant,den,forms,rays);assert pf==c['initial_preflight']
        assert pf['feasible'] and pf['int64_safe']
        final=quant
        if pf['zero_form_count']:
            final=None
            for k,r in enumerate(c['repairs'],1):
                assert 1<=k<=16 and r['k']==k
                moved=quant[:];moved[5]+=k
                test=prep.preflight(moved,den,forms,rays)
                assert r['numerators']==moved and r['preflight']==test
                if test['feasible'] and test['int64_safe'] and test['zero_form_count']==0:
                    assert k==len(c['repairs']);final=moved;break
            assert final is not None
        else:assert c['repairs']==[]
        assert final[:4]==fixed and final==c['final_numerators']==rec['numerators']==saved['numerators']
        assert key(final,den)not in seen and c['not_in_old_registry']
        assert rec['denominator']==saved['denominator']==den
        assert rec['command']==[reg['binary_path'],'full',reg['fan_path'],str(den),*map(str,final)]
        pf=prep.preflight(final,den,forms,rays)
        assert pf==c['final_preflight'] and pf['feasible'] and pf['int64_safe'] and not pf['zero_form_count']
        assert F(c['actual_step_gradient_dot'])==sum(grad[j]*F(final[j]-nums[j],den)for j in (4,5))
        assert out['status']=='generic_nonzero_supported' and out['precision_bits']==512
        assert out['full_group_terms']==967680 and out['computed_sign_pairs']==483840
        assert out['zero_node_terms']==0 and out['confluent_terms']%2==0
        g=iv(out['G']);gg=[iv(v)for v in out['gradient']];assert len(gg)==6 and g[0]>0
        assert all(b-a<=max(F(1),abs((a+b)/2))*F(2)**-256 for a,b in [g,*gg])
        native_euler=iv(out['Euler_x_grad_plus5G']);assert native_euler[0]<=0<=native_euler[1]
        euler=(5*g[0]+sum(F(n,den)*a for n,(a,b)in zip(final,gg)),
               5*g[1]+sum(F(n,den)*b for n,(a,b)in zip(final,gg)))
        assert euler[0]<=0<=euler[1]
        oldg=iv(old['G']);diff=(g[0]-oldg[1],g[1]-oldg[0])
        assert saved['strict_descent']==(g[1]<oldg[0])
        assert saved['new_minus_previous_interval']==list(map(str,diff))
        v=saved['check']
        assert v['G_interval']==list(map(str,g)) and v['gradient_intervals']==[list(map(str,x))for x in gg]
        assert v['independent_Euler_interval']==list(map(str,euler))
        assert v['rejects_g0_strictly']==(g[1]<target)
        assert v['strictly_above_g0']==(g[0]>target)
        x=[F(n,den)for n in final]
        distance=max(min(abs(x[a]-x[b])for a,b in itertools.combinations(t,2))/2 for t in itertools.combinations(range(6),3))
        assert distance>0 and F(v['exact_dist_infinity_to_S2_lower_bound'])==distance
        bounds=(F((g[0]*10**20).__floor__(),10**20),F((g[1]*10**20).__ceil__(),10**20))
        assert bounds[0]<g[0]<=g[1]<bounds[1]
        rows.append({'index':i,'numerators':final,'denominator':den,'G_safe_rational_bounds':list(map(str,bounds)),
                     'strictly_above_target':g[0]>target,'strict_descent':g[1]<oldg[0],
                     'strict_increase':g[0]>oldg[1],'difference_interval':list(map(str,diff)),
                     'distance_to_S2_infinity_lower':str(distance),'repair_attempts':len(c['repairs']),
                     'min_abs_zeroform':pf['min_abs_primitive_form'],'last_two_gradient_intervals':[list(map(str,v))for v in gg[4:]]})
        if g[1]<target:expected_status='STRICT_G0_COUNTEREXAMPLE_IMMEDIATE_STOP'
        elif not g[1]<oldg[0]:expected_status='NO_CERTIFIED_STRICT_DESCENT_IMMEDIATE_STOP'
        elif i==4:expected_status='FOUR_STEPS_COMPLETE_GLOBAL_FLOOR_UNRESOLVED'
        seen.add(key(final,den));nums=final;old=out;oldpath=HERE/rec['stdout']
    assert results['status']==expected_status
    assert len(rows)==2 and rows[0]['strict_descent'] and rows[1]['strict_increase']
    assert all(r['strictly_above_target'] for r in rows)
    assert not(HERE/'candidate_03.json').exists() and not(HERE/'native_outputs/03.stdout.json').exists()
    return {'status':'PASS_EXACT_TWO_CALL_STOP_FIXED_STEP_RETURN','pid':os.getpid(),'wall_seconds':time.monotonic()-start,
            'new_native_calls':0,'original_native_calls':2,'no_seed_recompute':True,'no_extra_candidate':True,
            'source_hashes_unchanged':True,'points':rows,'target_refuted':False,'global_floor_proved':False,
            'local_or_global_minimum_proved':False,'new_C6_bound':False,'checker_sha256':sha(__file__)}

if __name__=='__main__':
    assert sys.argv[1:]in([],['--record'])
    r=check()
    if sys.argv[1:]==['--record']:
        assert not(HERE/'checked_results.json').exists()
        (HERE/'checked_results.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
