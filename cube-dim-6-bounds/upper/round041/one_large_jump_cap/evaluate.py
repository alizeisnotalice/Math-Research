#!/usr/bin/env python3
"""Six fixed structural/scalar configurations; never evaluates a spatial kernel."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,combinations
from collections import Counter
import json,hashlib,time,datetime,os,traceback
from flint import arb,fmpq,fmpq_mat,ctx
ctx.prec=512
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,o):(HERE/name).write_text(json.dumps(o,indent=2)+'\n')
def A(x):
    x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def enclose(x):
    assert x.is_finite()
    a,e=x.lower().man_exp();b,f=x.upper().man_exp()
    lo=Q(int(a))*Q(2)**int(e);hi=Q(int(b))*Q(2)**int(f)
    assert 0<=hi-lo<=max(Q(1),abs((lo+hi)/2))*Q(2)**-450
    return {'lower_mantissa':str(a),'lower_binary_exponent':int(e),'upper_mantissa':str(b),'upper_binary_exponent':int(f),'display_only':str(x)}
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def det3(g):
    a,b,c=g;return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def run():
    assert not(HERE/'run_receipt.json').exists(),'one batch only'
    reg=read(HERE/'registration.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    assert sha(HERE/'evaluate.py')==reg['code_sha256']['evaluate.py']
    deadline=datetime.datetime.fromisoformat(reg['deadline_utc']).timestamp()
    assert time.time()<deadline
    receipt={'status':'RUNNING','pid':os.getpid(),'started_epoch':time.time(),'attempted_configurations':0,'completed_configurations':0,'registration_sha256':sha(HERE/'registration.json'),'native_calls':0,'spatial_grid_points':0,'matlab_calls':0}
    results={'status':'INCOMPLETE','precision_bits':512,'configurations':[],'new_I_certificate':False,'new_C6_bound':False}
    save('run_receipt.json',receipt)
    tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=6 and time.time()<deadline
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,data):
        results['configurations'].append({'id':i,**data});receipt['completed_configurations']=i
        save('results.json',results);save('run_receipt.json',receipt)
    try:
        begin(1)
        body=[(1,)+s for s in product((-1,1),repeat=5)];v0=(1,)*6;other=[v for v in body if v!=v0]
        line=[];pairs=[]
        for i,v in enumerate(other):line.append({'index':i,'direction':v,'distance_squared':str(Q(6)-Q(dot(v0,v)**2,6))})
        for i,j in combinations(range(31),2):
            v,w=other[i],other[j];s=dot(v,w);a=dot(v0,v);b=dot(v0,w)
            d2=36-s*s;d3=det3([[6,a,b],[a,6,s],[b,s,6]])
            distance=Q(d3,d2)
            # Independent FLINT 2x2 projection formula for the same geometry.
            gram=fmpq_mat([[6,s],[s,6]]);col=fmpq_mat([[a],[b]])
            projected=(col.transpose()*gram.inv()*col)[0,0]
            assert distance==Q(str(fmpq(6)-projected))>0
            pairs.append({'indices':[i,j],'pair_Gram':d2,'triple_Gram':d3,'distance_squared':str(distance)})
        assert len(line)==31 and len(pairs)==465
        assert min(Q(r['distance_squared'])for r in line)==Q(10,3)
        assert min(Q(r['distance_squared'])for r in pairs)==2
        save('geometry_certificate.json',{'body_without_ones':other,'lines':line,'pairs':pairs})
        done(1,{'name':'exact_nonones_line_and_pair_separation','line_count':31,'pair_count':465,'line_min_squared':'10/3','pair_min_squared':'2','line_distance_histogram':dict(Counter(x['distance_squared']for x in line)),'pair_distance_histogram':dict(Counter(x['distance_squared']for x in pairs)),'independent_fmpq_projection_check':True,'geometry_certificate_sha256':sha(HERE/'geometry_certificate.json')})
        begin(2)
        peaks=read(HERE.parents[1]/'round033/exact_projection_peaks/results.json');pc={r['id']:r for r in peaks['cases']}
        classes=[(1,1,1),(1,1,-1),(1,-1,1),(1,-1,-1)]
        names=['4110','3210','2220','3111','2211'];counts=Counter()
        for rows in combinations(other,3):
            columns=Counter(tuple(rows[k][j]*rows[0][j]for k in range(3))for j in range(6))
            name=''.join(map(str,sorted([columns[c]for c in classes],reverse=True)))
            assert name in names;counts[name]+=1
            assert det3([[dot(a,b)for b in rows]for a in rows])>0
        expected=dict(zip(names,[435,1740,435,580,1305]));assert dict(counts)==expected and sum(counts.values())==4495
        weighted=Q(0);oldweighted=Q(0)
        for name in names:
            oldn=pc['r3_'+name]['triple_orbit_count'];rr=Q(pc['r3_'+name]['R_exact'])
            assert counts[name]==Q(29,32)*oldn
            weighted+=counts[name]*rr;oldweighted+=oldn*rr
        assert weighted==Q(29,32)*oldweighted
        crem_numerical_coefficient=24576*weighted
        done(2,{'name':'all_nonones_triples','unordered_counts':expected,'unordered_total':4495,'ordered_total':26970,'sum_nprime_R_exact':str(weighted),'Crem_times_pi6_exact':str(crem_numerical_coefficient),'restricted_to_all_triples_ratio':'29/32','all_triples_independent':True})
        begin(3)
        pi=arb.pi();rho=A(Q(1,64));eps=A(Q(1,1000));eta=A(6).sqrt()*rho
        c0=A(6).sqrt()/2-eta;c1=A(Q(10,3)).sqrt()/2-eta;c2=A(Q(1,2)).sqrt()-eta
        assert c0>0 and c1>0 and c2>0
        vol=4*A(6).sqrt()*pi**2/45*rho**5
        done(3,{'name':'one_cap_and_three_positive_gaps','rho':'1/64','epsilon':'1/1000','threshold':'1/2','eta':enclose(eta),'c0':enclose(c0),'c1':enclose(c1),'c2':enclose(c2),'cap_volume':enclose(vol),'time_denominator_coefficients':[992,960,928],'projected_frame_lower_eigenvalue':26,'operator_upper_eigenvalue':32})
        begin(4)
        R1=Q(pc['r1']['R_exact']);Cbig=A(R1)/(8*pi**6)
        amplitude=Cbig**A(Q(1,6))*pi**A(Q(-5,6))
        good=A(Q(6,7))*amplitude*eps**A(Q(7,6))
        assert Q(R1,8)/(32*R1)==Q(1,256)
        done(4,{'name':'selected_large_parameter_short_cap_budget','R1_exact':str(R1),'Cbig':enclose(Cbig),'Cbig_over_old_wholeface_Cf_exact':'1/256','Lq_time_amplitude':enclose(amplitude),'good_branch_short_norm_upper':enclose(good),'new_quotient_peak_evaluations':0})
        begin(5)
        Crem=A(crem_numerical_coefficient)/pi**6
        M=Crem/(c0*c0*c1*c1*c2*c2)
        bad=M*vol**A(Q(5,6))*eps
        done(5,{'name':'all_complementary_parameter_events_paid','Crem':enclose(Crem),'uniform_bad_density_upper':enclose(M),'bad_branch_short_norm_upper':enclose(bad),'discarded_parameter_events':0,'Q2_used':False,'extra_time_parameter':False,'removed_directions':3,'ones_density_integrated_then_restored':True})
        begin(6)
        total=good+bad;target=A(Q(4,48750));ratio=total/target
        full=enclose(total);lo=Q(int(full['lower_mantissa']))*Q(2)**full['lower_binary_exponent'];hi=Q(int(full['upper_mantissa']))*Q(2)**full['upper_binary_exponent']
        simple=Q((hi*10**18).__ceil__(),10**18);assert hi<simple
        done(6,{'name':'one_fully_paid_finite_time_cap','short_cap_norm_upper':full,'simple_rational_upper':str(simple),'conditional_whole_face_Bcrit':'4/48750','cap_short_upper_over_Bcrit':enclose(ratio),'strictly_below_Bcrit':bool(total<target),'bad_branch_fraction_of_sum':enclose(bad/total),'entire_face_integral_claimed':False,'time_tail_above_epsilon_paid':False})
        results['status']='COMPLETE_SIX_CONFIGURATIONS_TRUE_SHORT_CAP_UPPER';save('results.json',results)
        receipt['status']='COMPLETED';code=0
    except BaseException as e:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(e);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-tic;receipt['finished_epoch']=time.time();receipt['no_active_children']=True
        if(HERE/'results.json').exists():receipt['results_sha256']=sha(HERE/'results.json')
        save('run_receipt.json',receipt)
    print(json.dumps({'status':receipt['status'],'configurations':receipt['completed_configurations'],'elapsed_seconds':receipt['elapsed_seconds'],'error':receipt.get('error')}))
    return code
if __name__=='__main__':
    assert sys.argv[1:]==['--run'];raise SystemExit(run())
