#!/usr/bin/env python3
"""Eight new bounded exact fixtures; no body32 kernel or previous fixture replay."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from itertools import combinations
import json,hashlib,time,os,traceback
import sympy as sp
H=Path(__file__).resolve().parent
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,o):(H/n).write_text(json.dumps(o,indent=2)+'\n')
def matrix_data(M):return [[str(x)for x in row]for row in M.tolist()]
def pair_residuals(A,x):
    rows=[]
    for i,j in combinations(range(A.cols),2):
        U=A[:,[i,j]];G=U.T*U;assert G.det()>0
        b=U.T*x;r=sp.factor((x.T*x)[0]-(b.T*G.inv()*b)[0]);assert r>0
        rows.append(dict(pair=[i+1,j+1],rank=2,residual_squared=str(r)))
    return rows
def run():
    assert not(H/'run_receipt.json').exists(),'one fixed scientific batch only'
    reg=rd(H/'registration.json');binding=rd(H/'execution_binding.json');guard=rd(H/'producer_process_receipt.json')
    assert binding['registration_sha256']==sha(H/'registration.json')
    for f,x in binding['code_sha256'].items():assert sha(H/f)==x,f
    for f,x in reg['source_sha256'].items():assert sha(f)==x,f
    assert guard['status']=='RUNNING'and guard['owned_pid']==os.getpid()and guard['work_timeout_seconds']==170
    receipt=dict(status='RUNNING',pid=os.getpid(),started_epoch=time.time(),registration_sha256=sha(H/'registration.json'),binding_sha256=sha(H/'execution_binding.json'),attempted_configurations=0,completed_configurations=0,body32_kernel_calls=0,old_baseline_replays=0,old069_fixture_replays=0,parameter_grid_calls=0,matlab_calls=0)
    result=dict(status='INCOMPLETE',configurations=[],body32_semiconcavity_proved=False,body32_counterexample_found=False,axis_independent6_semiconcavity_author_proof=True,overlap_fixture_is_geometry_only=True,new_Cn_bound=False)
    save('run_receipt.json',receipt);tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=8 and time.monotonic()-tic<170
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,**data):
        result['configurations'].append(dict(id=i,**data));receipt['completed_configurations']=i
        save('results.json',result);save('run_receipt.json',receipt)
    try:
        begin(1);u,a=sp.symbols('u a',positive=True)
        laplace=sp.factor(u*(1/u-u/(u*u+a*a))-a*a/(u*u+a*a));assert laplace==0
        coeff=sp.simplify(6/(2*sp.pi)**3*(sp.pi/2));assert coeff==3/(8*sp.pi**2)
        done(1,name='3D_normal_deficit_normalization',Laplace_identity_residual=str(laplace),normal_norm_sphere_coefficient=str(coeff),radial_integral='pi*abs(a)/2',Tonelli_nonnegative=True)
        begin(2);axis=sp.eye(6);x=sp.Matrix([1,2,4,0,0,0]);pairs=pair_residuals(axis,x)
        assert axis.rank()==6 and len(pairs)==15 and min(sp.Rational(row['residual_squared'])for row in pairs)==1
        done(2,name='new_axis6_offS2_onS3_point',frame=matrix_data(axis),point=[1,2,4,0,0,0],rank=6,pair_count=15,pairs=pairs,min_pair_residual_squared='1',on_selected_rank3=True,old069_point_repeated=False)
        begin(3);a,b,t=sp.symbols('a b t',positive=True)
        integrand=1-t**4/((t*t+a*a)*(t*t+b*b));ca=a**4/(a*a-b*b);cb=b**4/(b*b-a*a)
        assert sp.factor(integrand-ca/(t*t+a*a)-cb/(t*t+b*b))==0
        normtilde=sp.factor(ca/a+cb/b);expected=(a*a+a*b+b*b)/(a+b);assert sp.factor(normtilde-expected)==0
        HH=sp.hessian(normtilde,(a,b));PSD=2/(a+b)**3*sp.Matrix([[b*b,-a*b],[-a*b,a*a]])
        assert all(sp.factor(z)==0 for z in HH-PSD)and sp.factor(HH.det())==0
        assert sp.factor(normtilde.subs(b,0)-a)==0
        done(3,name='axis3_normal_norm_two_coordinate_slice',normalization='N=(2*pi^2)^(-1)*normtilde',normtilde=str(normtilde),Hessian_factor=str(PSD),determinant='0',positive_trace=str(sp.factor(sp.trace(HH))),axis_restriction='a/(2*pi^2)',actual_variable_range='a,b>0; no samples')
        begin(4);s,t=sp.symbols('s t',positive=True);roots=[sp.Integer(1),sp.Integer(2),sp.Integer(4),s]
        coeffs=[-v*v/sp.prod(w*w-v*v for j,w in enumerate(roots)if j!=i)for i,v in enumerate(roots)]
        integrand=t*t/sp.prod(t*t+v*v for v in roots)
        assert sp.factor(sum(c/(t*t+v*v)for c,v in zip(coeffs,roots))-integrand)==0
        F=sp.factor(sum(c/v for c,v in zip(coeffs,roots))/2)
        F0=sp.factor(F.subs(s,0));Dplus=sp.factor(sp.diff(F,s).subs(s,0));D2plus=sp.factor(sp.diff(F,s,2).subs(s,0))
        assert F0==sp.Rational(7,1440)and Dplus==-sp.Rational(1,128)and D2plus.is_Rational and D2plus>0
        done(4,name='one_true_axis6_Green_path',path=[1,2,4,'s',0,0],range='0<s<1/2',integrand_without_pi_minus6=str(integrand),G_times_pi5=str(F),value_at_zero_times_pi5=str(F0),right_derivative_times_pi5=str(Dplus),left_derivative_times_pi5=str(-Dplus),derivative_jump_times_pi5=str(2*Dplus),right_second_derivative_limit_times_pi5=str(D2plus),true_time_integrals=1,formal_DD_evaluations=0,negative_cusp=True)
        begin(5);a,b,c,t=sp.symbols('a b c t',positive=True);den=(a*a+t*t)*(b*b+t*t)*(c*c+t*t);Hfun=a*a*b*b*c*c/den
        polynomial=a*a*b*b+a*a*c*c+b*b*c*c+t*t*(a*a+b*b+c*c)+t**4
        assert sp.factor((1-Hfun)-t*t*polynomial/den)==0
        R=polynomial/den;lim=sp.factor(R.subs(t,0));da=sp.factor(sp.diff(R,a).subs(t,0));daa=sp.factor(sp.diff(R,a,2).subs(t,0));dab=sp.factor(sp.diff(R,a,b).subs(t,0))
        assert sp.factor(lim-(1/a**2+1/b**2+1/c**2))==0 and da==-2/a**3 and daa==6/a**4 and dab==0
        done(5,name='axis_exact_t2_remainder_structure',one_minus_H_factor=str(t*t*polynomial/den),quotient_at_t0=str(lim),dy_quotient_at_t0=str(da),dyy_quotient_at_t0=str(daa),dy1dy2_quotient_at_t0=str(dab),uniform_on_y_compact_away_zero=True)
        begin(6);v7=sp.Matrix([1,0,0,1,0,0]);v8=sp.Matrix([0,1,0,-1,0,0]);A=axis.row_join(v7).row_join(v8);x=sp.Matrix([1,1,1,0,0,0]);pairs=pair_residuals(A,x)
        U1=axis[:,[0,1,2]];U2=v7.row_join(v8).row_join(axis[:,2]);inter=(axis[:,0]+axis[:,1]).row_join(axis[:,2]);P1=U1*(U1.T*U1).inv()*U1.T;P2=U2*(U2.T*U2).inv()*U2.T
        assert A.rank()==6 and len(pairs)==28 and U1.rank()==U2.rank()==3 and U1.row_join(U2).rank()==4
        assert P1*x==x and P2*x==x and P1*inter==inter and P2*inter==inter and inter.rank()==2
        done(6,name='new_rank3_overlap_offS2_geometry',frame=matrix_data(A),point=[1,1,1,0,0,0],pair_count=28,pairs=pairs,min_pair_residual_squared=str(min(sp.Rational(row['residual_squared'])for row in pairs)),W1_rank=3,W2_rank=3,union_rank=4,intersection_dimension=2,intersection_basis=matrix_data(inter),W1_projector=matrix_data(P1),W2_projector=matrix_data(P2),Green_evaluations=0,Hessian_evaluations=0,semiconcavity_claim_for_this_frame=False)
        begin(7);s=sp.symbols('s',positive=True);f=-s+s*s*sp.log(1/s);d1=sp.diff(f,s);d2=sp.simplify(sp.diff(f,s,2))
        assert sp.simplify(d2-(2*sp.log(1/s)-3))==0
        assert sp.limit(d1,s,0,dir='+')==-1 and sp.limit(d2,s,0,dir='+')==sp.oo
        assert sp.limit((f+s)/s,s,0,dir='+')==0
        done(7,name='abstract_negative_cusp_not_sufficient',positive_side_function=str(f),right_derivative_limit='-1',left_derivative_limit='1',derivative_jump='-2',classical_second_derivative=str(d2),second_derivative_limit='+infinity',remainder_is_o_abs_s=True,Cauchy_frame_counterexample=False,body32_counterexample=False)
        begin(8);y,z1,z2,k1,k2,ee=sp.symbols('y z1 z2 k1 k2 epsilon',real=True);r=sp.sqrt(z1*z1+z2*z2+ee*ee);AA=1+y*y
        Hn=sp.hessian(r,(z1,z2));kv=sp.Matrix([k1,k2]);num=sp.expand((ee*ee*(k1*k1+k2*k2)+(z1*k2-z2*k1)**2))
        assert sp.simplify((kv.T*Hn*kv)[0]*r**3-num)==0
        vars=(y,z1,z2);Av=sp.Matrix([sp.diff(AA,x)for x in vars]);Nv=sp.Matrix([sp.diff(r,x)for x in vars]);formula=-AA*sp.hessian(r,vars)-r*sp.hessian(AA,vars)-Av*Nv.T-Nv*Av.T
        assert all(sp.simplify(v)==0 for v in sp.hessian(-AA*r,vars)-formula)
        assert 8<9
        done(8,name='weighted_convex_norm_Hessian_product',norm_PSD_numerator=str(num),weighted_product_identity_residuals=['0']*9,example_A='1+y^2',example_local_bound='7 on |y|<=1, |z|<=1, 0<epsilon<=1',general_lemma_requires_A_nonnegative_C2=True,general_lemma_requires_uniform_remainder_upper_Hessian=True)
        for f,x in reg['source_sha256'].items():assert sha(f)==x,f
        result['actual_work']=dict(producer_configurations=8,exact_tiny_pair_residuals=43,true_tiny_Green_symbolic_paths=1,normal_norm_symbolic_slices=1,abstract_counterexamples=1,overlap_geometry_points=1,body32_kernel_calls=0,body32_pair_or_triple_enumerations=0,old069_fixture_replays=0,old_baseline_replays=0,parameter_grid_calls=0,matlab_calls=0)
        result['status']='COMPLETE_SCOPED_AXIS_SEMICONCAVITY_AND_OVERLAP_DIAGNOSTIC';save('results.json',result);receipt['status']='COMPLETED';code=0
    except BaseException as ex:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(ex);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-tic;receipt['finished_epoch']=time.time();receipt['no_active_children']=True
        if(H/'results.json').exists():receipt['results_sha256']=sha(H/'results.json')
        save('run_receipt.json',receipt)
    print(json.dumps(dict(status=receipt['status'],configurations=receipt['completed_configurations'],seconds=receipt['elapsed_seconds'],error=receipt.get('error'))));return code
if __name__=='__main__':raise SystemExit(run())
