#!/usr/bin/env python3
"""One <=20s standard-library candidate audit; one saved worst row only.

No AS, LP, point evaluation, integral, all-row reconstruction or root search.
Stored integer-root inequalities are checked directly with integer powers.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
from pathlib import Path
import signal
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
DIRECTIONS_SHA='98c0e2fcc7f5d448aef33acc24f2b7fe22977f7a2361b20dc7c2ade926a0fd6e'
SCALING_SHA='1f7b1a0184101a93f5d64e463377ee4e86374660fee3e3c2de7676ff46237159'
COST_SOURCE_SHA='f5de232467e0060ffd3f0495b028cc6328372ef67c086a0a7693d5f38e762b67'
SOLVE_SHA='e9074068c755b7635f7fcdfb2acf262e6a1d48e49143b4490e5d581978c753cd'
BUILD_SHA='a03f2f13909d95b140505c2c870cbfe901c2755e93d3c1666de715a603501f8f'
CANDIDATE_SHA='7081cfef6fb47e5d47c95d3740c7a6ba48d48bb6197928a84d8cb6340c4166d5'
STRUCTURED_SHA='00c50260612ec3e8b15e810af2fe62269c10a09db0a7cd136c624df8a90bc02f'
COST_SHA='0120e7a55961ae7858541fb95319a1e4454c3c4e953c603478e51e488b7f8b33'
COLUMNS=[2,3,4,5,198,199,201,202,222,244,267,268]


def require(condition,message):
    if not condition:raise ValueError(message)


def read(path):return json.loads(path.read_text())


def nearest_even(value):
    n,d=value.numerator,value.denominator;q,r=divmod(n,d)
    return q+int(2*r>d or (2*r==d and q%2!=0))


def mult(key,N=80):
    i,j=map(int,key.split(','))
    require(key==f'{i},{j}' and 0<=j<=i<=N,'Noncanonical/out-of-window coefficient key.')
    return 1 if i==j==0 else 4 if j==0 or i==j else 8


def main():
    output=ROOT/'candidate_saved_verification.json'
    if output.exists():raise SystemExit('Existing audit: no repeated verification.')
    start=monotonic();answer={'status':'started','budget_seconds':20,'worst_rows_reconstructed':0,
        'new_AS':0,'new_integrals':0,'new_LP':0,'new_points':0,'integer_root_searches':0,
        'scope':'Saved binary solver vector to prescribed contraction/quantization, exact structure/mass/root bracket and one recorded worst row. Does not recheck all 36720 rows or cached-point feasibility.'}
    def alarm(signum,frame):raise TimeoutError('20s saved-candidate audit deadline')
    def guard():
        if monotonic()-start>=20:raise TimeoutError('20s saved-candidate audit deadline')
    def digest(path):
        h=sha256()
        with path.open('rb') as stream:
            while True:
                guard();data=stream.read(1<<20)
                if not data:break
                h.update(data)
        return h.hexdigest()
    def pinned(relative,wanted):
        p=ROOT/relative;require(digest(p)==wanted,'Source SHA mismatch: '+relative);return p
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,20)
    try:
        r=read(ROOT/'solve_result.json');run=read(ROOT/'solve.run.json');br=read(ROOT/'build_result.json');brun=read(ROOT/'build.run.json')
        require(r['status']=='one_finite_lp_and_exact_local_check_closed' and r['candidate_exported'] is True,'No exported final candidate.')
        require(run['status']=='terminal' and run['exit_code']==0 and run['budget_seconds']==180 and run['elapsed_seconds']<=180,'Solve run incomplete.')
        require(run['source_sha256_before']==run['source_sha256_after']==r['source_sha256']==SOLVE_SHA,'Solve source execution mismatch.')
        pinned('solve.py',SOLVE_SHA);require(digest(ROOT/'solve.log')==run['log_sha256'],'Solve log SHA mismatch.')
        require(digest(ROOT/'build_result.json')==r['build_result_sha256'] and br['complete'] is True,'Build-result binding mismatch.')
        require(brun['status']=='terminal' and brun['exit_code']==0 and brun['budget_seconds']==300 and brun['elapsed_seconds']<=300,'Build run incomplete.')
        require(brun['source_sha256_before']==brun['source_sha256_after']==br['source_sha256']==BUILD_SHA,'Build source execution mismatch.')
        pinned('build_directions.py',BUILD_SHA);require(digest(ROOT/'build.log')==brun['log_sha256'],'Build log SHA mismatch.')
        for relative,wanted in (('candidate.json',CANDIDATE_SHA),('structured_kernel.json',STRUCTURED_SHA),('exact_cost.json',COST_SHA),
                ('parameters.json',r['parameters_sha256']),('solver_result.json',r['solver_result_sha256']),
                ('solver.m',r['solver_source_sha256']),('matlab.log',r['matlab_log_sha256']),('model.local.mat',r['matrix_sha256'])):
            pinned(relative,wanted)
        require(r['candidate_sha256']==CANDIDATE_SHA and r['exact_cost_sha256']==COST_SHA,'Candidate/cost result binding mismatch.')
        base=read(pinned('../round01/scaled_certificate.json',BASELINE_SHA))
        directions=read(pinned('../round14/directions.json',DIRECTIONS_SHA))['directions']
        scaling=read(pinned('../round01/scaling_arithmetic.json',SCALING_SHA))
        pinned('../round10/enclose_cost.py',COST_SOURCE_SHA)
        p=read(ROOT/'parameters.json');candidate=read(ROOT/'candidate.json');structured=read(ROOT/'structured_kernel.json');cost=read(ROOT/'exact_cost.json')
        solver=read(ROOT/'solver_result.json')
        require(solver['input_sha256']==r['matrix_sha256'] and solver['status']=='returned' and solver['exitflag']==1 and
            solver['calls_started']==solver['calls_completed']==r['LP_started']==r['LP_completed']==r['LP_requested']==1,'Solver execution mismatch.')
        x=solver['x'];require(len(x)==26 and all(isinstance(v,(float,int)) and math.isfinite(v) for v in x),'Invalid solver vector.')
        require([d['old_column'] for d in directions]==COLUMNS and p['direction_names']==[d['name'] for d in directions],'Direction order mismatch.')
        c1,h,alpha,R=map(Q,(base['c'],base['h'],base['alpha'],base['R']))
        require((alpha,R,h,base['N'])==(Q(47,40),Q(33,20),Q(1,32),80),'Baseline family mismatch.')
        lam=Q(1023,1024);grid=1<<50
        # Explicit IEEE-double integer ratios, rather than decimal-string conversion.
        xf=[Q(*float(v).as_integer_ratio()) for v in x[:14]]
        contracted=[(1-lam)*c1+lam*xf[0],(1-lam)+lam*xf[1]]+[lam*v for v in xf[2:]]
        g0=contracted[0]-contracted[1]*c1;g=Q(nearest_even(g0*grid),grid)
        s=Q(nearest_even(contracted[1]*grid),grid)
        expected=[s*c1+g,s]+[Q(nearest_even(v*grid),grid) for v in contracted[2:]]
        expected+=list(map(abs,expected[2:14]));q=list(map(Q,p['q']))
        require(q==expected and len(q)==26,'Prescribed contraction/unique quantization mismatch.')
        require(Q(p['lambda'])==Q(r['fixed_lambda'])==lam and Q(p['g'])==Q(r['structural_g'])==g and
            Q(r['structural_g_before_rounding'])==g0,'Saved contraction/g provenance mismatch.')
        require(g>=0 and q[0]-q[1]*c1==g and c1/2<=q[0]<=3*c1/2 and Q(1,2)<=q[1]<=Q(3,2),'Structural or c/s bounds failed.')
        for j,d in enumerate(directions):
            bound=Q(d['bound']);require(abs(q[2+j])<=bound and q[14+j]==abs(q[2+j]) and 0<=q[14+j]<=bound,'Signed/auxiliary bounds failed.')
        require(set(candidate)=={'alpha','R','h','N','c','a'} and all(candidate[k]==structured[k] for k in candidate),'Plain/structured schema mismatch.')
        require(all(candidate[k]==base[k] for k in ('alpha','R','h','N')) and Q(structured['beta'])==alpha and Q(candidate['c'])==q[0],
            'Candidate family/c changed.')
        st=structured['structure'];require(st['baseline_sha256']==BASELINE_SHA and Q(st['s'])==q[1] and Q(st['g'])==g,'Structure parameter mismatch.')
        T={};masses=[]
        for j,d in enumerate(directions):
            guard();dm=Q(0)
            for key,value in d['coefficients'].items():
                v=Q(value);dm+=h*h*mult(key)*v/2
                require(not v or h*(sum(map(int,key.split(',')))+4)<=Q(3,4),'Direction support exceeds 3/4.')
                T[key]=T.get(key,Q(0))+q[2+j]*v
            require(dm==Q(d['half_mass']),'Direction exact D4 mass mismatch.');masses.append(dm)
        T={k:v for k,v in T.items() if v}
        require({k:Q(v) for k,v in st['T'].items()}==T and st['T_support_bound']=='3/4','Expanded T structure mismatch.')
        a={key:q[1]*Q(v) for key,v in base['a'].items()}
        for key,value in T.items():a[key]=a.get(key,Q(0))+value
        a={k:v for k,v in a.items() if v};saved_a={k:Q(v) for k,v in candidate['a'].items()}
        require(a==saved_a and a=={k:Q(v) for k,v in structured['a'].items()},'Exact s*a1+T expansion mismatch.')
        full_sum=sum((mult(k)*v for k,v in a.items()),Q(0));offset=h*h*full_sum/2
        smass=h*h*sum((mult(k)*Q(v) for k,v in base['a'].items()),Q(0))/2
        require(offset==q[1]*smass+sum((q[2+j]*masses[j] for j in range(12)),Q(0)),'D4 full/affine cost mismatch.')
        require(cost['certificate_sha256']==STRUCTURED_SHA and cost['cost_source_sha256']==COST_SOURCE_SHA and
            Q(cost['alpha'])==Q(cost['beta'])==alpha and Q(cost['R'])==R,'Cost source/family binding mismatch.')
        require(cost['root_degree']==40 and cost['root_power']==33,'Wrong stored radical.')
        scale,root=int(cost['root_scale']),int(cost['root_floor'])
        require(scale==10**60 and root>=0,'Wrong root scale/floor.')
        rad_n,rad_d=R.numerator,R.denominator
        require(root**40*rad_d**33<=rad_n**33*scale**40<(root+1)**40*rad_d**33,'Stored 40th-root bracket failed.')
        factor=q[0]*alpha/(2-alpha)
        lower=factor*Q(root,scale)+offset;upper=factor*Q(root+1,scale)+offset
        require(Q(cost['sum_coefficients'])==full_sum and Q(cost['cost_factor'])==factor and Q(cost['cost_offset'])==offset,
            'Stored exact cost factors mismatch.')
        require(Q(cost['lower'])==lower and Q(cost['upper'])==upper and Q(cost['width'])==upper-lower,'Stored cost enclosure mismatch.')
        b1lower=Q(scaling['new_cost_lower'])
        require(Q(cost['baseline_lower'])==b1lower and Q(cost['strict_gain_lower'])==b1lower-upper,'Baseline/gain enclosure mismatch.')
        require(cost['improves_baseline']==(upper<b1lower) and cost['below_user_target']==(upper<Q(17,5)) and
            cost['user_target']=='17/5' and upper<b1lower,'Cost target flags mismatch.')
        worst=r['worst_exact_row']
        require(worst=={'node_id':658,'region_key':[658,1],'parent':[52,1],'bernstein_index':[2,0]},'Frozen worst-row identity changed.')
        row_path=pinned('rows_new.json',r['new_rows_sha256'])
        require(br['new_rows_sha256']==br['rows_new_sha256']==r['new_rows_sha256'],'Build/solve new-row SHA mismatch.')
        matches=[];scanned=0
        # The pinned streaming writer emits one compact complete row per line.
        with row_path.open() as stream:
            for raw in stream:
                guard();scanned+=1
                if '"leaf_node_id":658,' in raw and '"region_key":[658,1]' in raw and '"bernstein_index":[2,0]' in raw:
                    matches.append(json.loads(raw.strip().rstrip(',')))
        require(len(matches)==1,'Saved worst row not uniquely found.');row=matches[0]
        require(len(row['t_midpoints'])==len(row['t_penalties'])==12,'Worst row direction count mismatch.')
        residual=q[0]*Q(row['c_column'])+q[1]*Q(row['s_column'])+sum((q[2+j]*Q(row['t_midpoints'][j])-
            q[14+j]*Q(row['t_penalties'][j]) for j in range(12)),Q(0))-Q(row['baseline_slack'])/2
        require(residual==Q(r['exact_min_slack']) and residual>=0,'One recorded worst-row residual mismatch.')
        require(r['exact_rows_checked']==36720 and r['exact_bad_count']==0 and r['point_rows_checked']==103148 and r['point_bad_count']==0,
            'Recorded overall gate summary incomplete.')
        for batch in ('float_point_checks','point_checks'):
            require(all(r[batch][kind]['strengthened_bad']==r[batch][kind]['original_bad']==0 for kind in ('K','A')),
                'Recorded cache gate failed.')
        guard();answer.update(status='passed_saved_candidate_arithmetic_audit',worst_rows_reconstructed=1,
            worst_row_identity=worst,worst_row_residual=str(residual),source_lines_scanned=scanned,
            quantized_parameters_checked=26,expanded_nonzero_coefficients=len(a),stored_root_bracket_checks=1,
            structural_g=str(g),cost_lower=str(lower),cost_upper=str(upper),strict_gain_lower=str(b1lower-upper),
            improves_baseline=True,below_user_target=upper<Q(17,5),claim_359087_over_100000=upper<Q(359087,100000),
            candidate_sha256=CANDIDATE_SHA,structured_sha256=STRUCTURED_SHA,exact_cost_sha256=COST_SHA,
            solve_result_sha256=digest(ROOT/'solve_result.json'),solve_run_sha256=digest(ROOT/'solve.run.json'),
            build_result_sha256=digest(ROOT/'build_result.json'),build_run_sha256=digest(ROOT/'build.run.json'),
            parameters_sha256=digest(ROOT/'parameters.json'),solver_result_sha256=digest(ROOT/'solver_result.json'),
            solve_source_sha256=SOLVE_SHA,build_source_sha256=BUILD_SHA,wrapper_sha256=digest(ROOT/'run_stage.py'))
    except Exception as exc:answer.update(status='failed',reason=f'{type(exc).__name__}: {exc}')
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        answer.update(elapsed_seconds=monotonic()-start,verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
        tmp=output.with_suffix('.json.tmp');tmp.write_text(json.dumps(answer,indent=2,allow_nan=False)+'\n');tmp.replace(output)
    print(json.dumps(answer));return 0 if answer['status'].startswith('passed_') else 1


if __name__=='__main__':sys.exit(main())
