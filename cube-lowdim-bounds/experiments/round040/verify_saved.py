#!/usr/bin/env python3
"""Independent standard-library saved-value audit; <=90s, no integration/Arb.

This rebuilds Fraction arithmetic and execution provenance, not the returned
Hab enclosures or their analytic validity. It is not a second integral method.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SOURCE_SHA='0c6531e82741ec785b4a0f66ba84a4630871c73cf84d9ba797858d589af2babf'
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
FIXED_SHA='4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05'
MODEL_SHA='8dd22932cefd30806acf819ffae37f10d9497e50e2e798361e10f0d9e207db26'
COVER_SHA='6028e48efe444e27c5259858cf1eeb8b41e0081b6e39f06102439ddb198c4cac'
PHYSICAL=['6677/4096','3339/2048','5/128','161/4096']
LOCAL=['21/128','11/64','1/4','33/128']


def require(condition,message):
    if not condition:raise ValueError(message)


def digest(path):return sha256(path.read_bytes()).hexdigest()


def read(path):return json.loads(path.read_text())


def main():
    output=ROOT/'saved_verification.json'
    if output.exists():raise SystemExit('Existing verification: no repeated audit.')
    start=monotonic()
    answer={'status':'started','budget_seconds':90,'rows_checked':0,
        'scope':'Exact saved Hab interval endpoint inequalities, one-sided correction, reused coefficients, row arithmetic and execution ledger only; no independent integration.',
        'new_integrals':0,'new_model_evaluations':0,'new_points':0,'LP_calls':0}
    def alarm(signum,frame):raise TimeoutError('90s saved arithmetic audit deadline')
    def guard():
        if monotonic()-start>=90:raise TimeoutError('90s saved arithmetic audit deadline')
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,90)
    try:
        j=read(ROOT/'result.json');m=read(ROOT/'signed_model.json');run=read(ROOT/'signed.run.json')
        require(run['status']=='terminal' and run['exit_code']==0,'Run did not finish successfully.')
        require(run['budget_seconds']==60 and 0<=run['elapsed_seconds']<=60,'Run budget mismatch.')
        require(run['command_relative']==['python','-B','signed_once.py'],'Unexpected executed command.')
        require(isinstance(run['pid'],int) and run['pid']>0,'Missing run PID.')
        source_hash=digest(ROOT/'signed_once.py')
        require(run['source_sha256_before']==run['source_sha256_after']==j['source_sha256']==source_hash==SOURCE_SHA,
            'Execution source binding failed.')
        require(run['log_sha256']==digest(ROOT/'signed.log'),'Run log hash mismatch.')
        require(j['status']=='completed' and j['in_flight']==[],'Numeric task did not complete.')
        require(0<=j['elapsed_seconds']<=run['elapsed_seconds'],'Internal elapsed time exceeds outer run.')
        require(j['signed_model_sha256']==digest(ROOT/'signed_model.json'),'Signed model hash mismatch.')
        require(j['calls_sha256']==digest(ROOT/'calls.jsonl'),'Call trace hash mismatch.')
        expected_inputs={
            'round17/fixed_checker/hlcert.py':FIXED_SHA,
            'round20/models.local/node_00661.json':MODEL_SHA,
            'round20/leaf_cover.json':COVER_SHA,
            'round01/scaled_certificate.json':BASELINE_SHA}
        require(len(j['inputs'])==4 and {v['file']:v['sha256'] for v in j['inputs']}==expected_inputs,
            'Input manifest mismatch.')
        for relative,wanted in expected_inputs.items():
            guard();require(digest(ROOT.parent/relative)==wanted,'Input hash mismatch: '+relative)
        old=read(ROOT.parent/'round20/models.local/node_00661.json')
        cover=read(ROOT.parent/'round20/leaf_cover.json');base=read(ROOT.parent/'round01/scaled_certificate.json')
        node=cover['nodes'][661]
        require((node['id'],node['parent'],node['depth'],node['model_sha256'])==(661,[52,1],4,MODEL_SHA),
            'Source cover node identity mismatch.')
        require(node['local']==old['local']==m['local']==j['local']==LOCAL and
            node['physical']==old['physical']==m['physical']==j['physical']==PHYSICAL,'Geometry mismatch.')
        require(old['node_id']==m['node_id']==j['node_id']==661,'Wrong fixed node.')
        require(old['baseline_sha256']==m['baseline_sha256']==BASELINE_SHA and
            old['fixed_hlcert_sha256']==m['fixed_hlcert_sha256']==FIXED_SHA,'Model source binding mismatch.')
        require(m['source_model_sha256']==MODEL_SHA and m['source_cover_sha256']==COVER_SHA,'Reused source hashes mismatch.')
        require(Q(base['alpha'])==Q(47,40) and Q(base['R'])==Q(33,20) and Q(base['h'])==Q(1,32) and base['N']==80,
            'Baseline family mismatch.')
        c1,h,R=Q(base['c']),Q(base['h']),Q(base['R']);u0,u1,v0,v1=map(Q,PHYSICAL);x0,x1,y0,y1=map(Q,LOCAL)
        require(c1>0 and (u0,u1,v0,v1)==(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1)),
            'Coefficient or exact geometry mismatch.')
        require(0<u0<u1<R and 0<v0<v1<R and u0+v0>R,'Box not in fixed smooth exterior branch.')
        require(j['precision_bits']==m['precision_bits']==192 and m['one_sided_only'] is True,'Precision/one-sided label mismatch.')
        require(old['Psi']['unscaled_by_c'] is True and old['Psi']['Bernstein_degree']==[3,3],
            'Reused Psi model is not the unscaled degree-(3,3) array.')
        require(m['reused_T1_Bernstein']==old['Psi']['bernstein'],'T1 coefficients were changed.')
        require(j['derivatives_saved']==m['derivatives'] and len(m['derivatives'])==2,'Derivative saving mismatch.')
        absolute_uppers=[]
        for row,(a,b) in zip(m['derivatives'],((u0,v0),(v0,u0))):
            guard();require((Q(row['a']),Q(row['b']),row['kind'])==(a,b,'Hab'),'Wrong Hab request.')
            lo,hi,au=map(Q,(row['lower'],row['upper'],row['absolute_upper']))
            require(lo<=hi and au>=max(abs(lo),abs(hi)),'Invalid saved absolute interval upper bound.')
            absolute_uppers.append(au)
        muv=sum(absolute_uppers,Q(0));du,dv=(u1-u0)/2,(v1-v0)/2;deduction=muv*du*dv
        require(Q(m['Muv_upper'])==muv and Q(m['du'])==du and Q(m['dv'])==dv,'Muv/halfwidth arithmetic mismatch.')
        require(Q(m['Psi_one_sided_lower_deduction'])==Q(j['Psi_one_sided_lower_deduction'])==deduction,
            'One-sided correction arithmetic mismatch.')
        old_ep,es=Q(old['Psi']['remainder_upper']),Q(old['AS']['remainder_upper'])
        require(0<=deduction and 0<=old_ep and 0<=es,'Negative upper deduction.')
        require(Q(m['old_Psi_absolute_remainder_upper'])==Q(j['old_Psi_absolute_remainder_upper'])==old_ep and
            Q(m['AS_remainder_upper'])==es,'Old remainder reference mismatch.')
        require(m['deduction_not_greater_than_old_absolute_remainder']==(deduction<=old_ep),
            'Old/new deduction comparison mismatch.')
        require(len(m['rows'])==len(old['rows'])==len(old['Psi']['bernstein'])==len(old['AS']['bernstein'])==16,
            'Expected 16 rows and coefficients.')
        slacks=[];old_slacks=[];diffs=[]
        for ell,row in enumerate(m['rows']):
            guard();index=[ell//4,ell%4];pb,sb,old_row=old['Psi']['bernstein'][ell],old['AS']['bernstein'][ell],old['rows'][ell]
            require(row['bernstein_index']==pb['index']==sb['index']==old_row['bernstein_index']==index,'Coefficient index mismatch.')
            pl,pu,sl,su=map(Q,(pb['lower'],pb['upper'],sb['lower'],sb['upper']))
            require(pl<=pu and sl<=su,'Reversed reused coefficient interval.')
            cc,ss=pl-deduction,sl-es;slack=c1*cc+ss;old_slack=c1*(pl-old_ep)+ss
            require((Q(old_row['c_column']),Q(old_row['s_column']),Q(old_row['baseline_slack']))==(pl-old_ep,ss,old_slack),
                'Old absolute row arithmetic mismatch.')
            require((Q(row['c_column']),Q(row['s_column']),Q(row['baseline_slack']),Q(row['old_absolute_slack']),Q(row['difference']))==
                (cc,ss,slack,old_slack,slack-old_slack),'New one-sided row arithmetic mismatch.')
            slacks.append(slack);old_slacks.append(old_slack);diffs.append(slack-old_slack);answer['rows_checked']+=1
        all_positive=all(v>0 for v in slacks)
        require(j['accepted']==all_positive and j['exact_rows_checked']==16,'Acceptance flag mismatch.')
        minimum=min(slacks);minimum_index=list(divmod(slacks.index(minimum),4))
        require(Q(j['minimum_slack'])==minimum and j['minimum_index']==minimum_index,'Minimum row summary mismatch.')
        require(j['counts']=={kind:{'started':2,'completed':2,'failed':0,'interrupted':0} for kind in ('Hint','acb_integral')},
            'Incomplete or unexpected call counts.')
        trace=[json.loads(line) for line in (ROOT/'calls.jsonl').read_text().splitlines()]
        require(len(trace)==8 and j['calls_started_total']==4,'Expected four calls/eight trace events.')
        require(j['last_event']==trace[-1],'Last event mismatch.')
        ids={event['id'] for event in trace};require(ids=={0,1,2,3},'Unexpected event IDs.')
        for cid in range(4):
            guard();events=[event for event in trace if event['id']==cid]
            require(len(events)==2 and events[0]['status']=='started' and events[1]['status']=='completed',
                'Event pair missing or not completed.')
            before,after=events
            require(all(after[key]==value for key,value in before.items() if key!='status'),
                'Call metadata changed between start and finish.')
            require(after['kind']==('Hint' if cid%2==0 else 'acb_integral') and 0<=after['elapsed_seconds']<=60,
                'Wrong call order/kind/time.')
            require(after['Hint_id']==cid//2,'Integral/Hint parent mismatch.')
            if cid%2==0:
                a,b=((u0,v0),(v0,u0))[cid//2]
                require((Q(after['a']),Q(after['b']),after['derivative_kind'])==(a,b,'Hab'),'Wrong trace derivative input.')
        starts=[(e['id'],e['status']) for e in trace]
        require(starts==[(0,'started'),(1,'started'),(1,'completed'),(0,'completed'),
            (2,'started'),(3,'started'),(3,'completed'),(2,'completed')],'Invalid call nesting.')
        require(all(j[key]==0 for key in ('new_AS_calls','full_Psi_model_calls','new_const_calls','F0_gradient_calls',
            'point_calls','subdivisions','T_models','LP_calls','full_checker_calls')),'Forbidden new work recorded.')
        guard()
        answer.update(status='passed_saved_arithmetic_and_execution_audit',accepted=all_positive,
            minimum_slack=str(minimum),minimum_index=minimum_index,minimum_old_absolute_slack=str(min(old_slacks)),
            minimum_gain=str(min(diffs)),rows_improved=sum(v>0 for v in diffs),rows_worsened=sum(v<0 for v in diffs),
            Hab_intervals_checked=2,absolute_upper_inequalities_checked=2,Muv_upper=str(muv),
            Psi_one_sided_lower_deduction=str(deduction),old_Psi_absolute_remainder_upper=str(old_ep),
            trace_events=8,Hint_calls=2,integral_calls=2,
            result_sha256=digest(ROOT/'result.json'),signed_model_sha256=digest(ROOT/'signed_model.json'),
            run_sha256=digest(ROOT/'signed.run.json'),log_sha256=digest(ROOT/'signed.log'),
            calls_sha256=digest(ROOT/'calls.jsonl'),source_sha256=source_hash,wrapper_sha256=digest(ROOT/'run_once.py'))
    except Exception as exc:answer.update(status='failed',reason=f'{type(exc).__name__}: {exc}')
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        answer.update(elapsed_seconds=monotonic()-start,verifier_sha256=digest(Path(__file__)))
        temporary=output.with_suffix('.json.tmp');temporary.write_text(json.dumps(answer,indent=2,allow_nan=False)+'\n');temporary.replace(output)
    print(json.dumps(answer));return 0 if answer['status'].startswith('passed_') else 1


if __name__=='__main__':sys.exit(main())
