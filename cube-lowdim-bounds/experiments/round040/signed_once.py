#!/usr/bin/env python3
"""Two Hab integrals; reuse saved T1/AS Bernstein balls. One-sided lower model."""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS,PRECISION=60.,192
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
FIXED_SHA='4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05'
MODEL_SHA='8dd22932cefd30806acf819ffae37f10d9497e50e2e798361e10f0d9e207db26'
COVER_SHA='6028e48efe444e27c5259858cf1eeb8b41e0081b6e39f06102439ddb198c4cac'
LOCAL=('21/128','11/64','1/4','33/128')
PHYSICAL=('6677/4096','3339/2048','5/128','161/4096')


class BudgetStop(Exception):pass


def dump(path,value):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n');tmp.replace(path)


def endpoint(value):
    if not value.is_finite() or not value.is_exact():raise ValueError('Expected finite exact endpoint.')
    m,e=value.man_exp();m,e=int(m),int(e)
    return Q(m*(1<<e)) if e>=0 else Q(m,1<<(-e))


def bounds(value):
    if not value.is_finite():raise ValueError('Nonfinite Arb result.')
    lo,hi=endpoint(value.lower()),endpoint(value.upper())
    if lo>hi:raise ValueError('Reversed interval.')
    return lo,hi


def main():
    if any((ROOT/name).exists() for name in ('result.json','signed_model.json','calls.jsonl')):
        raise SystemExit('Existing evidence: no repeat/resume.')
    start=monotonic();path=ROOT/'result.json';trace_path=ROOT/'calls.jsonl'
    counts={kind:{key:0 for key in ('started','completed','failed','interrupted')} for kind in ('Hint','acb_integral')}
    result={'status':'started','counts':counts,'in_flight':[],'calls_started_total':0,
        'limits':{'seconds_hard':60,'Arb_bits':192,'Hint':2,'acb_integral':2,'exact_rows':16},
        'new_AS_calls':0,'full_Psi_model_calls':0,'new_const_calls':0,'F0_gradient_calls':0,
        'point_calls':0,'subdivisions':0,'T_models':0,'LP_calls':0,'full_checker_calls':0,
        'interpretation':'One-sided lower correction from nonnegative diagonal curvature, not an absolute error bound.'}
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_60s_deadline')
    def save():result['elapsed_seconds']=monotonic()-start;dump(path,result)
    def verify(relative,wanted):
        guard();p=ROOT.parent/relative;actual=sha256(p.read_bytes()).hexdigest()
        if actual!=wanted:raise ValueError('Pinned SHA mismatch: '+relative)
        result.setdefault('inputs',[]).append({'file':relative,'sha256':actual});return p
    save();M=None;kernel=None;real_acb=None;original_hint=None;trace=None
    try:
        source=verify('round17/fixed_checker/hlcert.py',FIXED_SHA)
        model_path=verify('round20/models.local/node_00661.json',MODEL_SHA)
        cover_path=verify('round20/leaf_cover.json',COVER_SHA)
        base_path=verify('round01/scaled_certificate.json',BASELINE_SHA)
        old=json.loads(model_path.read_text());cover=json.loads(cover_path.read_text());base=json.loads(base_path.read_text())
        node=cover['nodes'][661]
        if (node['id'],node['parent'],node['depth'],node['local'],node['physical'],node['model_sha256'])!=(
                661,[52,1],4,list(LOCAL),list(PHYSICAL),MODEL_SHA):raise ValueError('Frozen cover node mismatch.')
        if (old['node_id'],old['parent'],old['depth'],old['local'],old['physical'],old['baseline_sha256'],old['fixed_hlcert_sha256'])!=(
                661,[52,1],4,list(LOCAL),list(PHYSICAL),BASELINE_SHA,FIXED_SHA):raise ValueError('Frozen model identity mismatch.')
        c1,h,R=map(Q,(base['c'],base['h'],base['R']))
        if (Q(base['alpha']),h,R,base['N'])!=(Q(47,40),Q(1,32),Q(33,20),80) or c1<=0:
            raise ValueError('Baseline parameters changed.')
        u0,u1,v0,v1=map(Q,PHYSICAL);x0,x1,y0,y1=map(Q,LOCAL)
        if (u0,u1,v0,v1)!=(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1)) or \
                not (0<u0<u1<R and 0<v0<v1<R and u0+v0>R):raise ValueError('Fixed same-branch geometry failed.')
        if old['Psi']['unscaled_by_c'] is not True or old['Psi']['Bernstein_degree']!=[3,3] or \
                any(len(old[key]['bernstein'])!=16 for key in ('Psi','AS')) or len(old['rows'])!=16:
            raise ValueError('Saved first-order Bernstein schema changed.')
        es,old_ep=Q(old['AS']['remainder_upper']),Q(old['Psi']['remainder_upper'])
        if es<0 or old_ep<0:raise ValueError('Invalid saved remainder.')
        du,dv=(u1-u0)/2,(v1-v0)/2
        spec=importlib.util.spec_from_file_location('round21_fixed_hlcert',source)
        M=importlib.util.module_from_spec(spec);sys.modules[spec.name]=M;spec.loader.exec_module(M)
        M.ctx.prec=PRECISION
        guard();kernel=M.Kern(base);guard()
        result.update(precision_bits=M.ctx.prec,python=sys.version,node_id=661,physical=list(PHYSICAL),local=list(LOCAL))
        trace=trace_path.open('x');active_hint=[]
        def emit(event):trace.write(json.dumps(event,ensure_ascii=False,allow_nan=False)+'\n');trace.flush()
        def recorded(kind,fn,metadata):
            guard()
            if counts[kind]['started']>=2:raise BudgetStop(kind+'_call_cap')
            cid=result['calls_started_total'];result['calls_started_total']+=1;ts=monotonic()
            event={'id':cid,'kind':kind,'status':'started',**metadata}
            counts[kind]['started']+=1;result['in_flight'].append(event.copy());emit(event);save()
            try:
                value=fn()
                if not value.is_finite():raise ValueError('Nonfinite '+kind+' return.')
                event.update(status='completed',returned_display=value.str(70));counts[kind]['completed']+=1
                return value
            except BudgetStop as exc:
                event.update(status='interrupted',reason=str(exc));counts[kind]['interrupted']+=1;raise
            except Exception as exc:
                event.update(status='failed',reason=f'{type(exc).__name__}: {exc}');counts[kind]['failed']+=1;raise
            finally:
                event['elapsed_seconds']=monotonic()-ts
                result['in_flight']=[e for e in result['in_flight'] if e['id']!=cid]
                result['last_event']=event;emit(event);save()
                if event['status']=='completed':guard()
        real_acb=M.acb;original_hint=kernel._Hint
        class AcbProxy:
            def __call__(self,*args,**kwargs):return real_acb(*args,**kwargs)
            def integral(self,*args,**kwargs):
                return recorded('acb_integral',lambda:real_acb.integral(*args,**kwargs),
                    {'Hint_id':active_hint[-1] if active_hint else None})
        def hint(a,b,kind):
            hid=counts['Hint']['started'];active_hint.append(hid)
            try:return recorded('Hint',lambda:original_hint(a,b,kind),
                {'Hint_id':hid,'a':str(a),'b':str(b),'derivative_kind':kind})
            finally:active_hint.pop()
        M.acb=AcbProxy();kernel._Hint=hint
        derivative_rows=[]
        for a,b in ((u0,v0),(v0,u0)):
            guard();value=kernel._Hint(M.Q(a),M.Q(b),'Hab')
            lo,hi=bounds(value);_,abs_upper=bounds(abs(value))
            derivative_rows.append({'a':str(a),'b':str(b),'kind':'Hab','lower':str(lo),'upper':str(hi),
                'absolute_upper':str(abs_upper),'display':value.str(70)})
            result['derivatives_saved']=derivative_rows;save()
        guard();muv=sum((Q(row['absolute_upper']) for row in derivative_rows),Q(0));deduction=muv*du*dv
        model={'node_id':661,'physical':list(PHYSICAL),'local':list(LOCAL),'precision_bits':PRECISION,
            'source_model_sha256':MODEL_SHA,'source_cover_sha256':COVER_SHA,'fixed_hlcert_sha256':FIXED_SHA,
            'baseline_sha256':BASELINE_SHA,'derivatives':derivative_rows,'Muv_upper':str(muv),
            'du':str(du),'dv':str(dv),'Psi_one_sided_lower_deduction':str(deduction),
            'old_Psi_absolute_remainder_upper':str(old_ep),'AS_remainder_upper':str(es),
            'deduction_not_greater_than_old_absolute_remainder':deduction<=old_ep,
            'one_sided_only':True,'reused_T1_Bernstein':old['Psi']['bernstein'],'rows':[]}
        for ell in range(16):
            guard();pb,sb,old_row=old['Psi']['bernstein'][ell],old['AS']['bernstein'][ell],old['rows'][ell]
            index=[ell//4,ell%4]
            if pb['index']!=index or sb['index']!=index or old_row['bernstein_index']!=index:
                raise ValueError('Saved Bernstein indexing mismatch.')
            pl,pu,sl,su=map(Q,(pb['lower'],pb['upper'],sb['lower'],sb['upper']))
            if pl>pu or sl>su:raise ValueError('Reversed saved coefficient interval.')
            scol=sl-es;old_slack=c1*(pl-old_ep)+scol
            if Q(old_row['c_column'])!=pl-old_ep or Q(old_row['s_column'])!=scol or \
                    Q(old_row['baseline_slack'])!=old_slack:raise ValueError('Old saved row arithmetic mismatch.')
            ccol=pl-deduction;slack=c1*ccol+scol
            model['rows'].append({'bernstein_index':index,'c_column':str(ccol),'s_column':str(scol),
                'baseline_slack':str(slack),'old_absolute_slack':str(old_slack),'difference':str(slack-old_slack)})
        guard();model_path=ROOT/'signed_model.json';dump(model_path,model)
        worst=min(model['rows'],key=lambda row:Q(row['baseline_slack']))
        result.update(status='completed',accepted=all(Q(row['baseline_slack'])>0 for row in model['rows']),
            exact_rows_checked=16,minimum_slack=worst['baseline_slack'],minimum_index=worst['bernstein_index'],
            Psi_one_sided_lower_deduction=str(deduction),old_Psi_absolute_remainder_upper=str(old_ep),
            signed_model_sha256=sha256(model_path.read_bytes()).hexdigest())
    except BudgetStop as exc:result.update(status='interrupted_saved',reason=str(exc),accepted=False)
    except Exception as exc:result.update(status='failed',reason=f'{type(exc).__name__}: {exc}',accepted=False)
    finally:
        if M is not None and real_acb is not None:M.acb=real_acb
        if kernel is not None and original_hint is not None:kernel._Hint=original_hint
        if trace is not None:trace.close()
    for kind in counts:
        assert counts[kind]['started']<=2
        assert counts[kind]['started']==sum(counts[kind][key] for key in ('completed','failed','interrupted'))
    result['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    if trace_path.exists():result['calls_sha256']=sha256(trace_path.read_bytes()).hexdigest()
    save();print(json.dumps({'status':result['status'],'accepted':result.get('accepted',False),
        'counts':counts,'elapsed_seconds':result['elapsed_seconds']}),flush=True)


if __name__=='__main__':main()
