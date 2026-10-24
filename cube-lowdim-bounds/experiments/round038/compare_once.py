#!/usr/bin/env python3
"""One fixed exterior Psi model; reuse node615 AS exactly. No new AS/const/LP.

Instrumentation delegates unchanged arguments/callbacks to the pinned numerical
functions. The module-local acb proxy changes only construction/call accounting.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS,PRECISION=60.,192
FIXED_SHA='4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05'
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
OLD_MODEL_SHA='ca863a29c6cd9fc479c1e198ae8c6af9f470d12b655f7e2a9120e55323dc4800'
OLD_COVER_SHA='1b545178105b1b344557e1b0540f0ce5af126555df62cd3cd021ca95556ab6b4'
LOCAL=('27/256','7/64','35/128','9/32')
PHYSICAL=('13339/8192','3335/2048','163/4096','41/1024')


class BudgetStop(Exception):pass


def dump(path,value):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temp.replace(path)


def load(path):
    spec=importlib.util.spec_from_file_location('round19_fixed_hlcert',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
    spec.loader.exec_module(module);return module


def endpoint(value):
    if not value.is_finite() or not value.is_exact():raise ValueError('Nonfinite/inexact endpoint.')
    m,e=value.man_exp();m,e=int(m),int(e)
    return Q(m*(1<<e)) if e>=0 else Q(m,1<<(-e))


def bounds(value):
    if not value.is_finite():raise ValueError('Nonfinite Arb result.')
    lo,hi=endpoint(value.lower()),endpoint(value.upper())
    if lo>hi:raise ValueError('Reversed interval.')
    return lo,hi


def main():
    if any((ROOT/name).exists() for name in ('result.json','new_model.json','calls.jsonl')):
        raise SystemExit('Existing evidence; no repeat/resume.')
    start=monotonic();path=ROOT/'result.json';trace_path=ROOT/'calls.jsonl'
    counts={kind:{k:0 for k in ('started','completed','failed','interrupted')}
        for kind in ('Psi_model','Hint','acb_integral')}
    result={'status':'started','limits':{'seconds_hard':60,'Arb_bits':192,
        'Psi_models':1,'Hint_calls':12,'acb_integral_calls':12},'counts':counts,
        'calls':[],'in_flight':[],'new_AS_calls':0,'new_const_calls':0,'point_calls':0,
        'subdivisions':0,'T_models':0,'LP_calls':0,'full_checker_calls':0,
        'component_policy':'completed means finite return; nonfinite/exception is failed; no retry.',
        'instrumentation':'Delegates original acb construction, integral callback/args and bound _Hint unchanged.'}
    def save():
        result['elapsed_seconds']=monotonic()-start;dump(path,result)
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_60s_deadline')
    def verify(relative,wanted):
        guard();p=ROOT.parent/relative;actual=sha256(p.read_bytes()).hexdigest()
        if actual!=wanted:raise ValueError('Pinned SHA mismatch: '+relative)
        result.setdefault('inputs',[]).append({'file':relative,'sha256':actual});return p
    def emit(row):
        with trace_path.open('a') as stream:
            stream.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n')
    save()
    try:
        fixed=verify('round17/fixed_checker/hlcert.py',FIXED_SHA)
        baseline_path=verify('round01/scaled_certificate.json',BASELINE_SHA)
        old_path=verify('round18/models.local/node_00615.json',OLD_MODEL_SHA)
        cover_path=verify('round18/leaf_cover.json',OLD_COVER_SHA)
        old=json.loads(old_path.read_text());cover=json.loads(cover_path.read_text())
        node=cover['nodes'][615]
        if (node['id'],node['parent'],node['depth'],node['local'],node['physical'],node['model_sha256'])!=(
                615,[52,1],4,list(LOCAL),list(PHYSICAL),OLD_MODEL_SHA):
            raise ValueError('Source cover node mismatch.')
        if (old['node_id'],old['parent'],old['depth'],old['local'],old['physical'],old['baseline_sha256'])!=(
                615,[52,1],4,list(LOCAL),list(PHYSICAL),BASELINE_SHA):
            raise ValueError('Source model identity mismatch.')
        base=json.loads(baseline_path.read_text());alpha,R,h,c1=map(Q,(base['alpha'],base['R'],base['h'],base['c']))
        if (alpha,R,h,base['N'])!=(Q(47,40),Q(33,20),Q(1,32),80) or c1<=0:
            raise ValueError('Baseline parameters changed.')
        u0,u1,v0,v1=map(Q,PHYSICAL);x0,x1,y0,y1=map(Q,LOCAL)
        if not (0<u0<u1<R and 0<v0<v1<R and u0+v0>R):
            raise ValueError('Fixed same-branch exterior geometry failed.')
        if (u0,u1,v0,v1)!=(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1)):
            raise ValueError('Physical/local geometry mismatch.')
        old_bs=old['AS']['bernstein'];es=Q(old['AS']['remainder_upper'])
        if len(old_bs)!=16 or len(old['rows'])!=16 or es<0:raise ValueError('Invalid saved AS model.')
        for ell,b in enumerate(old_bs):
            if b['index']!=[ell//4,ell%4] or Q(b['lower'])>Q(b['upper']):
                raise ValueError('Saved AS indexing/interval invalid.')
            row=old['rows'][ell]
            if row['bernstein_index']!=b['index'] or Q(row['c_column'])!=Q(old['const']['lower']) or \
                    Q(row['s_column'])!=Q(b['lower'])-es or \
                    Q(row['baseline_slack'])!=c1*Q(row['c_column'])+Q(row['s_column']):
                raise ValueError('Saved const row arithmetic mismatch.')
        M=load(fixed)  # Pinned import sets 96 bits; all objects below use the reset 192 bits.
        M.ctx.prec=PRECISION
        guard();kernel=M.Kern(base);guard()
        result.update(precision_bits=M.ctx.prec,python=sys.version,
            node_id=615,local=list(LOCAL),physical=list(PHYSICAL),baseline_c=str(c1))
        limits={'Psi_model':1,'Hint':12,'acb_integral':12}
        active_hint=[]

        def recorded(kind,fn,finite,metadata):
            guard()
            if counts[kind]['started']>=limits[kind]:raise BudgetStop(kind+'_cap')
            call_id=len(result['calls']);started=monotonic()
            row={'id':call_id,'kind':kind,'status':'started',**metadata}
            result['calls'].append(row);result['in_flight'].append(call_id)
            counts[kind]['started']+=1;save();emit(row.copy())
            try:
                value=fn()
                if not finite(value):
                    row['returned_display']=str(value)
                    raise ValueError('Nonfinite '+kind+' return.')
                counts[kind]['completed']+=1;row['status']='completed'
                if kind!='Psi_model':row['returned_display']=value.str(70)
                row['elapsed_seconds']=monotonic()-started;save();guard();return value
            except BudgetStop as exc:
                if row['status']=='started':counts[kind]['interrupted']+=1;row['status']='interrupted'
                row['reason']=str(exc);raise
            except Exception as exc:
                if row['status']=='started':counts[kind]['failed']+=1;row['status']='failed'
                row['reason']=f'{type(exc).__name__}: {exc}';raise
            finally:
                row['elapsed_seconds']=monotonic()-started
                result['in_flight'].remove(call_id);emit(row.copy());save()

        real_acb=M.acb
        class AcbProxy:
            def __call__(self,*args,**kwargs):return real_acb(*args,**kwargs)
            def integral(self,*args,**kwargs):
                return recorded('acb_integral',lambda:real_acb.integral(*args,**kwargs),
                    lambda value:value.is_finite(),{'Hint_id':active_hint[-1] if active_hint else None})
        original_hint=kernel._Hint
        def hint(a,b,kind):
            hint_id=counts['Hint']['started']
            active_hint.append(hint_id)
            try:
                return recorded('Hint',lambda:original_hint(a,b,kind),lambda value:value.is_finite(),
                    {'Hint_id':hint_id,'a':str(a),'b':str(b),'derivative_kind':kind})
            finally:active_hint.pop()
        def finite_model(value):
            B,E=value
            return B is not None and E is not None and E.is_finite() and len(B)==4 and \
                all(len(row)==4 and all(v.is_finite() for v in row) for row in B)
        M.acb=AcbProxy();kernel._Hint=hint
        try:
            B,EPsi=recorded('Psi_model',lambda:kernel._outside_taylor(*map(M.Q,PHYSICAL)),
                finite_model,{'node_id':615,'unscaled_by_c':True,'Bernstein_degree':[3,3]})
        finally:M.acb=real_acb;kernel._Hint=original_hint
        guard();epl,epu=bounds(EPsi)
        if epl<0:raise ValueError('Psi remainder is not nonnegative.')
        model={'node_id':615,'physical':list(PHYSICAL),'local':list(LOCAL),
            'source_hlcert_sha256':FIXED_SHA,'source_AS_model_sha256':OLD_MODEL_SHA,
            'source_cover_sha256':OLD_COVER_SHA,'baseline_sha256':BASELINE_SHA,
            'precision_bits':PRECISION,'unscaled_by_c':True,'Bernstein_degree':[3,3],
            'Psi_remainder_lower':str(epl),'Psi_remainder_upper':str(epu),
            'Psi_remainder_display':EPsi.str(70),'AS_remainder_upper':str(es),
            'original_const':old['const'],'Psi_bernstein':[],'rows':[]}
        for ell in range(16):
            guard();i,j=ell//4,ell%4;lo,hi=bounds(B[i][j])
            model['Psi_bernstein'].append({'index':[i,j],'lower':str(lo),'upper':str(hi),'display':B[i][j].str(70)})
            ccol=lo-epu;scol=Q(old_bs[ell]['lower'])-es;slack=c1*ccol+scol
            old_slack=Q(old['rows'][ell]['baseline_slack'])
            model['rows'].append({'bernstein_index':[i,j],'c_column':str(ccol),'s_column':str(scol),
                'baseline_slack':str(slack),'old_const_slack':str(old_slack),'difference':str(slack-old_slack)})
        guard();model_path=ROOT/'new_model.json';dump(model_path,model)
        worst=min(model['rows'],key=lambda row:Q(row['baseline_slack']))
        result.update(status='completed',accepted=all(Q(row['baseline_slack'])>0 for row in model['rows']),
            minimum_slack=worst['baseline_slack'],minimum_index=worst['bernstein_index'],
            exact_rows_checked=16,new_model_sha256=sha256(model_path.read_bytes()).hexdigest())
    except BudgetStop as exc:result.update(status='interrupted_saved',reason=str(exc),accepted=False)
    except Exception as exc:result.update(status='failed',reason=f'{type(exc).__name__}: {exc}',accepted=False)
    for kind in counts:
        assert counts[kind]['started']==sum(counts[kind][key] for key in ('completed','failed','interrupted'))
        assert counts[kind]['started']<={'Psi_model':1,'Hint':12,'acb_integral':12}[kind]
    result['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    if trace_path.exists():result['calls_sha256']=sha256(trace_path.read_bytes()).hexdigest()
    save();print(json.dumps({'status':result['status'],'accepted':result.get('accepted',False),
        'counts':counts,'elapsed_seconds':result['elapsed_seconds']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
