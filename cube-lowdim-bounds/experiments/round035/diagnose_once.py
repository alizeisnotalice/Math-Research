#!/usr/bin/env python3
"""Two fixed points x two reviewed implementations, <=4 scalars/60s.

The AK point is outside r=R: first uses an analytic-aware local adaptation of
the original finite acb integral, second the public hypergeometric point formula.
The public integral callback is not invoked. No full checker/LP.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
PRECISION,MAX_CALLS,SECONDS=192,4,60.
KERNEL_SHA='f46a31bb2b2356651a2e255ad63f058ddd4ce857e14f1365d07d8cabd24de9dd'
STRUCTURED_SHA='01b310b798ea6fc3f3bb93f3da226a3a6af68b34a54e6dfa728a6995e7a6c67d'
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
REPORT_SHA='f2b655e4441cbbee021675df9e364e0517e4efc10288a4e680fd9e8b5da3d3d8'
ADAPTER_SHA='028f2d38b0fbbd61c890b56d62fa870976b7ae7b56b567c97dbcb691de7b305b'
ORIGINAL_SHA='4bd7e823bdcb336f8e20c188d68d6d8c1e4aa22ec52353885059c24445317db5'
INDEPENDENT_SHA='c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'
METHODS=('original_hlcert_adapter','independent_hl_core_point_methods')
K_POINT=(Q(47,64),Q(47,64))
A_POINT=(Q(106841,65536),Q(2585,65536))


class BudgetStop(Exception):pass


def dump(path,value):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n');temp.replace(path)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def original_APsi_outside_analytic(M,K,u,v):
    """Original exterior formula with the acb.integral analytic contract fixed locally.

    The real integration paths satisfy s+b>0 and a-s>=r-R>0.
    Complex balls still require analytic=analytic on every logarithm.
    Upstream code is pinned for provenance, but its callback is not called.
    """
    alpha,R=K.al,K.R;total=M.arb(0)
    U,V=M.Q(u),M.Q(v)
    for a,b in ((U,V),(V,U)):
        L=K.R_q-b
        if L<=0:continue
        if not b>0 or not a-L>0:raise ValueError('Exterior integration path meets a singularity.')
        aa,bb,alc=M.acb(a),M.acb(b),M.acb(alpha)
        Ra=M.acb(M.apow(R,-alpha))  # Positive real constant, outside the callback.
        def integrand(s,analytic):
            return ((-(s+bb).log(analytic=analytic)*alc).exp()-Ra)*(
                (-(aa-s).log(analytic=analytic)*(1+alc)).exp()+
                (-(aa+s).log(analytic=analytic)*(1+alc)).exp())
        value=M.acb.integral(integrand,M.acb(0),M.acb(L))
        total+=2*value.real
    return total


def original_AK_outside(adapter,u,v,exact_B):
    """Reviewed exact-B/D spline sum plus analytic-aware local exterior integral."""
    M,K=adapter.module,adapter.kernel
    if not u>0 or not v>0 or not u+v>Q(str(K.R_q)):
        raise ValueError('Exterior adapter requires positive coordinates and r>R.')
    A1,A2=M.arb(0),M.arb(0)
    for comp in K.comps:
        h=Q(str(comp.h_q));s,t=u/h,v/h
        bu={i:exact_B(s-i) for i in range(-comp.N,comp.N+1)}
        bv={j:exact_B(t-j) for j in range(-comp.N,comp.N+1)}
        bu={i:M.Q(b) for i,b in bu.items() if b};bv={j:M.Q(b) for j,b in bv.items() if b}
        one,two=M.arb(0),M.arb(0)
        for m,line in comp.D.items():
            Dv=sum((line[j+comp.N]*b for j,b in bv.items()),M.fmpq(0))
            Du=sum((line[i+comp.N]*b for i,b in bu.items()),M.fmpq(0))
            # exponent 3-alpha>0: exact zero terms contribute 0, not 0*log(0).
            if Dv and s!=m:one+=M.apow(M.arb(M.Q(abs(s-m))),comp.p)*M.arb(Dv)
            if Du and t!=m:two+=M.apow(M.arb(M.Q(abs(t-m))),comp.p)*M.arb(Du)
        A1+=comp.sc*one;A2+=comp.sc*two
    singular=K.c*original_APsi_outside_analytic(M,K,u,v)
    return singular+A1+A2,{'c_A_Psi':singular,'A1S':A1,'A2S':A2,'AS':A1+A2}


def main():
    if any((ROOT/name).exists() for name in ('result.json','points.json','calls.jsonl')):
        raise SystemExit('Existing diagnostic evidence; no repeated point evaluation.')
    start=monotonic();current=None
    counts={'started':0,'completed':0,'failed':0,'interrupted':0,'not_started':MAX_CALLS}
    result={'status':'started','counts':counts,'rows':[],'precision_bits':PRECISION,'maximum_calls':MAX_CALLS,
        'call_unit':'One complete K or AK scalar at the fixed point with one implementation',
        'LP_calls':0,'full_checker_calls':0,'new_grid_points':0,
        'independence_scope':'Two code paths share python-flint/Arb and reviewed spline formulas; exterior APsi uses direct finite integration versus hypergeometric representation.'}
    result_path=ROOT/'result.json';trace_path=ROOT/'calls.jsonl';dump(result_path,result);trace=trace_path.open('x')
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_60s_deadline')
    def save():
        counts['not_started']=MAX_CALLS-counts['started'];result['elapsed_seconds']=monotonic()-start;dump(result_path,result)
    def verify(path,wanted):
        guard()
        if sha256(path.read_bytes()).hexdigest()!=wanted:raise ValueError('Pinned source SHA mismatch: '+path.name)
    def emit(row):trace.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n');trace.flush()
    try:
        kernel_path=ROOT.parent/'round15/candidate.json';report_path=ROOT.parent/'round15/first_report.json'
        verify(kernel_path,KERNEL_SHA);verify(report_path,REPORT_SHA)
        kernel=json.loads(kernel_path.read_text());report=json.loads(report_path.read_text())
        if (Q(kernel['alpha']),Q(kernel['R']),Q(kernel['h']),kernel['N'])!=(Q(47,40),Q(33,20),Q(1,32),80):
            raise ValueError('Fixed kernel analytic parameters changed.')
        structured_path=ROOT.parent/'round15/structured_kernel.json';verify(structured_path,STRUCTURED_SHA)
        baseline_path=ROOT.parent/'round01/scaled_certificate.json';verify(baseline_path,BASELINE_SHA)
        structured=json.loads(structured_path.read_text());baseline=json.loads(baseline_path.read_text())
        if any(structured[key]!=kernel[key] for key in ('alpha','R','h','N','c','a')):
            raise ValueError('Plain/structured kernel field binding changed.')
        structure=structured['structure'];scale=Q(structure['s']);coef=structured['a'];T=structure['T']
        if structure['baseline_sha256']!=BASELINE_SHA or scale<=0 or Q(kernel['c'])<=0:
            raise ValueError('Positive c/s component extraction requires the frozen structure.')
        if any(sum(map(int,key.split(',')))>20 for key,value in T.items() if Q(value)):
            raise ValueError('Frozen T support is not inside r<=3/4.')
        if any(Q(coef.get(key,'0'))!=scale*Q(baseline['a'].get(key,'0'))+Q(T.get(key,'0')) for key in set(coef)|set(T)|set(baseline['a'])):
            raise ValueError('Exact a=s*a1+T binding failed.')
        kboxes=[tuple(Q(x) for x in row[:4]) for row in report['K_fail']]
        if len(kboxes)!=4 or not all(box[0]<=K_POINT[0]<=box[1] and box[2]<=K_POINT[1]<=box[3] for box in kboxes):
            raise ValueError('Fixed K point is not the recorded four-box junction.')
        afails=report['A_fail']
        worst_index,worst=min(enumerate(afails),key=lambda pair:(Q(str(pair[1][4])),tuple(Q(x) for x in pair[1][:4]),pair[0]))
        abox=tuple(Q(x) for x in worst[:4]);center=((abox[0]+abox[1])/2,(abox[2]+abox[3])/2)
        if center!=A_POINT or abox!=(Q(13355,8192),Q(53421,32768),Q(323,8192),Q(1293,32768)):
            raise ValueError('Fixed exterior A point/source box changed; no alternative point.')
        R=Q(kernel['R'])
        if not (1<sum(K_POINT)<R and sum(A_POINT)>R and all(x>0 for x in (*K_POINT,*A_POINT))):
            raise ValueError('Fixed threshold/inside-outside branch assertion failed.')
        points=[{'id':'K_outer','quantity':'K','u':str(K_POINT[0]),'v':str(K_POINT[1]),'r':str(sum(K_POINT)),
            'threshold':'0','branch':'inside_truncation_outside_unit_diamond',
            'source_boxes':[[str(x) for x in box] for box in kboxes],'selection':'Fixed common junction, not a new minimum search.'},
            {'id':'A_outer','quantity':'AK','u':str(A_POINT[0]),'v':str(A_POINT[1]),'r':str(sum(A_POINT)),
            'threshold':'0','branch':'outside_truncation','source_box':[str(x) for x in abox],
            'source_index':worst_index,'source_lower_float':worst[4],'selection':'Fixed rational midpoint of the lowest reported failing A box.'}]
        description={'status':'fixed_before_any_point_call','kernel_sha256':KERNEL_SHA,'source_report_sha256':REPORT_SHA,
            'points':points,'maximum_calls':MAX_CALLS,'precision_bits':PRECISION}
        dump(ROOT/'points.json',description)
        result.update(kernel_sha256=KERNEL_SHA,source_report_sha256=REPORT_SHA,
            points_sha256=sha256((ROOT/'points.json').read_bytes()).hexdigest(),points=points)
        helper_path=ROOT/'dependency_round04_diagnose.py'
        if not helper_path.exists():helper_path=ROOT.parent/'round04/diagnose.py'
        verify(helper_path,ADAPTER_SHA);helper=load('round16_reviewed_adapter',helper_path)
        public=ROOT/'public_src'
        if not (public/'hlcert.py').is_file():public=ROOT.parent/'upstream-v4/checker'
        original_path,independent_path=public/'hlcert.py',public/'independent/hl_core.py'
        verify(original_path,ORIGINAL_SHA);verify(independent_path,INDEPENDENT_SHA)
        original=load('round16_original_hlcert',original_path);independent=load('round16_independent_core',independent_path)
        from flint import arb,fmpq,ctx
        import flint
        ctx.prec=PRECISION
        first=helper.OriginalAdapter(original,kernel_path);second=helper.IndependentAdapter(independent,kernel_path)
        result.update(flint_version=flint.__version__,python=sys.version,
            structured_kernel_sha256=STRUCTURED_SHA,baseline_sha256=BASELINE_SHA,s_exact=str(scale),
            implementation_sha256={'round04_adapter':ADAPTER_SHA,'hlcert.py':ORIGINAL_SHA,'hl_core.py':INDEPENDENT_SHA},
            exterior_method='First: local original finite formula, two coordinate acb.integral calls, all three logarithms forward analytic=analytic; public Apsi_outside_point is not invoked. Second: APsi_point exterior hypergeometric formula. Neither calls the inside Apsi branch.')
        save();guard();values={}
        for point in points:
            u,v=Q(point['u']),Q(point['v'])
            for method in METHODS:
                guard()
                if counts['started']>=MAX_CALLS:raise BudgetStop('four_scalar_call_cap')
                current={'call':counts['started']+1,'point_id':point['id'],'quantity':point['quantity'],
                    'u':str(u),'v':str(v),'r':point['r'],'threshold':'0','branch':point['branch'],
                    'implementation':method,'kernel_sha256':KERNEL_SHA,'source_report_sha256':REPORT_SHA,
                    'status':'started','started_seconds':monotonic()-start}
                counts['started']+=1;result['rows'].append(current);emit(current);save()
                try:
                    if method==METHODS[0]:
                        value,parts=first.K(u,v) if point['quantity']=='K' else original_AK_outside(first,u,v,helper.exact_B)
                    else:value,parts=getattr(second,point['quantity'])(u,v)
                    if point['quantity']=='K':
                        # Same evaluation only: T=0 here, so S=s*S1. No baseline point evaluation.
                        cq=Q(kernel['c']);cA=arb(fmpq(cq.numerator,cq.denominator))
                        sA=arb(fmpq(scale.numerator,scale.denominator))
                        parts['Psi']=parts['c_Psi']/cA
                        parts['S1_from_same_evaluation']=parts['S']/sA
                    guard()
                    if not value.is_finite() or not all(part.is_finite() for part in parts.values()):
                        raise ValueError('Nonfinite point interval; no sign conclusion.')
                    ball=value.str(70);recorded=arb(ball)
                    if recorded.upper()<0:classification='strict_negative_at_fixed_point'
                    elif recorded.lower()>=0:classification='nonnegative_at_this_point_only'
                    else:classification='inconclusive_interval_crosses_zero'
                    current.update(status='completed',ball=ball,lower_endpoint_ball=recorded.lower().str(70),
                        upper_endpoint_ball=recorded.upper().str(70),classification=classification,
                        parts={name:part.str(70) for name,part in parts.items()},completed_seconds=monotonic()-start)
                    if point['quantity']=='K':
                        current['component_derivation']='At r=47/32>3/4, exact structure has T=0. Psi=(cPsi)/c and S1=S/s use already computed terms; no baseline K or ratio -S1/Psi evaluated.'
                except BudgetStop as exc:
                    counts['interrupted']+=1;current.update(status='interrupted',reason=str(exc));emit(current);save();current=None;raise
                except Exception as exc:
                    counts['failed']+=1;current.update(status='failed_unknown',reason=f'{type(exc).__name__}: {exc}');emit(current);save();current=None;raise
                counts['completed']+=1;values[point['id'],method]=recorded;emit(current);save();current=None
        comparisons=[]
        for point in points:
            a,b=(values[point['id'],method] for method in METHODS);difference=a-b
            comparisons.append({'point_id':point['id'],'intervals_overlap':not(a.upper()<b.lower() or b.upper()<a.lower()),
                'difference_contains_zero':not(difference.lower()>0 or difference.upper()<0),'difference_ball':difference.str(70)})
        result.update(status='four_fixed_point_calls_closed',cross_implementation_comparisons=comparisons,
            continuity_argument='K is continuous away from the origin; the fixed exterior AK point is away from axes/R and has a locally continuous APsi branch plus finite continuous spline generator. Any strict upper<0 implies an open negative neighborhood; for AK the reviewed away-origin distribution bridge then excludes the super-solution condition. No neighborhood size or new C2 bound computed.')
    except BudgetStop as exc:result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:result.update(status='failed_preparation_or_diagnostic',reason=f'{type(exc).__name__}: {exc}')
    finally:trace.close()
    assert counts['started']<=MAX_CALLS
    assert counts['started']==counts['completed']+counts['failed']+counts['interrupted']
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),trace_sha256=sha256(trace_path.read_bytes()).hexdigest())
    save();print(json.dumps({'status':result['status'],'counts':counts,'seconds':result['elapsed_seconds']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
