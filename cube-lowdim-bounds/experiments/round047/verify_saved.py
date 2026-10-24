"""Independent <=30s saved sparse-dual audit; standard library only.

No generating-code import, optimization, kernel/Arb value, integral or root search.
The true coefficient intervals and exact box imply the checked lower bound.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isfinite
import json
from pathlib import Path
from time import monotonic

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'saved_verification.json'
SEARCH_SHA='715678d73fa80e7eb108d9642a6b02fde242773be7feec8d44931359510731d2'
SOLVER_SHA='d37a5f3194e8fce144f7eef94a2c2f84de03579180a68769abecf1c161ddcd8c'
WRAPPER_SHA='2622068da96b9d24df5f42c38e03093504bf1bcf219313d66533150c176808a8'
BASE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
DIR_SHA='98c0e2fcc7f5d448aef33acc24f2b7fe22977f7a2361b20dc7c2ade926a0fd6e'
POINT_SHA='513c8d953695d9ddd1c7a097ee9af2ab624f79b875330d02e0d42ad2dcc70ebd'
COLUMNS=(2,3,4,5,198,199,201,202,222,244,267,268)
NEW_POINTS=((F(23497,65536),F(10819,65536)),(F(47105,131072),F(21905,131072)))


def integer_ties_even(value):
    """Exact floor/remainder rule, independent of the generator's round(Fraction)."""
    q,r=divmod(value.numerator,value.denominator)
    return q+(2*r>value.denominator or (2*r==value.denominator and q%2!=0))


def main():
    if OUT.exists():raise SystemExit('Existing saved audit: no repeat or overwrite.')
    start=monotonic();state=dict(status='started',complete=False,budget_seconds=30,
        support_rows_checked=0,point_rows_checked=0,structural_rows_checked=0,coefficients_checked=0,
        new_point_values=0,new_Arb_calls=0,new_LP=0,new_integrals=0,integer_root_searches=0)
    def save():
        state['elapsed_seconds']=monotonic()-start
        temp=OUT.with_suffix('.json.tmp');temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(OUT)
    def guard():
        if monotonic()-start>=30:raise TimeoutError('30s saved-data budget; no extension.')
    def read(path):guard();j=json.loads(path.read_text());guard();return j
    def digest(path):
        h=sha256()
        with path.open('rb') as stream:
            while True:
                guard();block=stream.read(1<<20)
                if not block:break
                h.update(block)
        return h.hexdigest()
    def vector(values,n=26):
        assert len(values)==n
        return list(map(F,values))
    def mass(coeffs,h):
        value=F(0)
        for key,a in coeffs.items():
            guard();i,j=map(int,key.split(','));assert 0<=j<=i<=80
            multiplicity=1 if i==j==0 else 4 if j==0 or i==j else 8
            value+=h*h*multiplicity*F(a)/2
        return value
    def desc_identity(desc):
        return desc['kind'],F(desc['u']),F(desc['v']),F(desc['target'])
    save()
    try:
        result,run,cert,sol=(read(ROOT/name) for name in ('result.json','search.run.json','certificate.json','solver_result.json'))
        assert run['status']=='terminal' and run['exit_code']==0 and 0<=run['elapsed_seconds']<300
        assert digest(ROOT/'search.py')==result['source_sha256_before']==result['source_sha256_after']==run['source_sha256_before']==run['source_sha256_after']==SEARCH_SHA
        assert digest(ROOT/'solver.m')==result['solver_source_sha256']==SOLVER_SHA and digest(ROOT/'run_once.py')==WRAPPER_SHA
        assert run['log_sha256']==digest(ROOT/'search.log') and result['trace_sha256']==digest(ROOT/'calls.jsonl')
        assert result['status']=='one_lp_and_sparse_exact_lower_bound_closed' and result['candidate_exported'] is False
        assert result['limits']==dict(seconds=300,column_values=392,LP=1,LP_soft_seconds=20,bits=192)
        assert result['inputs_unchanged'] and result['input_sha256_before']==result['input_sha256_after']
        for path,wanted in result['input_sha256_before'].items():assert digest(ROOT/path)==wanted
        assert digest(ROOT.parent/'round01/scaled_certificate.json')==BASE_SHA
        assert digest(ROOT.parent/'round14/directions.json')==DIR_SHA
        assert digest(ROOT.parent/'round06/point_manifest.json')==POINT_SHA
        assert cert['baseline_sha256']==BASE_SHA and cert['directions_sha256']==DIR_SHA
        assert result['certificate_sha256']==digest(ROOT/'certificate.json')
        assert result['point_enclosures_sha256']==digest(ROOT/'point_enclosures.json')
        model_sha=digest(ROOT/'model.local.mat')
        assert cert['model_sha256']==result['model_sha256']==sol['input_sha256']==model_sha
        assert cert['solver_result_sha256']==result['solver_result_sha256']==digest(ROOT/'solver_result.json')
        assert result['matlab_log_sha256']==digest(ROOT/'matlab.log')
        assert sol['status']=='returned' and sol['exitflag']==result['exitflag']==1 and result['matlab_exit_code']==0
        assert sol['calls_started']==sol['calls_completed']==1
        base=read(ROOT.parent/'round01/scaled_certificate.json');directions=read(ROOT.parent/'round14/directions.json')['directions']
        points=read(ROOT.parent/'round06/point_manifest.json')
        assert len(points['K'])==26604 and len(points['A'])==24970
        alpha,R,h,c1=map(F,(base['alpha'],base['R'],base['h'],base['c']))
        assert (alpha,R,h,base['N'])==(F(47,40),F(33,20),F(1,32),80) and c1>0
        assert len(directions)==12 and tuple(d['old_column'] for d in directions)==COLUMNS
        bounds=[];masses=[]
        for j,d in enumerate(directions):
            expected_bound=F(1,2) if j<4 else F(1,4)
            assert F(d['bound'])==expected_bound
            bounds.append(expected_bound);mm=mass(d['coefficients'],h)
            assert mm==F(d['half_mass'])
            assert all(not F(a) or h*(sum(map(int,k.split(',')))+4)<=F(3,4) for k,a in d['coefficients'].items())
            masses.append(mm)
        lower=[c1/2,F(1,2)]+[-b for b in bounds]+[F(0)]*12
        upper=[3*c1/2,F(3,2)]+bounds+bounds
        assert vector(cert['lower_bounds'])==lower and vector(cert['upper_bounds'])==upper
        rc=cert['root_certificate'];root,scale=int(rc['root_floor']),int(rc['root_scale'])
        assert (rc['root_degree'],rc['root_power'],F(rc['R']),scale)==(40,33,R,10**60)
        assert root**40*R.denominator**33<=R.numerator**33*scale**40<(root+1)**40*R.denominator**33
        q_lower=alpha*F(root,scale)/(2-alpha);q_upper=alpha*F(root+1,scale)/(2-alpha)
        f_lower=[q_lower,mass(base['a'],h)]+masses+[F(0)]*12
        f_upper=[q_upper]+f_lower[1:]
        assert vector(cert['cost_lower_coefficients'])==f_lower and vector(cert['cost_upper_coefficients'])==f_upper
        family=result['exact_family']
        assert (F(family['alpha']),F(family['R']),F(family['h']),family['N'],F(family['c1']))==(alpha,R,h,80,c1)
        assert vector(family['lower_bounds'])==lower and vector(family['upper_bounds'])==upper
        assert vector(family['cost_lower_coefficients'])==f_lower and vector(family['cost_upper_coefficients'])==f_upper
        # Reconstruct the one frozen sparse selection directly from saved IEEE values.
        raw=sol['lambda']['ineqlin'];assert len(raw)==51601
        quantized=[]
        for value in raw:
            guard();number=float(value);assert isfinite(number)
            n,d=number.as_integer_ratio();rounded=integer_ties_even(F(n*2**40,d))
            quantized.append(F(max(0,rounded),2**40))
        support=sorted(sorted((i for i,y in enumerate(quantized) if y>0),key=lambda i:(-quantized[i],i))[:26])
        assert 0<len(support)<=26 and result['support_row_indices']==support
        assert cert['maximum_support']==26 and [r['index'] for r in cert['rows']]==support
        # Saved primitive coefficients and trace, not any regenerated kernel values.
        erows=read(ROOT/'point_enclosures.json')['rows'];enclosures={}
        for record in erows:
            guard();key=(record['kind'],F(record['u']),F(record['v']))
            assert key not in enclosures and key[0] in ('K','A') and len(record['coefficients'])==14
            kind,u,v=key;assert u>=0 and v>=0 and u+v>0
            assert F(record['target'])==(F(int(u+v<=1)) if kind=='K' else F(0))
            if kind=='A':
                assert u>0 and v>0 and u+v!=R
                if u+v>R:
                    for x,y in ((u,v),(v,u)):
                        if y<R:assert 0<y/(x+y)<R/(x+y)<1 and x-(R-y)>0
            for br in record['coefficients']:assert F(br['lower'])<=F(br['upper'])
            enclosures[key]=record
        for u,v in NEW_POINTS:assert ('A',u,v) in enclosures
        events=[json.loads(s) for s in (ROOT/'calls.jsonl').read_text().splitlines()]
        assert len(events)%2==0 and len(events)//2<=392
        seen=set()
        for i in range(len(events)//2):
            guard();a,b=events[2*i:2*i+2]
            assert a['status']=='started' and b['status']=='completed' and a['call']==b['call']==i+1
            assert all(b[k]==v for k,v in a.items() if k!='status')
            key=(a['kind'],F(a['u']),F(a['v']));col=a['basis_column'];assert 0<=col<14
            assert (key,col) not in seen;seen.add((key,col))
            br=enclosures[key]['coefficients'][col]
            assert all(b[k]==br[k] for k in ('lower','upper','ball'))
        assert len(seen)==14*len(enclosures)
        counts=result['counts'];ncalls=len(seen)
        assert counts==dict(columns_started=ncalls,columns_completed=ncalls,columns_failed=0,columns_interrupted=0,
            LP_requested=1,LP_started=1,LP_completed=1,integral_calls=0,box_model_calls=0)
        assert result['in_flight'] is None
        for j in range(14):
            path=ROOT/'basis.local'/(str(j)+'.json');assert digest(path)==result['basis_sha256'][str(j)]
            data=read(path)
            expected=base if j<2 else dict(alpha=base['alpha'],R=base['R'],h=base['h'],N=80,c='0',a=directions[j-2]['coefficients'])
            assert data==expected
        metadata=read(ROOT/'model_metadata.json')
        assert (metadata['point_rows'],metadata['structural_rows'],metadata['rows'],metadata['variables'])==(51576,25,51601,26)
        assert metadata['model_sha256']==model_sha and len(metadata['row_descriptions'])==51576
        maxabs=[max(abs(l),abs(u)) for l,u in zip(lower,upper)]
        r=f_lower.copy();rhs_weight=F(0)
        for cr in cert['rows']:
            guard();idx=cr['index'];y=F(cr['multiplier']);assert y==quantized[idx] and y>0
            origin=cr['origin']
            if idx<51576:
                if idx<26604:kind='K';u,v=map(F,points[kind][idx]);target=F(int(u+v<=1))
                elif idx<51574:kind='A';u,v=map(F,points[kind][idx-26604]);target=F(0)
                else:kind='A';u,v=NEW_POINTS[idx-51574];target=F(0)
                assert origin['type']=='necessary_point' and desc_identity(origin['point'])==(kind,u,v,target)
                assert desc_identity(metadata['row_descriptions'][idx])==(kind,u,v,target)
                assert origin['enclosure_key']==[kind,str(u),str(v)]
                record=enclosures[kind,u,v]
                alo=[-F(z['upper']) for z in record['coefficients']]+[F(0)]*12
                ahi=[-F(z['lower']) for z in record['coefficients']]+[F(0)]*12
                b=-target;state['point_rows_checked']+=1
            else:
                k=idx-51576;assert origin==dict(type='exact_structure',structural_row=k)
                alo=[F(0)]*26
                if k<24:
                    j=k//2;alo[2+j]=F(1 if k%2==0 else -1);alo[14+j]=F(-1)
                else:assert k==24;alo[0]=F(-1);alo[1]=c1
                ahi=alo.copy();b=F(0);state['structural_rows_checked']+=1
            assert vector(cr['A_lower'])==alo and vector(cr['A_upper'])==ahi and F(cr['b'])==b
            center=[(a+z)/2 for a,z in zip(alo,ahi)];radius=[(z-a)/2 for a,z in zip(alo,ahi)]
            assert all(v>=0 for v in radius)
            error=sum((rr*mm for rr,mm in zip(radius,maxabs)),F(0));relaxed=b+error
            assert vector(cr['A_midpoint'])==center and vector(cr['A_radius'])==radius
            assert F(cr['box_error'])==error and F(cr['relaxed_b'])==relaxed
            rhs_weight+=y*relaxed
            r=[rr+y*a for rr,a in zip(r,center)]
            state['support_rows_checked']+=1;state['coefficients_checked']+=26
        terms=[min(rr*l,rr*u) for rr,l,u in zip(r,lower,upper)]
        bound=sum(terms,F(0))-rhs_weight
        assert vector(cert['residual'])==r and vector(cert['box_terms'])==terms
        assert F(cert['weighted_relaxed_rhs'])==rhs_weight and F(cert['lower_bound'])==bound and F(result['exact_lower_bound'])==bound
        assert F(cert['target'])==F(17,5)
        assert bound>F(8977,2500)>F(17,5)
        assert cert['excludes_cost_at_most_target']==result['excludes_cost_at_most_target']==(bound>F(17,5))
        guard();state.update(status='passed_saved_sparse_dual_audit',complete=True,
            lower_bound=str(bound),target='17/5',strict_comparison_8977_over_2500=True,strict_comparison_17_over_5=True,excludes_cost_at_most_target=bound>F(17,5),
            restricted_scope='Fixed alpha/R, 12 directions, full exact parameter box and c>=s*c1 only; not a lower bound for true C2 or all kernels.',
            primitive_column_values_checked=ncalls,point_enclosures_checked=len(enclosures),trace_events_checked=len(events),
            root_brackets_checked=1,coefficient_error_rule='Ahat*x<=b+sum(delta*maxabs(box)); cost coefficient uses exact root lower for c>=0.',
            result_sha256=digest(ROOT/'result.json'),certificate_sha256=digest(ROOT/'certificate.json'),
            trace_sha256=result['trace_sha256'],verifier_sha256=digest(Path(__file__)))
    except TimeoutError as exc:state.update(status='stopped_30_second_budget',complete=False,reason=str(exc))
    except Exception as exc:state.update(status='failed_saved_sparse_dual_audit',complete=False,reason=f'{type(exc).__name__}: {exc}')
    save();print(json.dumps(state))
    if not state['complete']:raise SystemExit(1)


if __name__=='__main__':main()
