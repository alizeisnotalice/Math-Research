#!/usr/bin/env python3
"""Fixed 12-box Bernstein rows, <=1536 models/300s; import does no experiment.

First 24 psi/S1 models gate the already certified baseline. A failure ends this
stage without any T model, refinement retry, LP or search invocation.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS,MAX_MODELS,PRECISION=300.,1536,192
PARENTS=((7,5),(8,6),(10,6))
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
ROUND10_RESULT_SHA='56b60ccc272d722c534853f2dd460d303649e37a96a79ccfc8f938a3f3113ee5'
CORE_SHA='c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'
CHECK_SHA='0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8'
DEFINITION_SHA='b5d9f5de7259b5c430ef6f0cb6e7e2c75310a41e2731497739e877d2e54cefd5'
REFINEMENT_SHA='3018b6bf86029ea434d820210c97cf8dec25e05827b5e7ae8ed0e5da9879ec2e'


class BudgetStop(Exception):pass


def dump(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    temporary.replace(path)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def exact_point(value):
    """An exact finite Arb point -> exact dyadic Fraction, without display parsing."""
    if not value.is_finite() or not value.is_exact():raise ValueError('Not an exact finite Arb endpoint.')
    mantissa,exponent=value.man_exp();m,e=int(mantissa),int(exponent)
    return Q(m*(1<<e)) if e>=0 else Q(m,1<<(-e))


def rational_interval(value):
    if not value.is_finite():raise ValueError('Nonfinite Arb model coefficient.')
    lower,upper=exact_point(value.lower()),exact_point(value.upper())
    if lower>upper:raise ValueError('Reversed exact outward endpoints.')
    return lower,upper


def main():
    result_path,cells_path=ROOT/'build_result.json',ROOT/'cells.json'
    if result_path.exists() or cells_path.exists() or (ROOT/'models.local').exists():
        raise SystemExit('Existing build evidence: no automatic repetition/resume.')
    start=monotonic();current=None
    counts={k:0 for k in ('models_started','models_completed','models_failed','models_interrupted',
        'APsi_started','APsi_completed','AS_started','AS_completed','Bernstein_conversions_completed',
        'specs_started','specs_completed','baseline_rows_checked','T_models_started','T_models_completed')}
    result={'status':'started','phase':'input_integrity','counts':counts,'in_flight':None,
        'limits':{'seconds_hard_wrapper':300,'Arb_precision_bits':192,'APsi_models':12,
                  'AS_models':1524,'total_models':1536,'baseline_models':24,'T_models_after_gate':1512,
                  'boxes':12,'Bernstein_rows':192,'LP_calls':0},
        'baseline_passed':False,'full_rows_complete':False,'models':[]}
    cells={'schema':'exact_bernstein_cells_v1','status':'building_baseline_only','complete':False,
        'baseline_passed':False,'precision_bits':192,'boxes':[],'rows':[],
        'operator_alpha':'47/40','singular_beta':'47/40','R':'33/20','h':'1/32',
        'local_coordinate_convention':'physical u=h*(p+x),v=h*(q+y); model xi/eta in [-1,1]',
        'formula':'c*(Lpsi-Epsi)+s*(LS1-ES1)+sum(mj*tj-(radj+Ej)*zj)>=0; zj>=abs(tj)'}
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_300s_deadline')
    def save():result['elapsed_seconds']=monotonic()-start;dump(result_path,result)
    def checkpoint():
        dump(cells_path,cells);result['cells_sha256']=sha256(cells_path.read_bytes()).hexdigest();save()
    def verify(path,wanted):
        guard();actual=sha256(path.read_bytes()).hexdigest()
        if actual!=wanted:raise ValueError('Pinned source SHA mismatch: '+path.name)
        result.setdefault('inputs',[]).append({'file':str(path.relative_to(ROOT.parent)),'sha256':actual})
        return actual
    save()
    try:
        baseline_path=ROOT.parent/'round01/scaled_certificate.json'
        verify(baseline_path,BASELINE_SHA);baseline=json.loads(baseline_path.read_text())
        if (Q(baseline['alpha']),Q(baseline['R']),Q(baseline['h']),baseline['N'])!=(Q(47,40),Q(33,20),Q(1,32),80):
            raise ValueError('Fixed baseline analytic parameters changed.')
        c1=Q(baseline['c']);h=Q(baseline['h']);radius=Q(baseline['R'])
        previous_path=ROOT.parent/'round10/result.json';verify(previous_path,ROUND10_RESULT_SHA)
        previous=json.loads(previous_path.read_text());order=previous['direction_order']
        coarse=[(i,j) for i in range(10) for j in range(i+1) if i+j<=9]
        prefix=['absolute_c','absolute_s']+[f'coarse_{i}_{j}' for i,j in coarse]
        if order[:32]!=prefix or len(order)!=128:raise ValueError('Frozen inner-support 128 directions changed.')
        fine=[]
        for name in order[32:]:
            label,i,j=name.split('_');ij=(int(i),int(j))
            if label!='fine' or not 0<=ij[1]<=ij[0] or sum(ij)>48:raise ValueError('Fine support/order malformed.')
            fine.append(ij)
        if len(set(fine))!=96:raise ValueError('Frozen fine class duplication.')
        cells['parameter_order']=order
        cells['LP_order']=order+[f'abs_{name}' for name in order[2:]]
        cells['baseline_sha256']=BASELINE_SHA
        cells['source_direction_result_sha256']=ROUND10_RESULT_SHA
        public=ROOT/'public_src'
        if not (public/'hl_check.py').is_file():public=ROOT.parent/'upstream-v4/checker/independent'
        core_path,check_path=public/'hl_core.py',public/'hl_check.py'
        verify(core_path,CORE_SHA);verify(check_path,CHECK_SHA)
        core=load('hl_core',core_path);checker=load('round12_pinned_hl_check',check_path)
        if checker.Spec is not core.Spec:raise ValueError('Checker imported a different hl_core module.')
        from flint import arb,arb_mat,ctx
        import flint
        ctx.prec=PRECISION
        checker.MC3=arb_mat([[core.A(z) for z in row] for row in checker.power_to_bern_matrix(3,Q(-1),Q(1))])
        checker.MC3T=checker.MC3.transpose()
        result.update(flint_version=flint.__version__,python=sys.version,
            MC3_rebuilt_after_precision=True,cache_policy='Each direction owns a separate AKChecker; no cross-mlist phi or cross-direction D reuse.')
        model_dir=ROOT/'models.local';model_dir.mkdir();spec_dir=ROOT/'specs.local';spec_dir.mkdir()
        half_intervals=((Q(0),Q(1,2)),(Q(1,2),Q(1)))
        for parent in PARENTS:
            p,q=parent
            for xx in half_intervals:
                for yy in half_intervals:
                    box={'id':len(cells['boxes']),'parent':[p,q],'local':[str(x) for x in (*xx,*yy)],
                        'physical':[str(h*(p+xx[0])),str(h*(p+xx[1])),str(h*(q+yy[0])),str(h*(q+yy[1]))],
                        'models':{}}
                    physical=[Q(x) for x in box['physical']]
                    if not (0<physical[2]<=physical[3]<physical[0] and 0<physical[0]+physical[2] and physical[1]+physical[3]<radius):
                        raise ValueError('Frozen half-box outside approved inside branch.')
                    cells['boxes'].append(box)
        if len(cells['boxes'])!=12:raise ValueError('Twelve-box schedule changed.')
        store={}

        def new_checker(path):
            guard();counts['specs_started']+=1;save()
            spec=core.Spec(path);ak=checker.AKChecker(spec)
            counts['specs_completed']+=1;save();guard();return ak

        def model(box,name,ak,method,spec_path,is_T=False):
            nonlocal current
            guard()
            if counts['models_started']>=MAX_MODELS:raise BudgetStop('1536_model_cap')
            current={'box_id':box['id'],'direction':name,'method':method,'status':'started',
                'started_seconds':monotonic()-start,
                'spec_file':'../round01/scaled_certificate.json' if spec_path==baseline_path else str(spec_path.relative_to(ROOT)),
                'spec_sha256':sha256(spec_path.read_bytes()).hexdigest()}
            counts['models_started']+=1
            counts['APsi_started' if method=='APsi_taylor' else 'AS_started']+=1
            if is_T:counts['T_models_started']+=1
            result['in_flight']=dict(current);save()
            try:
                p,q=box['parent'];x0,x1,y0,y1=(Q(x) for x in box['local'])
                if not (0<=x0<x1<=1 and 0<=y0<y1<=1):raise ValueError('Invalid local box.')
                P,E=getattr(ak,method)(p,q,x0,x1,y0,y1)
                guard();el,eu=rational_interval(E)
                if el<0:raise ValueError('Remainder is not certified nonnegative.')
                B=checker.MC3*P*checker.MC3T;bern=[]
                for i in range(4):
                    for j in range(4):
                        guard();lo,hi=rational_interval(B[i,j]);mid=(lo+hi)/2;rad=(hi-lo)/2
                        bern.append({'index':[i,j],'ball_display':B[i,j].str(70),'lower':str(lo),'upper':str(hi),
                            'midpoint':str(mid),'radius':str(rad)})
                current.update(status='completed',power_balls=[[P[i,j].str(70) for j in range(4)] for i in range(4)],
                    bernstein=bern,remainder_ball_display=E.str(70),remainder_upper=str(eu),
                    endpoint_conversion='Exact finite Arb.lower/upper -> man_exp dyadic Fraction, never display parsing.',
                    completed_seconds=monotonic()-start)
                path=model_dir/f'box_{box["id"]:02d}_{name}.json';dump(path,current)
                reference={'file':str(path.relative_to(ROOT)),'sha256':sha256(path.read_bytes()).hexdigest()}
                box['models'][name]=reference;store[box['id'],name]=current
                counts['models_completed']+=1
                counts['APsi_completed' if method=='APsi_taylor' else 'AS_completed']+=1
                counts['Bernstein_conversions_completed']+=1
                if is_T:counts['T_models_completed']+=1
                result['models'].append({'box_id':box['id'],'direction':name,**reference})
                result['in_flight']=None;current=None;save();return store[box['id'],name]
            except BudgetStop:
                counts['models_interrupted']+=1
                if current is not None:current['status']='interrupted';result['in_flight']=dict(current)
                raise
            except Exception:
                counts['models_failed']+=1
                if current is not None:current['status']='failed';result['in_flight']=dict(current)
                raise

        result['phase']='24_baseline_models';baseline_ak=new_checker(baseline_path);checkpoint()
        for box in cells['boxes']:
            model(box,'psi',baseline_ak,'APsi_taylor',baseline_path)
            model(box,'S1',baseline_ak,'AS_model',baseline_path);checkpoint()
        baseline_rows=[]
        for box in cells['boxes']:
            psi,S=store[box['id'],'psi'],store[box['id'],'S1']
            ep,es=Q(psi['remainder_upper']),Q(S['remainder_upper'])
            for ell in range(16):
                guard();lp=Q(psi['bernstein'][ell]['lower'])-ep
                ls=Q(S['bernstein'][ell]['lower'])-es
                slack=c1*lp+ls
                baseline_rows.append({'box_id':box['id'],'bernstein_index':psi['bernstein'][ell]['index'],
                    'c_column':str(lp),'s_column':str(ls),'exact_baseline_slack':str(slack),'passed':slack>=0})
                counts['baseline_rows_checked']+=1
        if len(baseline_rows)!=192:raise ValueError('Baseline row schedule changed.')
        bad=[i for i,row in enumerate(baseline_rows) if not row['passed']]
        cells['baseline_rows']=baseline_rows
        cells['baseline_passed']=not bad;result['baseline_passed']=not bad
        result['baseline_gate']={'rows':192,'bad_indices':bad,'bad_count':len(bad),
            'minimum_slack':str(min(Q(row['exact_baseline_slack']) for row in baseline_rows)),
            'scope':'Conservative Bernstein-row acceptance, not a new evaluation of the true baseline minimum.'}
        checkpoint();guard()
        if bad:
            cells['status']='partial_baseline_only';result.update(status='baseline_enclosure_failed',
                phase='terminal_no_T_models_no_search',reason='Certified baseline not accepted by these fixed sufficient rows; no refinement or row deletion.')
        else:
            result['phase']='1512_T_models_after_baseline_gate';cells['status']='building_T_models';checkpoint()
            definition_path=ROOT/'dependency_round06_search.py'
            if not definition_path.exists():definition_path=ROOT.parent/'round06/search.py'
            refinement_path=ROOT.parent/'round05/refinement_identity.json'
            verify(definition_path,DEFINITION_SHA);verify(refinement_path,REFINEMENT_SHA)
            D=load('round12_refinement_definitions',definition_path);mask=D.refinement(refinement_path)
            maps=[D.refined_direction(ij,mask) for ij in coarse]+[{ij:Q(1)} for ij in fine]
            for name,direction in zip(order[2:],maps):
                guard()
                if not direction or any(i+j>48 for (i,j),v in direction.items() if v):raise ValueError('T support invalid.')
                path=spec_dir/(name+'.json')
                dump(path,{'alpha':'47/40','R':'33/20','h':'1/32','N':80,'c':'0',
                    'a':{f'{i},{j}':str(v) for (i,j),v in sorted(direction.items()) if v}})
                ak=new_checker(path)
                for box in cells['boxes']:model(box,name,ak,'AS_model',path,is_T=True)
                checkpoint()
            result['phase']='assemble_exact_192_rows'
            for bindex,base_row in enumerate(baseline_rows):
                guard();box_id=base_row['box_id'];ell=bindex%16
                row=[Q(base_row['c_column']),Q(base_row['s_column'])];errors=[]
                for name in order[2:]:
                    item=store[box_id,name];entry=item['bernstein'][ell]
                    row.append(Q(entry['midpoint']));errors.append(-(Q(entry['radius'])+Q(item['remainder_upper'])))
                row+=errors
                if len(row)!=254:raise ValueError('Exact row dimension changed.')
                cells['rows'].append({'box_id':box_id,'bernstein_index':base_row['bernstein_index'],
                    'coefficients':[str(x) for x in row],'lower_bound':'0','baseline_slack':base_row['exact_baseline_slack']})
            if len(cells['rows'])!=192 or counts['models_completed']!=1536:raise ValueError('Incomplete fixed model cover.')
            cells.update(status='complete_baseline_and_direction_rows',complete=True)
            result.update(status='complete_exact_cell_rows',full_rows_complete=True,phase='terminal_ready_for_search')
    except BudgetStop as exc:
        result.update(status='interrupted_saved_checkpoint',reason=str(exc));cells['status']='incomplete_interrupted'
    except Exception as exc:
        result.update(status='failed_preparation_or_model',reason=f'{type(exc).__name__}: {exc}');cells['status']='incomplete_failed'
    assert counts['models_started']<=MAX_MODELS
    assert counts['models_started']==counts['models_completed']+counts['models_failed']+counts['models_interrupted']
    assert counts['APsi_started']<=12 and counts['AS_started']<=1524 and counts['T_models_started']<=1512
    result['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();checkpoint()
    print(json.dumps({'status':result['status'],'counts':counts,'baseline_passed':result['baseline_passed'],
        'seconds':result['elapsed_seconds']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
