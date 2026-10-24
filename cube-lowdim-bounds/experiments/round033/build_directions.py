#!/usr/bin/env python3
"""12 frozen signed directions on 1017 saved leaves; <=12204 AS/300s.

No import-time experiment, no baseline/Psi model, geometry reconstruction or LP.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS,MAX_MODELS,PRECISION=300.,12204,192
INDICES=(2,3,4,5,38,39,41,42,62,84,107,108)
NAMES=('coarse_0_0','coarse_1_0','coarse_1_1','coarse_2_0','fine_8_1','fine_8_0',
       'fine_9_1','fine_9_0','fine_5_3','fine_6_4','fine_4_4','fine_5_5')
OLD_COLUMNS=(2,3,4,5,198,199,201,202,222,244,267,268)
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
COVER_SHA='56fa8c066dad18e4acf7e5b3fb3f55cb0c99a4abeec468eac4f3ad0ae7a44f4a'
VERIFICATION_SHA='313374d3cca0f98a8116b0f6eb9fa1c6aa58e1cb2a0801a3364bea574cfcf1c9'
PARAMETER_SHA='fb9bfcc167dd84adfb7293a7691e4581c2d71bd0d38287707f5420f78cce9178'
ROUND10_SHA='56b60ccc272d722c534853f2dd460d303649e37a96a79ccfc8f938a3f3113ee5'
HELPER_SHA='114298a0fc80a9b995ec9cbfb170abcade2e995122ec69f6b7a7b377bb0481ea'
DEFINITION_SHA='b5d9f5de7259b5c430ef6f0cb6e7e2c75310a41e2731497739e877d2e54cefd5'
REFINEMENT_SHA='3018b6bf86029ea434d820210c97cf8dec25e05827b5e7ae8ed0e5da9879ec2e'
CORE_SHA='c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'
CHECK_SHA='0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8'


class BudgetStop(Exception):pass


def dump(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,ensure_ascii=False,allow_nan=False)+'\n');temporary.replace(path)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def main():
    if any((ROOT/name).exists() for name in ('build_result.json','directions.json','rows.json','models.local')):
        raise SystemExit('Existing direction build evidence: no repeat/resume.')
    start=monotonic();current=None;active_stream=None
    counts={key:0 for key in ('models_started','models_completed','models_failed','models_interrupted',
        'AS_calls_started','AS_calls_returned','Bernstein_conversions_completed','baseline_rows_reused',
        'specs_started','specs_completed','directions_completed')}
    result={'status':'started','complete':False,'counts':counts,'in_flight':None,
        'directions_sha256':None,'rows_sha256':None,'completed_direction_files':[],
        'limits':{'hard_seconds':300,'Arb_bits':192,'T_AS_models':12204,'directions':12,'leaves':1017},
        'baseline_model_calls':0,'Psi_model_calls':0,'new_point_values':0,'LP_calls':0}
    directions={'complete':False,'baseline_sha256':BASELINE_SHA,'leaf_cover_sha256':COVER_SHA,
        'parameter_order':['c','s']+list(NAMES)+['abs_'+name for name in NAMES],
        'directions':[],'selection':'Top twelve old absolute coordinate amplitudes, not mass-weighted ranking; old-order output.'}
    rows={'complete':False,'baseline_sha256':BASELINE_SHA,'leaf_cover_sha256':COVER_SHA,
        'direction_names':list(NAMES),'cell_rhs':'baseline_slack/2','rows':[]}
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_300s_deadline')
    def save():result['elapsed_seconds']=monotonic()-start;dump(ROOT/'build_result.json',result)
    def verify(path,wanted):
        guard();actual=sha256(path.read_bytes()).hexdigest()
        if actual!=wanted:raise ValueError('Pinned input SHA mismatch: '+path.name)
        result.setdefault('inputs',[]).append({'file':str(path.relative_to(ROOT.parent)),'sha256':actual})
    save()
    try:
        baseline_path=ROOT.parent/'round01/scaled_certificate.json';verify(baseline_path,BASELINE_SHA)
        baseline=json.loads(baseline_path.read_text());c1=Q(baseline['c'])
        if (Q(baseline['alpha']),Q(baseline['R']),Q(baseline['h']),baseline['N'])!=(Q(47,40),Q(33,20),Q(1,32),80):
            raise ValueError('Fixed analytic baseline parameters changed.')
        cover_path=ROOT.parent/'round13/leaf_cover.json';verify(cover_path,COVER_SHA)
        verification_path=ROOT.parent/'round13/cover_verification.json';verify(verification_path,VERIFICATION_SHA)
        cover=json.loads(cover_path.read_text());verification=json.loads(verification_path.read_text())
        if not cover['complete'] or cover['baseline_sha256']!=BASELINE_SHA or cover['pending'] or cover['unresolved'] or cover['in_flight']:
            raise ValueError('Incomplete or wrong saved baseline cover.')
        if verification['cover_sha256']!=COVER_SHA or verification['status']!='exact_slab_cover_and_all_positive_rows_verified':
            raise ValueError('Saved independent exact cover verification invalid.')
        leaves=cover['leaves']
        if len(leaves)!=1017 or len({leaf['node_id'] for leaf in leaves})!=1017:raise ValueError('Frozen leaf set changed.')
        for leaf in leaves:
            if len(leaf['rows'])!=16:raise ValueError('Missing baseline rows.')
            for entry in leaf['rows']:
                guard();slack=Q(entry['baseline_slack'])
                if c1*Q(entry['c_column'])+Q(entry['s_column'])!=slack or slack<=0:
                    raise ValueError('Saved exact positive baseline row failed arithmetic.')
                rows['rows'].append({'leaf_node_id':leaf['node_id'],'bernstein_index':entry['bernstein_index'],
                    'c_column':entry['c_column'],'s_column':entry['s_column'],'baseline_slack':entry['baseline_slack'],
                    't_midpoints':[],'t_penalties':[]})
                counts['baseline_rows_reused']+=1
        if len(rows['rows'])!=16272:raise ValueError('16272 reused row schedule changed.')
        previous_path=ROOT.parent/'round10/result.json';verify(previous_path,ROUND10_SHA)
        parameter_path=ROOT.parent/'round10/beta_01_floating_parameters.json';verify(parameter_path,PARAMETER_SHA)
        previous=json.loads(previous_path.read_text());parameters=json.loads(parameter_path.read_text())
        order=previous['direction_order'];values=[Q(x) for x in parameters['parameters']]
        if len(order)!=128 or len(values)!=128 or parameters['beta']!='47/40':raise ValueError('Old parameter source changed.')
        ranked=sorted(range(2,128),key=lambda i:(-abs(values[i]),i))[:12]
        if tuple(sorted(ranked))!=INDICES or tuple(order[i] for i in INDICES)!=NAMES:
            raise ValueError('Frozen top-amplitude twelve directions do not match old data.')
        definition_path=ROOT/'dependency_round06_search.py'
        if not definition_path.exists():definition_path=ROOT.parent/'round06/search.py'
        refinement_path=ROOT.parent/'round05/refinement_identity.json'
        verify(definition_path,DEFINITION_SHA);verify(refinement_path,REFINEMENT_SHA)
        D=load('round14_refinement_definitions',definition_path);mask=D.refinement(refinement_path)
        helper_path=ROOT/'dependency_export_cover.py'
        if not helper_path.exists():helper_path=ROOT.parent/'round13/export_cover.py'
        verify(helper_path,HELPER_SHA);H=load('round14_exact_arb_helpers',helper_path)
        public=ROOT/'public_src'
        if not (public/'hl_check.py').is_file():public=ROOT.parent/'upstream-v4/checker/independent'
        verify(public/'hl_core.py',CORE_SHA);verify(public/'hl_check.py',CHECK_SHA)
        core=load('hl_core',public/'hl_core.py');checker=load('round14_pinned_hl_check',public/'hl_check.py')
        if checker.Spec is not core.Spec:raise ValueError('Checker imported another hl_core module.')
        from flint import arb_mat,ctx
        import flint
        ctx.prec=PRECISION
        checker.MC3=arb_mat([[core.A(x) for x in line] for line in checker.power_to_bern_matrix(3,Q(-1),Q(1))])
        checker.MC3T=checker.MC3.transpose()
        result.update(flint_version=flint.__version__,python=sys.version,MC3_rebuilt_after_precision=True,
            cache_policy='One separate Spec/AKChecker per signed direction, never baseline/Psi reconstruction.')
        model_dir=ROOT/'models.local';spec_dir=ROOT/'specs.local';model_dir.mkdir();spec_dir.mkdir()
        for index,name,old_column in zip(INDICES,NAMES,OLD_COLUMNS):
            label,i,j=name.split('_');ij=int(i),int(j)
            if label=='coarse':
                mapping=D.refined_direction(ij,mask);spacing=D.COARSE_Q;bound=Q(1,2)
                expected_column=2+D.COARSE.index(ij)
            else:
                mapping={ij:Q(1)};spacing=D.H_Q;bound=Q(1,4)
                expected_column=192+order[32:].index(name)
            if expected_column!=old_column or any(i+j>48 for (i,j),v in mapping.items() if v):
                raise ValueError('Old column mapping or inner support invalid.')
            mass=spacing**2*len(D.orbit(*ij))/2
            refined_mass=sum((D.H_Q**2*len(D.orbit(*kl))*v/2 for kl,v in mapping.items()),Q(0))
            if refined_mass!=mass:raise ValueError('Exact signed-direction half mass identity failed.')
            item={'name':name,'old_column':old_column,'source_round10_index':index,
                'coefficients':{f'{i},{j}':str(v) for (i,j),v in sorted(mapping.items()) if v},
                'half_mass':str(mass),'bound':str(bound),'old_parameter':str(values[index])}
            directions['directions'].append(item)
        if len(directions['directions'])!=12:raise ValueError('Twelve-column metadata changed.')
        call_keys=set()
        for dindex,item in enumerate(directions['directions']):
            guard();name=item['name'];spec_path=spec_dir/(name+'.json')
            dump(spec_path,{'alpha':'47/40','R':'33/20','h':'1/32','N':80,'c':'0','a':item['coefficients']})
            item['spec_sha256']=sha256(spec_path.read_bytes()).hexdigest()
            counts['specs_started']+=1;save();spec=core.Spec(spec_path);ak=checker.AKChecker(spec)
            counts['specs_completed']+=1;guard()
            path=model_dir/(name+'.jsonl');active_stream=path.open('x')
            for leaf_index,leaf in enumerate(leaves):
                guard();key=(name,leaf['node_id'])
                if key in call_keys:raise ValueError('Duplicate T model call.')
                call_keys.add(key)
                if counts['models_started']>=MAX_MODELS:raise BudgetStop('12204_model_cap')
                current={'direction':name,'direction_index':dindex,'leaf_node_id':leaf['node_id'],
                    'parent':leaf['parent'],'local':leaf['local'],'status':'started'}
                counts['models_started']+=1;counts['AS_calls_started']+=1;result['in_flight']=dict(current);save()
                try:
                    p,q=leaf['parent'];x0,x1,y0,y1=(Q(x) for x in leaf['local'])
                    if not (0<=x0<x1<=1 and 0<=y0<y1<=1):raise ValueError('Invalid saved leaf local geometry.')
                    P,E=ak.AS_model(p,q,x0,x1,y0,y1);counts['AS_calls_returned']+=1;guard()
                    elo,eup=H.interval(E)
                    if elo<0:raise ValueError('Remainder not certified nonnegative; no clamp.')
                    B=checker.MC3*P*checker.MC3T;bern=[];mids=[];penalties=[]
                    for ii in range(4):
                        for jj in range(4):
                            guard();lo,hi=H.interval(B[ii,jj]);mid=(lo+hi)/2;rad=(hi-lo)/2;penalty=rad+eup
                            if penalty<0:raise ValueError('Negative signed error penalty.')
                            bern.append({'index':[ii,jj],'lower':str(lo),'upper':str(hi),
                                'midpoint':str(mid),'radius':str(rad),'penalty':str(penalty),
                                'ball_display':B[ii,jj].str(70)})
                            mids.append(str(mid));penalties.append(str(penalty))
                    current.update(status='completed',spec_sha256=item['spec_sha256'],
                        remainder_upper=str(eup),remainder_ball_display=E.str(70),bernstein=bern,
                        power_balls=[[P[ii,jj].str(70) for jj in range(4)] for ii in range(4)])
                    active_stream.write(json.dumps(current,ensure_ascii=False,allow_nan=False)+'\n');active_stream.flush()
                    for ell,(mid,penalty) in enumerate(zip(mids,penalties)):
                        row=rows['rows'][16*leaf_index+ell]
                        if row['leaf_node_id']!=leaf['node_id'] or row['bernstein_index']!=bern[ell]['index']:
                            raise ValueError('Signed coefficient-to-baseline row alignment changed.')
                        row['t_midpoints'].append(mid);row['t_penalties'].append(penalty)
                    counts['models_completed']+=1;counts['Bernstein_conversions_completed']+=1
                    result['in_flight']=None;current=None;save();guard()
                except BudgetStop:
                    if current is not None:counts['models_interrupted']+=1;current['status']='interrupted';result['in_flight']=dict(current)
                    raise
                except Exception:
                    if current is not None:counts['models_failed']+=1;current['status']='failed';result['in_flight']=dict(current)
                    raise
            active_stream.close();active_stream=None
            item.update(models_file=str(path.relative_to(ROOT)),models_sha256=sha256(path.read_bytes()).hexdigest(),models=1017)
            result['completed_direction_files'].append({'name':name,'file':item['models_file'],'sha256':item['models_sha256']})
            counts['directions_completed']+=1;save();del ak,spec;guard()
        if counts['models_completed']!=12204 or any(len(row['t_midpoints'])!=12 or len(row['t_penalties'])!=12 for row in rows['rows']):
            raise ValueError('Incomplete twelve-direction row set.')
        directions['complete']=True;rows['complete']=True
        result.update(status='complete_twelve_signed_safety_columns',complete=True)
    except BudgetStop as exc:result.update(status='interrupted_saved_checkpoint',reason=str(exc))
    except Exception as exc:result.update(status='failed_preparation_or_model',reason=f'{type(exc).__name__}: {exc}')
    finally:
        if active_stream is not None:active_stream.close()
    assert counts['models_started']<=MAX_MODELS
    assert counts['models_started']==counts['models_completed']+counts['models_failed']+counts['models_interrupted']
    assert counts['AS_calls_started']==counts['models_started'] and counts['AS_calls_returned']<=counts['AS_calls_started']
    dump(ROOT/'directions.json',directions);dump(ROOT/'rows.json',rows)
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        directions_sha256=sha256((ROOT/'directions.json').read_bytes()).hexdigest(),
        rows_sha256=sha256((ROOT/'rows.json').read_bytes()).hexdigest())
    save();print(json.dumps({'status':result['status'],'complete':result['complete'],'counts':counts,
        'seconds':result['elapsed_seconds']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
