#!/usr/bin/env python3
"""Recover only three baseline knot covers; <=1500 boxes/300s, no T/LP.

One AS_model, one constant-psi bound, and only if needed one APsi_taylor per
unique box. All accepted leaves keep exact geometry and sixteen c/s rows.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS,MAX_BOXES,MAX_DEPTH,PRECISION=300.,1500,4,192
PARENTS=((7,5),(8,6),(10,6))
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
HISTORY_SHA='30a79ecad7546f78082d8e5299292396b54e118ad4d8df9ff97bc3df4dc11faf'
CORE_SHA='c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'
CHECK_SHA='0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8'
GRID=tuple(Q(x) for x in ('0','1/256','1/128','1/64','1/32','1/16','1/8','1/4','3/8','1/2',
    '5/8','3/4','7/8','15/16','31/32','63/64','127/128','255/256','1'))


class BudgetStop(Exception):pass


def dump(path,value,compact=False):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=None if compact else 2,ensure_ascii=False,allow_nan=False)+'\n')
    temporary.replace(path)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def exact_point(value):
    if not value.is_finite() or not value.is_exact():raise ValueError('Expected finite exact Arb endpoint.')
    mantissa,exponent=value.man_exp();m,e=int(mantissa),int(exponent)
    return Q(m*(1<<e)) if e>=0 else Q(m,1<<(-e))


def interval(value):
    if not value.is_finite():raise ValueError('Nonfinite model/constant bound.')
    lo,hi=exact_point(value.lower()),exact_point(value.upper())
    if lo>hi:raise ValueError('Reversed exact bounds.')
    return lo,hi


def main():
    result_path,cover_path=ROOT/'result.json',ROOT/'leaf_cover.json'
    if result_path.exists() or cover_path.exists() or (ROOT/'models.local').exists():
        raise SystemExit('Existing cover evidence: no repeat/resume or hidden baseline rerun.')
    start=monotonic();current=None
    counts={key:0 for key in ('boxes_started','boxes_completed','boxes_failed','boxes_interrupted',
        'AS_started','AS_completed','AS_failed','AS_interrupted','const_started','const_completed',
        'const_failed','const_interrupted','Taylor_started','Taylor_completed','Taylor_failed','Taylor_interrupted',
        'Bernstein_conversions_completed','exact_rows_checked','const_accepted','Taylor_accepted','splits')}
    result={'status':'started','counts':counts,'in_flight':None,
        'limits':{'wall_seconds_hard':300,'started_boxes':1500,'max_depth':4,'AS_model_calls':1500,
            'APsi_lower_const_calls':1500,'APsi_taylor_calls':1500,'total_component_calls':4500,
            'Arb_bits':192,'T_models':0,'LP_calls':0},'new_T_models':0,'LP_calls':0,
        'scope':'Reconstruct successful geometry missing from three historical baseline knot-cell summaries only.'}
    cover={'status':'preparing','complete':False,'baseline_sha256':BASELINE_SHA,'alpha':'47/40',
        'R':'33/20','h':'1/32','max_depth':4,'precision_bits':192,'parents':[list(p) for p in PARENTS],
        'grid':[str(x) for x in GRID],'roots':[],'nodes':[],'leaves':[],'pending':[],
        'unresolved':[],'in_flight':None,'parameter_order':['c','s'],
        'endpoint_conversion':'Exact finite Arb.lower/upper.man_exp to dyadic Fraction; display strings diagnostic only.'}
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_300s_deadline')
    def save():result['elapsed_seconds']=monotonic()-start;dump(result_path,result)
    def checkpoint():
        dump(cover_path,cover,compact=True);result['cover_sha256']=sha256(cover_path.read_bytes()).hexdigest();save()
    def verify(path,wanted):
        guard();actual=sha256(path.read_bytes()).hexdigest()
        if actual!=wanted:raise ValueError('Pinned input SHA mismatch: '+path.name)
        result.setdefault('inputs',[]).append({'file':str(path.relative_to(ROOT.parent)),'sha256':actual})
    save()
    try:
        baseline_path=ROOT.parent/'round01/scaled_certificate.json';verify(baseline_path,BASELINE_SHA)
        baseline=json.loads(baseline_path.read_text());c1=Q(baseline['c']);h=Q(baseline['h']);R=Q(baseline['R'])
        if (Q(baseline['alpha']),R,h,baseline['N'])!=(Q(47,40),Q(33,20),Q(1,32),80):
            raise ValueError('Fixed baseline alpha/R/h/N changed.')
        history_path=ROOT.parent/'scaled_second.json';verify(history_path,HISTORY_SHA)
        history=json.loads(history_path.read_text())
        if history['spec_sha256']!=BASELINE_SHA or not history['c4']['ok']:
            raise ValueError('Historical full certificate/baseline source invalid.')
        cover['historical_report_sha256']=HISTORY_SHA;cover['baseline_c']=str(c1)
        public=ROOT/'public_src'
        if not (public/'hl_check.py').is_file():public=ROOT.parent/'upstream-v4/checker/independent'
        core_path,check_path=public/'hl_core.py',public/'hl_check.py'
        verify(core_path,CORE_SHA);verify(check_path,CHECK_SHA)
        core=load('hl_core',core_path);checker=load('round13_pinned_hl_check',check_path)
        if checker.Spec is not core.Spec:raise ValueError('Checker loaded a different hl_core module.')
        if tuple(checker.GRID)!=GRID:raise ValueError('Public graded GRID changed.')
        from flint import arb_mat,ctx
        import flint
        ctx.prec=PRECISION
        checker.MC3=arb_mat([[core.A(z) for z in row] for row in checker.power_to_bern_matrix(3,Q(-1),Q(1))])
        checker.MC3T=checker.MC3.transpose()
        guard();spec=core.Spec(baseline_path);ak=checker.AKChecker(spec);guard()
        model_dir=ROOT/'models.local';model_dir.mkdir()
        result.update(flint_version=flint.__version__,python=sys.version,MC3_rebuilt_after_precision=True)
        geometry_keys=set();call_keys=set()

        def add_node(parent,local,depth,parent_id=None,root_id=None):
            p,q=parent;x0,x1,y0,y1=local
            if not (0<=x0<x1<=1 and 0<=y0<y1<=1):raise ValueError('Bad local box.')
            physical=(h*(p+x0),h*(p+x1),h*(q+y0),h*(q+y1))
            if not (physical[0]>0 and physical[2]>0 and physical[1]+physical[3]<R):
                raise ValueError('Box outside fixed smooth inside branch.')
            key=(*parent,*local)
            if key in geometry_keys:raise ValueError('Duplicate geometric box: no repeated model call.')
            geometry_keys.add(key);node_id=len(cover['nodes'])
            node={'id':node_id,'parent':list(parent),'local':[str(x) for x in local],
                'physical':[str(x) for x in physical],'depth':depth,'parent_id':parent_id,
                'root_id':node_id if root_id is None else root_id,'status':'pending'}
            cover['nodes'].append(node);return node_id

        for parent in PARENTS:
            for i in range(18):
                for j in range(18):
                    guard();node_id=add_node(parent,(GRID[i],GRID[i+1],GRID[j],GRID[j+1]),0)
                    cover['roots'].append(node_id)
        if len(cover['roots'])!=972:raise ValueError('972-root schedule changed.')
        cover['pending']=list(reversed(cover['roots']));cover['status']='reconstructing';checkpoint()

        def component(node,kind,fn):
            guard();key=(node['id'],kind)
            if key in call_keys:raise ValueError('Same component invoked twice on the same box.')
            call_keys.add(key);counts[kind+'_started']+=1
            if counts[kind+'_started']>1500:raise BudgetStop('component_call_cap')
            cover['in_flight']['stage']=kind;result['in_flight']=dict(cover['in_flight']);save()
            returned=False
            try:
                value=fn();returned=True;counts[kind+'_completed']+=1;guard();save();return value
            except BudgetStop:
                # A returned call is completed even if the following clock check ends the box.
                if counts[kind+'_started']>counts[kind+'_completed']+counts[kind+'_failed']+counts[kind+'_interrupted']:
                    counts[kind+'_interrupted']+=1
                raise
            except Exception:
                if not returned:counts[kind+'_failed']+=1
                raise

        def polynomial_record(P,E):
            guard();el,eu=interval(E)
            if el<0:raise ValueError('Remainder is not nonnegative.')
            B=checker.MC3*P*checker.MC3T;bern=[]
            for i in range(4):
                for j in range(4):
                    guard();lo,hi=interval(B[i,j])
                    bern.append({'index':[i,j],'lower':str(lo),'upper':str(hi),'ball_display':B[i,j].str(70)})
            counts['Bernstein_conversions_completed']+=1
            return {'power_balls':[[P[i,j].str(70) for j in range(4)] for i in range(4)],
                'bernstein':bern,'remainder_upper':str(eu),'remainder_ball_display':E.str(70)}

        def rows_for(bundle,mode):
            AS=bundle['AS'];es=Q(AS['remainder_upper']);rows=[]
            for ell in range(16):
                guard();scol=Q(AS['bernstein'][ell]['lower'])-es
                if mode=='const':ccol=Q(bundle['const']['lower'])
                else:ccol=Q(bundle['Taylor']['bernstein'][ell]['lower'])-Q(bundle['Taylor']['remainder_upper'])
                slack=c1*ccol+scol;counts['exact_rows_checked']+=1
                rows.append({'bernstein_index':AS['bernstein'][ell]['index'],'c_column':str(ccol),
                    's_column':str(scol),'baseline_slack':str(slack)})
            return rows

        capped=False
        while cover['pending']:
            guard()
            if counts['boxes_started']>=MAX_BOXES:capped=True;break
            node_id=cover['pending'].pop();node=cover['nodes'][node_id];current=node
            counts['boxes_started']+=1;node['status']='processing'
            cover['in_flight']={'node_id':node_id,'stage':'before_AS'};result['in_flight']=dict(cover['in_flight']);checkpoint()
            p,q=node['parent'];x0,x1,y0,y1=(Q(x) for x in node['local'])
            bundle={'node_id':node_id,'parent':node['parent'],'local':node['local'],'physical':node['physical'],
                'depth':node['depth'],'baseline_sha256':BASELINE_SHA,'alpha':'47/40','R':'33/20','h':'1/32'}
            P,E=component(node,'AS',lambda:ak.AS_model(p,q,x0,x1,y0,y1))
            bundle['AS']=polynomial_record(P,E)
            u0,u1,v0,v1=(Q(x) for x in node['physical'])
            constant=component(node,'const',lambda:spec.APsi_lower_const(u0,v0,u1,v1))
            lo,hi=interval(constant);bundle['const']={'lower':str(lo),'upper':str(hi),'ball_display':constant.str(70)}
            rows=rows_for(bundle,'const');mode='const'
            if not all(Q(row['baseline_slack'])>0 for row in rows):
                PP,EP=component(node,'Taylor',lambda:ak.APsi_taylor(p,q,x0,x1,y0,y1))
                bundle['Taylor']=polynomial_record(PP,EP)
                rows=rows_for(bundle,'Taylor');mode='Taylor'
            accepted=all(Q(row['baseline_slack'])>0 for row in rows)
            bundle.update(mode=mode,rows=rows,accepted=accepted)
            path=model_dir/f'node_{node_id:05d}.json'
            if path.exists():raise ValueError('Immutable node model already exists.')
            dump(path,bundle);model_sha=sha256(path.read_bytes()).hexdigest()
            node.update(model_file=str(path.relative_to(ROOT)),model_sha256=model_sha)
            if accepted:
                node.update(status='accepted',mode=mode);counts[mode+'_accepted']+=1
                cover['leaves'].append({'node_id':node_id,'parent':node['parent'],'local':node['local'],
                    'physical':node['physical'],'depth':node['depth'],'mode':mode,'rows':rows,
                    'model_file':node['model_file'],'model_sha256':model_sha})
            elif node['depth']<MAX_DEPTH:
                xm,ym=(x0+x1)/2,(y0+y1)/2
                children=[add_node((p,q),local,node['depth']+1,node_id,node['root_id']) for local in (
                    (x0,xm,y0,ym),(xm,x1,y0,ym),(x0,xm,ym,y1),(xm,x1,ym,y1))]
                node.update(status='split',children=children,split_method='exact_midpoint')
                cover['pending'].extend(reversed(children));counts['splits']+=1
            else:node['status']='unresolved_depth4';cover['unresolved'].append(node_id)
            counts['boxes_completed']+=1;cover['in_flight']=None;result['in_flight']=None;current=None;checkpoint()
        if capped:result['status']='box_cap_saved_pending';cover['status']='incomplete_box_cap'
        elif cover['unresolved']:result['status']='depth_cap_unresolved';cover['status']='incomplete_depth4'
        else:
            for parent in PARENTS:
                area=sum(((Q(leaf['local'][1])-Q(leaf['local'][0]))*(Q(leaf['local'][3])-Q(leaf['local'][2]))
                    for leaf in cover['leaves'] if tuple(leaf['parent'])==parent),Q(0))
                if area!=1:raise ValueError('Accepted area does not exactly cover one whole parent.')
            cover.update(status='geometry_complete',complete=True);result['status']='geometry_complete'
    except BudgetStop as exc:
        result.update(status='interrupted_saved_checkpoint',reason=str(exc));cover['status']='incomplete_interrupted'
        if current is not None:
            current['status']='interrupted';counts['boxes_interrupted']+=1
    except Exception as exc:
        result.update(status='failed_preparation_or_box',reason=f'{type(exc).__name__}: {exc}');cover['status']='incomplete_failed'
        if current is not None:current['status']='failed';counts['boxes_failed']+=1
    assert counts['boxes_started']<=MAX_BOXES
    assert counts['boxes_started']==counts['boxes_completed']+counts['boxes_failed']+counts['boxes_interrupted']
    for kind in ('AS','const','Taylor'):
        assert counts[kind+'_started']<=1500
        assert counts[kind+'_started']==sum(counts[kind+'_'+key] for key in ('completed','failed','interrupted'))
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),accepted_leaves=len(cover['leaves']),
        pending_count=len(cover['pending']),unresolved_count=len(cover['unresolved']),complete=cover['complete'])
    checkpoint();print(json.dumps({'status':result['status'],'counts':counts,'accepted_leaves':len(cover['leaves']),
        'seconds':result['elapsed_seconds']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
