#!/usr/bin/env python3
"""Continue only the frozen (52,1) pending tree, without recomputing old models.

New box: one AS, eight Hint calls, then four algebraic proof pieces / 64 rows. At the
first depth-four rejection stop. Numerical events are append-only; no call list
is repeatedly serialized into result.json. No const/T/point/LP/full checker.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from posixpath import normpath
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS,MAX_BOXES,MAX_DEPTH,PRECISION=300.,1500,4,192
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
COVER_SHA='e78b2bf4e6f5e00974a180cb85dfa9044ed42939930c4add48078e6c4064c8e7'
MODEL663_SHA='8b614d16e19e8dddd9cb16793595b99aff3acf680795710c50ced9c3f76c108e'
PIECES663_SHA='1921a3b1e009005ecc5e3bea6309cedc8932e5f8b2862401232974d27eef1581'
ALGEBRA_SHA='fd16ccf0309a4ae4627049ff730b2bb233be3948013d25cd514824febcc107b2'
FIXED_SHA='4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05'
CORE_SHA='c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'
CHECK_SHA='0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8'


class BudgetStop(Exception):pass


def dump(path,value,compact=False):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=None if compact else 2,ensure_ascii=False,allow_nan=False)+'\n')
    tmp.replace(path)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module);return module


def endpoint(value):
    if not value.is_finite() or not value.is_exact():raise ValueError('Expected finite exact endpoint.')
    m,e=value.man_exp();m,e=int(m),int(e)
    return Q(m*(1<<e)) if e>=0 else Q(m,1<<(-e))


def bounds(value):
    if not value.is_finite():raise ValueError('Nonfinite ball.')
    lo,hi=endpoint(value.lower()),endpoint(value.upper())
    if lo>hi:raise ValueError('Reversed interval.')
    return lo,hi


def main():
    if any((ROOT/name).exists() for name in ('result.json','leaf_cover.json','models.local','calls.jsonl')):
        raise SystemExit('Existing evidence: no repeat/resume.')
    start=monotonic();current=None;cover=None
    result_path,cover_path,trace_path=ROOT/'result.json',ROOT/'leaf_cover.json',ROOT/'calls.jsonl'
    counts={kind:{key:0 for key in ('started','completed','failed','interrupted')}
        for kind in ('AS','Psi','Hint','acb_integral')}
    boxes={key:0 for key in ('started','completed','failed','interrupted','accepted','splits')}
    proof_counts={kind:{key:0 for key in ('started','completed','failed','interrupted')} for kind in ('pieces','rows')}
    result={'proof_counts':proof_counts,'proof_in_flight':None,'status':'started','counts':counts,'new_boxes':boxes,'calls_started_total':0,'in_flight':[],
        'limits':{'seconds_hard':300,'Arb_bits':192,'new_started_boxes':1500,'max_original_root_depth':4,
            'AS':1500,'Psi':1500,'Hint':12000,'acb_integral':12000,'proof_pieces':6000,'new_rows':96000},
        'full_absolute_Psi_calls':0,'new_const_calls':0,'new_T_calls':0,'new_point_calls':0,'LP_calls':0,'full_checker_calls':0,
        'inherited_work':{'round18_AS_calls':385,'round18_const_calls':385,'round19_Psi_models':1,'round20_AS_calls':61,'round20_Psi_models':61,
            'round21_Hint_calls':2,'round21_integral_calls':2,'round22_AS_calls':2,'round22_Psi_models':2,
            'round22_Hint_calls':16,'round22_integral_calls':16,'round23_new_qualification_rows':64,
            'initial_accepted_tree_leaves':363,'initial_accepted_proof_regions':366,'initial_pending':216},
        'call_completion':'Finite numeric return; nonfinite/exception failed. Saved-component counts are separate.'}
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_300s_deadline')
    def save():result['elapsed_seconds']=monotonic()-start;dump(result_path,result)
    def checkpoint():
        if cover is not None:
            cover['accepted_tree_leaves']=len(cover['leaves'])
            cover['accepted_proof_regions']=sum(leaf.get('proof_region_count',1) for leaf in cover['leaves'])
            result['accepted_tree_leaves']=cover['accepted_tree_leaves']
            result['accepted_proof_regions']=cover['accepted_proof_regions']
            dump(cover_path,cover,True);result['cover_sha256']=sha256(cover_path.read_bytes()).hexdigest()
        save()
    def verify(relative,wanted):
        guard();path=ROOT.parent/relative;actual=sha256(path.read_bytes()).hexdigest()
        if actual!=wanted:raise ValueError('Pinned SHA mismatch: '+relative)
        result.setdefault('inputs',[]).append({'file':relative,'sha256':actual});return path
    save()
    trace=None;M=None;kernel=None;original_hint=None;real_acb=None
    try:
        old_cover=verify('round22/leaf_cover.json',COVER_SHA)
        old663_path=verify('round22/models.local/node_00663.json',MODEL663_SHA)
        pieces663_path=verify('round23/quadrants.json',PIECES663_SHA)
        algebra_path=verify('round23/quadrants_once.py',ALGEBRA_SHA)
        baseline_path=verify('round01/scaled_certificate.json',BASELINE_SHA)
        fixed=verify('round17/fixed_checker/hlcert.py',FIXED_SHA)
        public=ROOT/'public_src'
        if not (public/'hl_check.py').is_file():public=ROOT.parent/'upstream-v4/checker/independent'
        core_path,check_path=public/'hl_core.py',public/'hl_check.py'
        for path,wanted in ((core_path,CORE_SHA),(check_path,CHECK_SHA)):
            guard()
            if sha256(path.read_bytes()).hexdigest()!=wanted:raise ValueError('Pinned independent source mismatch.')
            result.setdefault('inputs',[]).append({'file':str(path.relative_to(ROOT.parent)),'sha256':wanted})
        cover=json.loads(old_cover.read_text());old663=json.loads(old663_path.read_text());pieces663=json.loads(pieces663_path.read_text())
        base=json.loads(baseline_path.read_text());c1,h,R=map(Q,(base['c'],base['h'],base['R']))
        if (Q(base['alpha']),R,h,base['N'])!=(Q(47,40),Q(33,20),Q(1,32),80) or c1<=0:
            raise ValueError('Fixed baseline parameters changed.')
        if cover['baseline_sha256']!=BASELINE_SHA or cover['parents']!=[[52,1]] or \
                len(cover['roots'])!=324 or len(cover['leaves'])!=362 or len(cover['pending'])!=216 or \
                cover['unresolved']!=[663] or cover['in_flight'] is not None or cover['max_depth']!=4:
            raise ValueError('Frozen continuation state changed.')
        old_nodes=cover['nodes'];old_nodes_count=len(old_nodes)
        if any(node['id']!=i for i,node in enumerate(old_nodes)):raise ValueError('Noncanonical node IDs.')
        if len(set(cover['pending']))!=216 or any(old_nodes[i]['status']!='pending' for i in cover['pending']):
            raise ValueError('Invalid inherited pending list.')
        for node in old_nodes:
            guard()
            # Archived final models are mandatory; local partial checkpoints are optional provenance.
            if node['status']!='pending' and 'model_file' not in node:
                raise ValueError('Completed historical node lacks a final model.')
            if 'model_file' in node:
                old_relative=node['model_file'];old_path=old_cover.parent/old_relative
                if sha256(old_path.read_bytes()).hexdigest()!=node['model_sha256']:
                    raise ValueError('Old final model hash mismatch.')
                node['model_file']=normpath('../round22/'+old_relative)
            if 'partial_model_file' in node:
                node['partial_model_file']=normpath('../round22/'+node['partial_model_file'])
                node['historical_partial_optional']=True
            for fk,sk in (('supplemental_Psi_file','supplemental_Psi_sha256'),('supplemental_signed_file','supplemental_signed_sha256'),('supplemental_piece_file','supplemental_piece_sha256')):
                if fk in node:
                    p=old_cover.parent/node[fk]
                    if sha256(p.read_bytes()).hexdigest()!=node[sk]:raise ValueError('Historical supplemental hash mismatch.')
                    node[fk]=normpath('../round22/'+node[fk])
            node['inherited_from_cover']='round22'
        for leaf in cover['leaves']:
            guard();node=old_nodes[leaf['node_id']]
            if node['status']!='accepted' or leaf['model_sha256']!=node['model_sha256']:
                raise ValueError('Inherited leaf membership changed.')
            leaf['model_file']=node['model_file'];leaf['inherited_from_cover']='round22';leaf['proof_region_count']=1
            for fk in ('supplemental_Psi_file','supplemental_signed_file','supplemental_piece_file'):
                if fk in node:leaf[fk]=node[fk]
            for row in leaf['rows']:
                if Q(row['baseline_slack'])!=c1*Q(row['c_column'])+Q(row['s_column']) or Q(row['baseline_slack'])<=0:
                    raise ValueError('Inherited accepted row invalid.')
        node663=old_nodes[663]
        if (node663['local'],node663['physical'])!=(old663['local'],old663['physical']) or \
                (pieces663['local'],pieces663['physical'])!=(old663['local'],old663['physical']) or \
                pieces663['source_model_sha256']!=MODEL663_SHA or pieces663['source_cover_sha256']!=COVER_SHA or \
                pieces663['baseline_sha256']!=BASELINE_SHA or pieces663['source_tree_depth']!=4 or \
                pieces663['source_tree_changed'] is not False or len(pieces663['pieces'])!=4:
            raise ValueError('Node663 piecewise provenance mismatch.')
        proof_pieces=[]
        for piece in pieces663['pieces']:
            guard()
            if piece['proof_piece'] is not True or piece['source_tree_depth']!=4 or len(piece['rows'])!=16:
                raise ValueError('Invalid saved proof piece.')
            for row in piece['rows']:
                if Q(row['baseline_slack'])!=c1*Q(row['c_column'])+Q(row['s_column']) or Q(row['baseline_slack'])<=0:
                    raise ValueError('Saved node663 accepted row invalid.')
            copied=dict(piece);copied.update(source_tree_node_id=663,origin='round23')
            proof_pieces.append(copied)
        node663.update(status='accepted',mode='outside_quadrants_reused',supplemental_piece_file='../round23/quadrants.json',
            supplemental_piece_sha256=PIECES663_SHA,previous_status='unresolved_depth4',proof_region_count=4)
        cover['leaves'].append({'node_id':663,'parent':[52,1],'local':node663['local'],'physical':node663['physical'],
            'depth':4,'mode':'outside_quadrants_reused','proof_pieces':proof_pieces,'proof_region_count':4,
            'model_file':node663['model_file'],'model_sha256':MODEL663_SHA,
            'supplemental_piece_file':'../round23/quadrants.json','supplemental_piece_sha256':PIECES663_SHA,
            'origin':'round22_AS_T1_plus_round23_pieces'})
        cover.update(status='extending',complete=False,unresolved=[],in_flight=None,branch='strict_exterior_leaf_internal_quadrants',
            source_cover_sha256=COVER_SHA,source_pieces663_sha256=PIECES663_SHA,initial_inherited_nodes=old_nodes_count,
            new_const_calls=0,new_T_models=0,LP_calls=0,one_sided_new_Psi=True,
            accepted_tree_leaves=363,accepted_proof_regions=366)
        algebra=load('round24_frozen_quadrant_algebra',algebra_path)
        checkpoint()
        core=load('hl_core',core_path);checker=load('round24_hl_check',check_path);M=load('round24_fixed_hlcert',fixed)
        from flint import arb_mat,ctx
        import flint
        ctx.prec=PRECISION
        checker.MC3=arb_mat([[core.A(z) for z in row] for row in checker.power_to_bern_matrix(3,Q(-1),Q(1))])
        checker.MC3T=checker.MC3.transpose()
        if checker.Spec is not core.Spec:raise ValueError('Independent core import mismatch.')
        guard();spec=core.Spec(baseline_path);ak=checker.AKChecker(spec);kernel=M.Kern(base);guard()
        result.update(flint_version=flint.__version__,python=sys.version,MC3_rebuilt_after_precision=True)
        model_dir=ROOT/'models.local';model_dir.mkdir();trace=trace_path.open('x')
        seen={(*node['parent'],*node['local']) for node in old_nodes};component_keys=set()
        limits={'AS':1500,'Psi':1500,'Hint':12000,'acb_integral':12000};local_calls={}
        active_hint=[]

        def emit(event):trace.write(json.dumps(event,ensure_ascii=False,allow_nan=False)+'\n');trace.flush()
        def recorded(kind,fn,finite,metadata=None):
            guard()
            if counts[kind]['started']>=limits[kind]:raise BudgetStop(kind+'_call_cap')
            if kind in ('AS','Psi'):
                key=(current['id'],kind)
                if key in component_keys:raise ValueError('Repeated box component.')
                component_keys.add(key)
            else:
                if local_calls[kind]>=8:raise BudgetStop('per_box_'+kind+'_cap')
                local_calls[kind]+=1
            cid=result['calls_started_total'];result['calls_started_total']+=1;ts=monotonic()
            event={'id':cid,'kind':kind,'node_id':current['id'],'status':'started',**(metadata or {})}
            counts[kind]['started']+=1;result['in_flight'].append(event.copy());emit(event);save()
            try:
                value=fn()
                if not finite(value):raise ValueError('Nonfinite '+kind+' return.')
                counts[kind]['completed']+=1;event['status']='completed';return value
            except BudgetStop as exc:
                counts[kind]['interrupted']+=1;event.update(status='interrupted',reason=str(exc));raise
            except Exception as exc:
                counts[kind]['failed']+=1;event.update(status='failed',reason=f'{type(exc).__name__}: {exc}');raise
            finally:
                event['elapsed_seconds']=monotonic()-ts
                result['in_flight']=[e for e in result['in_flight'] if e['id']!=cid]
                result['last_event']=event;emit(event);save()
                if event['status']=='completed':guard()

        real_acb=M.acb;original_hint=kernel._Hint
        class AcbProxy:
            def __call__(self,*args,**kwargs):return real_acb(*args,**kwargs)
            def integral(self,*args,**kwargs):
                return recorded('acb_integral',lambda:real_acb.integral(*args,**kwargs),lambda v:v.is_finite(),
                    {'Hint_id':active_hint[-1] if active_hint else None})
        def hint(a,b,kind):
            hid=counts['Hint']['started'];active_hint.append(hid)
            try:return recorded('Hint',lambda:original_hint(a,b,kind),lambda v:v.is_finite(),
                {'Hint_id':hid,'a':str(a),'b':str(b),'derivative_kind':kind})
            finally:active_hint.pop()
        M.acb=AcbProxy();kernel._Hint=hint

        def poly_record(P,E):
            guard();el,eu=bounds(E)
            if el<0:raise ValueError('Negative AS remainder.')
            B=checker.MC3*P*checker.MC3T
            bern=[]
            for i in range(4):
                for j in range(4):
                    guard();lo,hi=bounds(B[i,j]);bern.append({'index':[i,j],'lower':str(lo),'upper':str(hi)})
            return {'bernstein':bern,'remainder_upper':str(eu),
                'power_balls':[[P[i,j].str(70) for j in range(4)] for i in range(4)]}
        def ball_record(value):
            lo,hi=bounds(value)
            return {'lower':str(lo),'upper':str(hi),'display':value.str(70)}
        def signed_Psi(u0,u1,v0,v1):
            guard();uc,vc=M.Q((u0+u1)/2),M.Q((v0+v1)/2)
            f0=kernel._Hint(uc,vc,'H')+kernel._Hint(vc,uc,'H')
            fu=kernel._Hint(uc,vc,'Ha')+kernel._Hint(vc,uc,'Hb')
            fv=kernel._Hint(uc,vc,'Hb')+kernel._Hint(vc,uc,'Ha')
            mixed=[]
            for a,b in ((u0,v0),(v0,u0)):
                value=kernel._Hint(M.Q(a),M.Q(b),'Hab')
                br=ball_record(value);_,au=bounds(abs(value))
                mixed.append({'a':str(a),'b':str(b),'kind':'Hab',**br,'absolute_upper':str(au)})
            muv=sum((Q(row['absolute_upper']) for row in mixed),Q(0))
            du,dv=(u1-u0)/2,(v1-v0)/2;dua,dva=M.arb(M.Q(du)),M.arb(M.Q(dv))
            B1=[[f0+fu*(-dua if i==0 else dua)+fv*(-dva if j==0 else dva) for j in range(2)] for i in range(2)]
            details={'one_sided_lower_deduction':str(muv*du*dv),'Muv_upper':str(muv),'du':str(du),'dv':str(dv),
                'mixed_derivatives':mixed,'center_values':{'F0':ball_record(f0),'Fu':ball_record(fu),'Fv':ball_record(fv)}}
            return M._elev2(B1,3,3),details
        def psi_record(B,details):
            guard()
            if Q(details['one_sided_lower_deduction'])<0:raise ValueError('Negative signed deduction.')
            bern=[]
            for i in range(4):
                for j in range(4):
                    guard();lo,hi=bounds(B[i][j]);bern.append({'index':[i,j],'lower':str(lo),'upper':str(hi)})
            return {'bernstein':bern,**details,'unscaled_by_c':True,'Bernstein_degree':[3,3],'one_sided_only':True}
        def finite_AS(value):
            P,E=value
            return E.is_finite() and all(P[i,j].is_finite() for i in range(4) for j in range(4))
        def finite_Psi(value):
            B,details=value
            return B is not None and len(B)==4 and Q(details['one_sided_lower_deduction'])>=0 and                 all(len(row)==4 and all(x.is_finite() for x in row) for row in B)
        def add_node(local,parent_node):
            if parent_node['depth']>=MAX_DEPTH:raise ValueError('Forbidden deeper split.')
            key=(52,1,*map(str,local))
            if key in seen:raise ValueError('Repeated geometric node.')
            seen.add(key);x0,x1,y0,y1=local
            if not (0<=x0<x1<=1 and 0<=y0<y1<=1):raise ValueError('Invalid local box.')
            node={'id':len(cover['nodes']),'parent':[52,1],'local':list(map(str,local)),
                'physical':list(map(str,(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1)))),
                'depth':parent_node['depth']+1,'parent_id':parent_node['id'],'root_id':parent_node['root_id'],
                'status':'pending','origin':'round24'}
            cover['nodes'].append(node);return node['id']

        capped=False;depth_failed=False
        while cover['pending']:
            guard()
            if boxes['started']>=MAX_BOXES:capped=True;break
            nid=cover['pending'].pop();current=cover['nodes'][nid];boxes['started']+=1
            if current['status']!='pending' or current['parent']!=[52,1] or not 0<=current['depth']<=4:
                raise ValueError('Unapproved inherited/new pending node.')
            current['status']='processing';local_calls={'Hint':0,'acb_integral':0}
            cover['in_flight']={'node_id':nid,'stage':'AS'};checkpoint()
            x0,x1,y0,y1=map(Q,current['local']);u0,u1,v0,v1=map(Q,current['physical'])
            if (u0,u1,v0,v1)!=(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1)) or \
                    not (u0>0 and v0>0 and u0+v0>R):raise ValueError('Exterior geometry failed.')
            bundle={'node_id':nid,'parent':[52,1],'local':current['local'],'physical':current['physical'],
                'depth':current['depth'],'baseline_sha256':BASELINE_SHA,'fixed_hlcert_sha256':FIXED_SHA,
                'independent_core_sha256':CORE_SHA,'independent_check_sha256':CHECK_SHA}
            P,E=recorded('AS',lambda:ak.AS_model(52,1,x0,x1,y0,y1),finite_AS)
            bundle['AS']=poly_record(P,E)
            partial=model_dir/f'node_{nid:05d}.partial.json';dump(partial,bundle)
            current.update(partial_model_file=str(partial.relative_to(ROOT)),partial_model_sha256=sha256(partial.read_bytes()).hexdigest())
            result['AS_models_saved']=result.get('AS_models_saved',0)+1
            cover['in_flight']['stage']='Psi';checkpoint();guard()
            B,details=recorded('Psi',lambda:signed_Psi(u0,u1,v0,v1),finite_Psi)
            bundle['Psi']=psi_record(B,details);dump(partial,bundle)
            current['partial_model_sha256']=sha256(partial.read_bytes()).hexdigest()
            result['Psi_models_saved']=result.get('Psi_models_saved',0)+1
            cover['in_flight']['stage']='rows';checkpoint();guard()
            es=Q(bundle['AS']['remainder_upper']);D=Q(bundle['Psi']['one_sided_lower_deduction'])
            muv=Q(bundle['Psi']['Muv_upper']);uc,vc=(u0+u1)/2,(v0+v1)/2;xc,yc=(x0+x1)/2,(y0+y1)/2
            t1=[[Q(bundle['Psi']['bernstein'][4*i+j]['lower']) for j in range(4)] for i in range(4)]
            AS=[[Q(bundle['AS']['bernstein'][4*i+j]['lower']) for j in range(4)] for i in range(4)]
            rt,rs=algebra.tensor_quadrants(t1),algebra.tensor_quadrants(AS)
            pu,pv=((u0,uc),(uc,u1)),((v0,vc),(vc,v1))
            lx,ly=((x0,xc),(xc,x1)),((y0,yc),(yc,y1));pieces=[]
            for pid,(iu,iv) in enumerate(((0,0),(1,0),(0,1),(1,1))):
                guard()
                if proof_counts['pieces']['started']>=6000:raise BudgetStop('proof_piece_cap')
                proof_counts['pieces']['started']+=1
                result['proof_in_flight']={'node_id':nid,'piece_id':pid,'stage':'restriction'}
                a,b=pu[iu];c,d=pv[iv];same_sign=iu==iv
                q2=[[Q(0) if same_sign else muv*(u-uc)*(v-vc) for v in (c,d)] for u in (a,b)]
                corr=algebra.elevate_bilinear(q2)
                piece={'piece_id':pid,'proof_piece':True,'source_tree_node_id':nid,'source_tree_depth':current['depth'],
                    'physical':list(map(str,(a,b,c,d))),'local':list(map(str,(*lx[iu],*ly[iv]))),
                    'u_sign':'-' if iu==0 else '+','v_sign':'-' if iv==0 else '+','same_sign':same_sign,
                    'restricted_T1_lower':algebra.strings(rt[(iu,iv)]),'restricted_AS_lower':algebra.strings(rs[(iu,iv)]),
                    'bilinear_Bernstein_degree1':algebra.strings(q2),'bilinear_Bernstein_degree3':algebra.strings(corr),'rows':[]}
                for i in range(4):
                    for j in range(4):
                        guard()
                        if proof_counts['rows']['started']>=96000:raise BudgetStop('qualification_row_cap')
                        proof_counts['rows']['started']+=1
                        result['proof_in_flight']={'node_id':nid,'piece_id':pid,'stage':'row','index':[i,j]}
                        tt,ss=rt[(iu,iv)][i][j],rs[(iu,iv)][i][j]-es;cc=tt+corr[i][j]
                        slack=c1*cc+ss;uniform_slack=c1*(tt-D)+ss;gain=slack-uniform_slack
                        if gain!=c1*(D+corr[i][j]) or gain<0:raise ValueError('Piece signed dominance mismatch.')
                        piece['rows'].append({'bernstein_index':[i,j],'c_column':str(cc),'s_column':str(ss),
                            'baseline_slack':str(slack),'restricted_uniform_slack':str(uniform_slack),'gain':str(gain)})
                        proof_counts['rows']['completed']+=1
                piece['accepted']=all(Q(row['baseline_slack'])>0 for row in piece['rows'])
                pieces.append(piece);proof_counts['pieces']['completed']+=1;result['proof_in_flight']=None
            accepted=all(piece['accepted'] for piece in pieces)
            bundle.update(mode='outside_quadrants',proof_pieces=pieces,proof_region_count=4,accepted=accepted,
                proof_restriction_source_sha256=ALGEBRA_SHA)
            out=model_dir/f'node_{nid:05d}.json'
            if out.exists():raise ValueError('Immutable model collision.')
            dump(out,bundle);current.update(model_file=str(out.relative_to(ROOT)),model_sha256=sha256(out.read_bytes()).hexdigest())
            if accepted:
                current.update(status='accepted',mode='outside_quadrants',proof_region_count=4);boxes['accepted']+=1
                cover['leaves'].append({key:current[key] for key in ('parent','local','physical','depth','mode','model_file','model_sha256')})
                cover['leaves'][-1].update(node_id=nid,proof_pieces=pieces,proof_region_count=4,origin='round24',
                    proof_file=current['model_file'],proof_sha256=current['model_sha256'])
            elif current['depth']<MAX_DEPTH:
                xm,ym=(x0+x1)/2,(y0+y1)/2
                children=[add_node(local,current) for local in ((x0,xm,y0,ym),(xm,x1,y0,ym),(x0,xm,ym,y1),(xm,x1,ym,y1))]
                current.update(status='split',children=children,split_method='exact_midpoint')
                cover['pending'].extend(reversed(children));boxes['splits']+=1
            else:
                current['status']='unresolved_depth4';cover['unresolved'].append(nid);depth_failed=True
            boxes['completed']+=1;current=None;cover['in_flight']=None;checkpoint()
            if depth_failed:break
        if capped:result['status']='box_cap_saved_pending';cover['status']='incomplete_box_cap'
        elif cover['unresolved']:result['status']='depth_cap_unresolved';cover['status']='incomplete_depth4'
        else:
            if cover['pending'] or cover['in_flight']:raise ValueError('Incomplete terminal state.')
            area=sum(((Q(leaf['local'][1])-Q(leaf['local'][0]))*(Q(leaf['local'][3])-Q(leaf['local'][2])) for leaf in cover['leaves']),Q(0))
            if area!=1:raise ValueError('Accepted local area is not exactly one.')
            cover.update(status='geometry_complete',complete=True);result['status']='geometry_complete'
    except BudgetStop as exc:
        result.update(status='interrupted_saved',reason=str(exc))
        for cc in proof_counts.values():cc['interrupted']+=cc['started']-cc['completed']
        if cover is not None:cover['status']='incomplete_interrupted'
        if current is not None:current['status']='interrupted';boxes['interrupted']+=1
    except Exception as exc:
        result.update(status='failed',reason=f'{type(exc).__name__}: {exc}')
        for cc in proof_counts.values():cc['failed']+=cc['started']-cc['completed']
        if cover is not None:cover['status']='incomplete_failed'
        if current is not None:current['status']='failed';boxes['failed']+=1
    finally:
        if M is not None and real_acb is not None:M.acb=real_acb
        if kernel is not None and original_hint is not None:kernel._Hint=original_hint
        if trace is not None:trace.close()
    assert boxes['started']<=MAX_BOXES
    assert boxes['started']==sum(boxes[k] for k in ('completed','failed','interrupted'))
    for kind in counts:
        assert counts[kind]['started']==sum(counts[kind][k] for k in ('completed','failed','interrupted'))
        assert counts[kind]['started']<=({'AS':1500,'Psi':1500,'Hint':12000,'acb_integral':12000}[kind])
    for kind,cap in (('pieces',6000),('rows',96000)):
        assert proof_counts[kind]['started']<=cap
        assert proof_counts[kind]['started']==sum(proof_counts[kind][k] for k in ('completed','failed','interrupted'))
    if cover is not None:
        cover['accepted_tree_leaves']=len(cover['leaves'])
        cover['accepted_proof_regions']=sum(leaf.get('proof_region_count',1) for leaf in cover['leaves'])
        result['accepted_proof_regions']=cover['accepted_proof_regions']
    result.update(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        accepted_leaves=len(cover['leaves']) if cover is not None else 0,
        pending_count=len(cover['pending']) if cover is not None else None,
        unresolved_count=len(cover['unresolved']) if cover is not None else None,
        complete=cover['complete'] if cover is not None else False)
    if trace_path.exists():result['calls_sha256']=sha256(trace_path.read_bytes()).hexdigest()
    checkpoint();print(json.dumps({'status':result['status'],'new_boxes':boxes,'counts':counts,
        'accepted_leaves':result['accepted_leaves'],'elapsed_seconds':result['elapsed_seconds']}),flush=True)


if __name__=='__main__':main()
