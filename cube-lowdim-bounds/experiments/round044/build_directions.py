#!/usr/bin/env python3
"""7092 new T AS models, restricted to 15336 direction-region records.

Only parent (52,1) is newly modeled. Legacy 12204 models / 16272 rows are
immutable references. No baseline/Psi/point/LP evaluation. Importing does not run.
"""
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from posixpath import normpath
import resource
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS,PRECISION,MAX_MODELS=300.,192,7092
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
COVER_SHA='805a735aa4d2edddd0782a018e4fdb051259ac665ed8ac4a901896c663f29c62'
VERIFICATION_SHA='e6523c78a0c77c37d40e2e4aac7ff3fb3a6f3a0022940e911ca6ae0665d6e893'
DIRECTIONS_SHA='98c0e2fcc7f5d448aef33acc24f2b7fe22977f7a2361b20dc7c2ade926a0fd6e'
OLD_ROWS_SHA='2affd7976754422037d2c1b396523e227b028bb82e67fffa8fb39885e034db84'
ALGEBRA_SHA='fd16ccf0309a4ae4627049ff730b2bb233be3948013d25cd514824febcc107b2'
CORE_SHA='c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf'
CHECK_SHA='0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8'
NAMES=('coarse_0_0','coarse_1_0','coarse_1_1','coarse_2_0','fine_8_1','fine_8_0',
       'fine_9_1','fine_9_0','fine_5_3','fine_6_4','fine_4_4','fine_5_5')
COLUMNS=(2,3,4,5,198,199,201,202,222,244,267,268)


class BudgetStop(Exception):pass


def dump(path,value):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,ensure_ascii=False,allow_nan=False)+'\n');tmp.replace(path)


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);return module


def endpoint(value):
    if not value.is_finite() or not value.is_exact():raise ValueError('Expected finite exact endpoint.')
    m,e=value.man_exp();m,e=int(m),int(e)
    return Q(m*(1<<e)) if e>=0 else Q(m,1<<(-e))


def interval(value):
    if not value.is_finite():raise ValueError('Nonfinite ball.')
    lo,hi=endpoint(value.lower()),endpoint(value.upper())
    if lo>hi:raise ValueError('Reversed interval.')
    return lo,hi


def main():
    if any((ROOT/n).exists() for n in ('build_result.json','rows_new.json','models.local','regions.local','model_manifest.json')):
        raise SystemExit('Existing build evidence: no repeat/resume.')
    start=monotonic();current=None;active=[];origin_stream=None
    counts={k:0 for k in ('models_started','models_completed','models_failed','models_interrupted','AS_started','AS_returned','AS_failed',
        'region_records_started','region_records_completed','region_records_failed','region_records_interrupted',
        'baseline_rows_reused','new_rows_written','directions_completed','specs_started','specs_completed','source_files_verified')}
    result={'status':'started','complete':False,'counts':counts,'in_flight':None,
        'limits':{'seconds_hard':300,'Arb_bits':192,'new_AS_models':7092,'tree_leaves':591,'proof_regions':1278,
            'direction_region_records':15336,'new_rows':20448},
        'baseline_model_calls':0,'Psi_model_calls':0,'new_point_values':0,'LP_calls':0,
        'legacy_models_reused':12204,'legacy_rows_count':16272,'total_continuous_rows':36720,
        'legacy_rows_file':'../round14/rows.json','legacy_rows_sha256':OLD_ROWS_SHA,
        'directions_file':'../round14/directions.json','directions_sha256':DIRECTIONS_SHA,
        'new_rows_file':'rows_new.json','new_rows_sha256':None,'rows_new_sha256':None,
        'completed_direction_files':[]}
    def guard():
        if monotonic()-start>=SECONDS:raise BudgetStop('internal_300s_deadline')
    def save():result['elapsed_seconds']=monotonic()-start;dump(ROOT/'build_result.json',result)
    def digest(path):
        h=sha256()
        with path.open('rb') as stream:
            while True:
                guard();block=stream.read(1<<20)
                if not block:break
                h.update(block)
        return h.hexdigest()
    def verify(relative,wanted):
        path=ROOT.parent/relative;actual=digest(path)
        if actual!=wanted:raise ValueError('Pinned input SHA mismatch: '+relative)
        result.setdefault('inputs',[]).append({'file':relative,'sha256':actual});return path
    def line(stream,record):
        raw=json.dumps(record,ensure_ascii=False,allow_nan=False,separators=(',',':'))+'\n'
        stream.write(raw);stream.flush();return sha256(raw.encode()).hexdigest()
    save()
    manifest={'complete':False,'cover_sha256':COVER_SHA,'direction_source_sha256':DIRECTIONS_SHA,'directions':[]}
    try:
        baseline_path=verify('round01/scaled_certificate.json',BASELINE_SHA)
        cover_path=verify('round24/leaf_cover.json',COVER_SHA)
        verification_path=verify('round24/saved_verification.json',VERIFICATION_SHA)
        direction_path=verify('round14/directions.json',DIRECTIONS_SHA)
        verify('round14/rows.json',OLD_ROWS_SHA)
        algebra_path=verify('round23/quadrants_once.py',ALGEBRA_SHA)
        base=json.loads(baseline_path.read_text());c1,h=Q(base['c']),Q(base['h'])
        if (Q(base['alpha']),Q(base['R']),h,base['N'])!=(Q(47,40),Q(33,20),Q(1,32),80) or c1<=0:
            raise ValueError('Fixed baseline family changed.')
        cover=json.loads(cover_path.read_text());verification=json.loads(verification_path.read_text())
        if not cover['complete'] or cover['baseline_sha256']!=BASELINE_SHA or cover['pending'] or cover['unresolved'] or cover['in_flight']:
            raise ValueError('Incomplete/wrong saved baseline cover.')
        if verification['status']!='passed_saved_arithmetic_tree_frontier_audit' or not verification['complete'] or \
                verification['cover_sha256']!=COVER_SHA or verification['accepted_leaves']!=591 or \
                verification['accepted_proof_regions']!=1278 or verification['accepted_rows_checked']!=20448:
            raise ValueError('Independent full cover admission missing.')
        leaves=cover['leaves']
        if len(leaves)!=591 or len({v['node_id'] for v in leaves})!=591:raise ValueError('Frozen tree leaf set changed.')
        source_dir=json.loads(direction_path.read_text())
        directions=source_dir['directions']
        if not source_dir['complete'] or source_dir['baseline_sha256']!=BASELINE_SHA or \
                tuple(d['name'] for d in directions)!=NAMES or tuple(d['old_column'] for d in directions)!=COLUMNS:
            raise ValueError('Frozen twelve directions changed.')
        for d in directions:
            guard();name=d['name'];label=name.split('_')[0];mass=Q(0)
            for key,value in d['coefficients'].items():
                i,j=map(int,key.split(','));v=Q(value)
                if not (0<=j<=i<=80) or (v and h*(i+j+4)>Q(3,4)):raise ValueError('Direction support/canonical index invalid.')
                mult=1 if i==j==0 else 4 if j==0 or i==j else 8
                mass+=h*h*mult*v/2
            if mass!=Q(d['half_mass']) or Q(d['bound'])!=(Q(1,2) if label=='coarse' else Q(1,4)) or d['models']!=1017:
                raise ValueError('Direction mass/bound/legacy count mismatch.')
            old_model=ROOT.parent/'round14'/d['models_file']
            if digest(old_model)!=d['models_sha256']:raise ValueError('Legacy direction model SHA mismatch.')
        origin_stream=(ROOT/'source_models_manifest.jsonl').open('x');active.append(origin_stream);source_seen={}
        def source_file(relative,wanted):
            path=cover_path.parent/relative;key=normpath('round24/'+relative)
            if key in source_seen:
                if source_seen[key]!=wanted:raise ValueError('Conflicting saved source SHA.')
                return
            if digest(path)!=wanted:raise ValueError('Saved baseline model/supplemental SHA mismatch.')
            source_seen[key]=wanted;line(origin_stream,{'file':key,'sha256':wanted});counts['source_files_verified']+=1
        regions=[];leaf_regions=[]
        for leaf in leaves:
            guard();node=cover['nodes'][leaf['node_id']]
            if node['status']!='accepted' or node['local']!=leaf['local'] or node['physical']!=leaf['physical'] or leaf['parent']!=[52,1]:
                raise ValueError('Saved accepted membership/geometry mismatch.')
            source_file(leaf['model_file'],leaf['model_sha256'])
            for fk,sk in (('supplemental_Psi_file','supplemental_Psi_sha256'),('supplemental_signed_file','supplemental_signed_sha256'),
                    ('supplemental_piece_file','supplemental_piece_sha256'),('proof_file','proof_sha256')):
                if fk in leaf:source_file(leaf[fk],leaf[sk])
            x0,x1,y0,y1=map(Q,leaf['local']);u0,u1,v0,v1=map(Q,leaf['physical'])
            if (u0,u1,v0,v1)!=(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1)) or not (u0>0 and v0>0 and u0+v0>Q(33,20)):
                raise ValueError('Saved exterior leaf geometry invalid.')
            if 'proof_pieces' not in leaf:
                candidates=[{'piece_id':'whole','local':leaf['local'],'physical':leaf['physical'],'rows':leaf['rows']}]
                if leaf.get('proof_region_count',1)!=1:raise ValueError('Whole-region count mismatch.')
            else:
                candidates=leaf['proof_pieces']
                if leaf['proof_region_count']!=4 or len(candidates)!=4 or [p['piece_id'] for p in candidates]!=list(range(4)):
                    raise ValueError('Four-piece schema mismatch.')
            indices=[]
            for pindex,piece in enumerate(candidates):
                guard();pid=piece['piece_id']
                if pid!='whole':
                    iu,iv=((0,0),(1,0),(0,1),(1,1))[pid];xc,yc=(x0+x1)/2,(y0+y1)/2
                    expected=(*((x0,xc) if iu==0 else (xc,x1)),*((y0,yc) if iv==0 else (yc,y1)))
                    if tuple(map(Q,piece['local']))!=expected:raise ValueError('Piece is not the frozen half/half restriction.')
                px0,px1,py0,py1=map(Q,piece['local'])
                physical=(h*(52+px0),h*(52+px1),h*(1+py0),h*(1+py1))
                if tuple(map(Q,piece['physical']))!=physical or len(piece['rows'])!=16:raise ValueError('Piece physical geometry/rows mismatch.')
                for ell,row in enumerate(piece['rows']):
                    guard()
                    if row['bernstein_index']!=[ell//4,ell%4] or Q(row['baseline_slack'])!=c1*Q(row['c_column'])+Q(row['s_column']) or Q(row['baseline_slack'])<=0:
                        raise ValueError('Saved exact positive c/s row invalid.')
                    counts['baseline_rows_reused']+=1
                index=len(regions);indices.append(index)
                regions.append({'region_index':index,'region_key':[leaf['node_id'],pid],'leaf_node_id':leaf['node_id'],
                    'proof_piece_id':pid,'parent':[52,1],'local':piece['local'],'physical':piece['physical'],
                    '_rows':piece['rows']})
            leaf_regions.append(indices)
        if len(regions)!=1278 or counts['baseline_rows_reused']!=20448:raise ValueError('1278-region schedule changed.')
        origin_stream.close();active.remove(origin_stream);origin_stream=None
        result['source_models_manifest_sha256']=digest(ROOT/'source_models_manifest.jsonl')
        dump(ROOT/'regions.json',{'cover_sha256':COVER_SHA,'regions':[{k:v for k,v in region.items() if k!='_rows'} for region in regions]})
        result['regions_sha256']=digest(ROOT/'regions.json');save()
        public=ROOT/'public_src'
        if not (public/'hl_check.py').is_file():public=ROOT.parent/'upstream-v4/checker/independent'
        for p,w in ((public/'hl_core.py',CORE_SHA),(public/'hl_check.py',CHECK_SHA)):
            if digest(p)!=w:raise ValueError('Pinned independent formula source mismatch.')
        algebra=load('round25_frozen_restriction',algebra_path)
        core=load('hl_core',public/'hl_core.py');checker=load('round25_hl_check',public/'hl_check.py')
        if checker.Spec is not core.Spec:raise ValueError('Independent core module mismatch.')
        from flint import arb_mat,ctx
        import flint
        ctx.prec=PRECISION
        checker.MC3=arb_mat([[core.A(x) for x in row] for row in checker.power_to_bern_matrix(3,Q(-1),Q(1))])
        checker.MC3T=checker.MC3.transpose()
        result.update(flint_version=flint.__version__,python=sys.version,MC3_rebuilt_after_precision=True)
        for sub in ('models.local','regions.local','specs.local'):(ROOT/sub).mkdir()
        call_keys=set()
        for di,d in enumerate(directions):
            guard();name=d['name'];spec_path=ROOT/'specs.local'/(name+'.json')
            dump(spec_path,{'alpha':'47/40','R':'33/20','h':'1/32','N':80,'c':'0','a':d['coefficients']})
            spec_sha=digest(spec_path);counts['specs_started']+=1;save()
            spec=core.Spec(spec_path);ak=checker.AKChecker(spec);counts['specs_completed']+=1;guard()
            mp,rp=ROOT/'models.local'/(name+'.jsonl'),ROOT/'regions.local'/(name+'.jsonl')
            ms,rrs=mp.open('x'),rp.open('x');active.extend((ms,rrs))
            for li,leaf in enumerate(leaves):
                guard();key=(name,leaf['node_id'])
                if key in call_keys:raise ValueError('Repeated direction/tree-leaf model.')
                if counts['models_started']>=MAX_MODELS:raise BudgetStop('7092_model_cap')
                call_keys.add(key);counts['models_started']+=1;counts['AS_started']+=1
                current={'direction':name,'direction_index':di,'leaf_node_id':leaf['node_id'],'stage':'AS'}
                result['in_flight']=dict(current);save()
                p,q=leaf['parent'];x0,x1,y0,y1=map(Q,leaf['local'])
                try:
                    P,E=ak.AS_model(p,q,x0,x1,y0,y1);counts['AS_returned']+=1
                except Exception:
                    counts['AS_failed']+=1;raise
                guard();elo,eup=interval(E)
                if elo<0:raise ValueError('Negative T remainder; no clamp.')
                B=checker.MC3*P*checker.MC3T;lower=[];upper=[];bern=[]
                for i in range(4):
                    lrow=[];urow=[]
                    for j in range(4):
                        guard();lo,hi=interval(B[i,j]);lrow.append(lo);urow.append(hi)
                        bern.append({'index':[i,j],'lower':str(lo),'upper':str(hi)})
                    lower.append(lrow);upper.append(urow)
                model={'direction':name,'leaf_node_id':leaf['node_id'],'parent':leaf['parent'],'local':leaf['local'],'physical':leaf['physical'],
                    'source_cover_sha256':COVER_SHA,'spec_sha256':spec_sha,'remainder_lower':str(elo),'remainder_upper':str(eup),
                    'bernstein':bern,'power_balls':[[P[i,j].str(70) for j in range(4)] for i in range(4)]}
                model_line_sha=line(ms,model);current['stage']='exact_restriction';result['in_flight']=dict(current);save()
                if len(leaf_regions[li])==4:rl,ru=algebra.tensor_quadrants(lower),algebra.tensor_quadrants(upper)
                for ri in leaf_regions[li]:
                    guard();counts['region_records_started']+=1;region=regions[ri];pid=region['proof_piece_id']
                    if pid=='whole':ll,uu=lower,upper
                    else:
                        side=((0,0),(1,0),(0,1),(1,1))[pid];ll,uu=rl[side],ru[side]
                    record={'region_index':ri,'region_key':region['region_key'],'direction':name,'leaf_node_id':leaf['node_id'],
                        'source_model_line_sha256':model_line_sha,'remainder_upper':str(eup),'whole_leaf_E_unchanged':True,'bernstein':[]}
                    for i in range(4):
                        for j in range(4):
                            guard();lo,hi=ll[i][j],uu[i][j]
                            if lo>hi:raise ValueError('Restricted endpoint order failed.')
                            mid,rad=(lo+hi)/2,(hi-lo)/2;pen=rad+eup
                            if pen<0:raise ValueError('Negative signed T penalty.')
                            record['bernstein'].append({'index':[i,j],'lower':str(lo),'upper':str(hi),
                                'midpoint':str(mid),'radius':str(rad),'penalty':str(pen)})
                    line(rrs,record);counts['region_records_completed']+=1
                counts['models_completed']+=1;current=None;result['in_flight']=None;save();guard()
            ms.close();rrs.close();active.remove(ms);active.remove(rrs)
            item={'name':name,'old_column':d['old_column'],'spec_file':str(spec_path.relative_to(ROOT)),'spec_sha256':spec_sha,
                'new_models_file':str(mp.relative_to(ROOT)),'new_models_sha256':digest(mp),'new_models':591,
                'new_regions_file':str(rp.relative_to(ROOT)),'new_regions_sha256':digest(rp),'new_region_records':1278,
                'new_models_bytes':mp.stat().st_size,'new_regions_bytes':rp.stat().st_size,'spec_bytes':spec_path.stat().st_size,
                'legacy_models_file':'../round14/'+d['models_file'],'legacy_models_sha256':d['models_sha256'],'legacy_models':1017}
            manifest['directions'].append(item);result['completed_direction_files'].append(item);counts['directions_completed']+=1
            dump(ROOT/'model_manifest.json',manifest);save();del ak,spec;guard()
        if counts['models_completed']!=7092 or counts['region_records_completed']!=15336:raise ValueError('Incomplete model/region schedule.')
        streams=[(ROOT/d['new_regions_file']).open() for d in manifest['directions']];active.extend(streams)
        partial=ROOT/'rows_new.partial.json';row_stream=partial.open('x');active.append(row_stream)
        header={'complete':True,'baseline_sha256':BASELINE_SHA,'cover_sha256':COVER_SHA,'direction_names':list(NAMES),'cell_rhs':'baseline_slack/2'}
        row_stream.write(json.dumps(header,separators=(',',':'))[:-1]+',"rows":[\n');first=True
        for region in regions:
            guard();records=[json.loads(stream.readline()) for stream in streams]
            for di,record in enumerate(records):
                if record['region_index']!=region['region_index'] or record['region_key']!=region['region_key'] or record['direction']!=NAMES[di]:
                    raise ValueError('Direction/region stream alignment mismatch.')
            for ell,base_row in enumerate(region['_rows']):
                guard();mids=[];pens=[]
                for record in records:
                    b=record['bernstein'][ell];lo,hi,mid,rad,pen=map(Q,(b['lower'],b['upper'],b['midpoint'],b['radius'],b['penalty']))
                    if b['index']!=base_row['bernstein_index'] or lo>hi or mid!=(lo+hi)/2 or rad!=(hi-lo)/2 or pen!=rad+Q(record['remainder_upper']):
                        raise ValueError('Stored midpoint/radius/whole E arithmetic mismatch.')
                    mids.append(b['midpoint']);pens.append(b['penalty'])
                row={key:base_row[key] for key in ('bernstein_index','c_column','s_column','baseline_slack')}
                row.update(leaf_node_id=region['leaf_node_id'],region_index=region['region_index'],region_key=region['region_key'],
                    proof_piece_id=region['proof_piece_id'],parent=[52,1],t_midpoints=mids,t_penalties=pens)
                if not first:row_stream.write(',\n')
                row_stream.write(json.dumps(row,separators=(',',':')));first=False;counts['new_rows_written']+=1
        for stream in streams:
            if stream.readline():raise ValueError('Unexpected trailing direction-region record.')
            stream.close();active.remove(stream)
        if counts['new_rows_written']!=20448:raise ValueError('20448 new rows not complete.')
        row_stream.write('\n]}\n');row_stream.close();active.remove(row_stream);guard()
        partial.replace(ROOT/'rows_new.json');new_sha=digest(ROOT/'rows_new.json')
        manifest['complete']=True;dump(ROOT/'model_manifest.json',manifest)
        result.update(status='complete_twelve_directions_on_full_outer_cover',complete=True,new_rows_sha256=new_sha,rows_new_sha256=new_sha,
            model_manifest_sha256=digest(ROOT/'model_manifest.json'))
    except BudgetStop as exc:
        result.update(status='interrupted_saved_checkpoint',reason=str(exc))
        if current is not None:counts['models_interrupted']+=1
        counts['region_records_interrupted']+=counts['region_records_started']-counts['region_records_completed']
    except Exception as exc:
        result.update(status='failed_preparation_or_model',reason=f'{type(exc).__name__}: {exc}')
        if current is not None:counts['models_failed']+=1
        counts['region_records_failed']+=counts['region_records_started']-counts['region_records_completed']
    finally:
        for stream in active:
            if not stream.closed:stream.close()
    assert counts['models_started']<=7092
    assert counts['models_started']==sum(counts[k] for k in ('models_completed','models_failed','models_interrupted'))
    assert counts['AS_started']==counts['models_started'] and counts['AS_returned']+counts['AS_failed']==counts['AS_started']
    assert counts['region_records_started']<=15336
    assert counts['region_records_started']==sum(counts[k] for k in ('region_records_completed','region_records_failed','region_records_interrupted'))
    sizes={}
    for name in ('rows_new.json','rows_new.partial.json','regions.json','model_manifest.json','source_models_manifest.jsonl'):
        file=ROOT/name
        if file.is_file():sizes[name]=file.stat().st_size
    for directory in ('models.local','regions.local','specs.local'):
        path=ROOT/directory
        if path.is_dir():
            for file in path.iterdir():
                if file.is_file():sizes[str(file.relative_to(ROOT))]=file.stat().st_size
    peak_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result.update(artifact_bytes=sizes,total_artifact_bytes=sum(sizes.values()),platform=sys.platform,
        ru_maxrss=peak_rss,ru_maxrss_unit='bytes' if sys.platform=='darwin' else 'KiB on Linux; platform-native otherwise')
    result['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();save()
    print(json.dumps({'status':result['status'],'complete':result['complete'],'counts':counts,'elapsed_seconds':result['elapsed_seconds']}),flush=True)


if __name__=='__main__':main()
