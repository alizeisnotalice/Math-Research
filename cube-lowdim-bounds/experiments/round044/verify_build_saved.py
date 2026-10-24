"""Independent saved-value T audit, 90 seconds, standard library only.

No model/Arb/point/integral/LP evaluation. Restriction matrices are independently
obtained by Bernstein -> power -> affine substitution -> Bernstein, once per
half interval; no de Casteljau or builder import. Legacy rows are hash-pinned.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb
import json
from pathlib import Path
from time import monotonic

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'build_saved_verification.json'
SOURCE_SHA='a03f2f13909d95b140505c2c870cbfe901c2755e93d3c1666de715a603501f8f'
COVER_SHA='805a735aa4d2edddd0782a018e4fdb051259ac665ed8ac4a901896c663f29c62'
PREVIOUS_VERIFICATION_SHA='e6523c78a0c77c37d40e2e4aac7ff3fb3a6f3a0022940e911ca6ae0665d6e893'
BASE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
DIRECTIONS_SHA='98c0e2fcc7f5d448aef33acc24f2b7fe22977f7a2361b20dc7c2ade926a0fd6e'
OLD_ROWS_SHA='2affd7976754422037d2c1b396523e227b028bb82e67fffa8fb39885e034db84'
COLUMNS=(2,3,4,5,198,199,201,202,222,244,267,268)
NAMES=('coarse_0_0','coarse_1_0','coarse_1_1','coarse_2_0','fine_8_1','fine_8_0',
       'fine_9_1','fine_9_0','fine_5_3','fine_6_4','fine_4_4','fine_5_5')


def restriction_matrix(a,b):
    """Expand each old cubic Bernstein basis in power, substitute, then return to Bernstein."""
    out=[]
    for new_i in range(4):
        row=[]
        for old_i in range(4):
            coeff=Q(0)
            for power in range(old_i,4):
                old_power=Q(comb(3,old_i)*comb(3-old_i,power-old_i)*(-1)**(power-old_i))
                for new_power in range(min(power,new_i)+1):
                    coeff+=old_power*comb(power,new_power)*a**(power-new_power)*(b-a)**new_power*Q(comb(new_i,new_power),comb(3,new_power))
            row.append(coeff)
        assert all(v>=0 for v in row) and sum(row,Q(0))==1
        out.append(row)
    return out


def transform(B,Wu,Wv):
    intermediate=[[sum((w*B[k][j] for k,w in enumerate(row) if w),Q(0)) for j in range(4)] for row in Wu]
    return [[sum((w*intermediate[i][k] for k,w in enumerate(row) if w),Q(0)) for row in Wv] for i in range(4)]


def main():
    if OUT.exists():raise SystemExit('Existing audit: no repeat or overwrite.')
    start=monotonic()
    state={'status':'started','complete':False,'seconds_limit':90,'new_models_checked':0,
        'region_records_checked':0,'region_coefficients_checked':0,'new_rows_checked':0,'directions_checked':0,
        'new_Arb_calls':0,'new_models_evaluated':0,'new_points':0,'LP_calls':0}
    def save():
        state['elapsed_seconds']=monotonic()-start
        tmp=OUT.with_suffix('.json.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(OUT)
    def guard():
        if monotonic()-start>=90:raise TimeoutError('90 second saved-value verification limit; no extension.')
    def read(path):
        guard();value=json.loads(path.read_text());guard();return value
    def digest(path):
        h=sha256()
        with path.open('rb') as stream:
            while True:
                guard();block=stream.read(1<<20)
                if not block:break
                h.update(block)
        return h.hexdigest()
    save()
    try:
        run,result=read(ROOT/'build.run.json'),read(ROOT/'build_result.json')
        assert run['status']=='terminal' and run['exit_code']==0 and 0<=run['elapsed_seconds']<300
        assert run['source_sha256_before']==run['source_sha256_after']==result['source_sha256']==digest(ROOT/'build_directions.py')==SOURCE_SHA
        assert run['log_sha256']==digest(ROOT/'build.log')
        assert result['complete'] and result['status']=='complete_twelve_directions_on_full_outer_cover' and result['in_flight'] is None
        assert result['limits']==dict(seconds_hard=300,Arb_bits=192,new_AS_models=7092,tree_leaves=591,proof_regions=1278,direction_region_records=15336,new_rows=20448)
        assert all(result[k]==0 for k in ('baseline_model_calls','Psi_model_calls','new_point_values','LP_calls'))
        for inp in result['inputs']:assert digest(ROOT.parent/inp['file'])==inp['sha256']
        fixed={'round01/scaled_certificate.json':BASE_SHA,'round24/leaf_cover.json':COVER_SHA,
            'round24/saved_verification.json':PREVIOUS_VERIFICATION_SHA,'round14/directions.json':DIRECTIONS_SHA,'round14/rows.json':OLD_ROWS_SHA}
        for path,wanted in fixed.items():assert digest(ROOT.parent/path)==wanted
        base=read(ROOT.parent/'round01/scaled_certificate.json');cover=read(ROOT.parent/'round24/leaf_cover.json')
        admission=read(ROOT.parent/'round24/saved_verification.json')
        assert admission['complete'] and admission['cover_sha256']==COVER_SHA and admission['accepted_rows_checked']==20448
        assert cover['complete'] and not cover['pending'] and not cover['unresolved'] and cover['in_flight'] is None
        assert len(cover['leaves'])==591 and cover['accepted_proof_regions']==1278
        c1,h=Q(base['c']),Q(base['h'])
        assert (Q(base['alpha']),Q(base['R']),h,base['N'])==(Q(47,40),Q(33,20),Q(1,32),80) and c1>0
        direction_data=read(ROOT.parent/'round14/directions.json');directions=direction_data['directions']
        assert direction_data['complete'] and tuple(d['name'] for d in directions)==NAMES and tuple(d['old_column'] for d in directions)==COLUMNS
        assert result['directions_file']=='../round14/directions.json' and result['directions_sha256']==DIRECTIONS_SHA
        assert result['legacy_rows_file']=='../round14/rows.json' and result['legacy_rows_sha256']==OLD_ROWS_SHA
        assert result['legacy_models_reused']==12204 and result['legacy_rows_count']==16272 and result['total_continuous_rows']==36720
        assert digest(ROOT/'rows_new.json')==result['new_rows_sha256']==result['rows_new_sha256']
        assert digest(ROOT/'model_manifest.json')==result['model_manifest_sha256']
        assert digest(ROOT/'regions.json')==result['regions_sha256']
        assert digest(ROOT/'source_models_manifest.jsonl')==result['source_models_manifest_sha256']
        rowsfile=read(ROOT/'rows_new.json');rows=rowsfile['rows'];manifest=read(ROOT/'model_manifest.json');regions_saved=read(ROOT/'regions.json')['regions']
        assert rowsfile['complete'] and rowsfile['baseline_sha256']==BASE_SHA and rowsfile['cover_sha256']==COVER_SHA
        assert rowsfile['cell_rhs']=='baseline_slack/2' and tuple(rowsfile['direction_names'])==NAMES and len(rows)==20448
        assert manifest['complete'] and manifest['cover_sha256']==COVER_SHA and manifest['direction_source_sha256']==DIRECTIONS_SHA
        assert len(manifest['directions'])==12
        # The already independently verified baseline cover defines the exact immutable region schedule.
        regions=[];leaf_regions=[]
        for leaf in cover['leaves']:
            indices=[]
            pieces=leaf.get('proof_pieces',[dict(piece_id='whole',local=leaf['local'],physical=leaf['physical'],rows=leaf.get('rows'))])
            assert len(pieces)==leaf['proof_region_count']
            for p in pieces:
                ri=len(regions);indices.append(ri)
                region=dict(region_index=ri,region_key=[leaf['node_id'],p['piece_id']],leaf_node_id=leaf['node_id'],
                    proof_piece_id=p['piece_id'],parent=[52,1],local=p['local'],physical=p['physical'])
                assert regions_saved[ri]==region
                for ell,source in enumerate(p['rows']):
                    guard();row=rows[ri*16+ell]
                    assert row['region_index']==ri and row['region_key']==region['region_key'] and row['leaf_node_id']==leaf['node_id']
                    assert row['proof_piece_id']==p['piece_id'] and row['parent']==[52,1]
                    assert row['bernstein_index']==source['bernstein_index']==list(divmod(ell,4))
                    assert all(row[k]==source[k] for k in ('c_column','s_column','baseline_slack'))
                    assert Q(row['baseline_slack'])==c1*Q(row['c_column'])+Q(row['s_column'])>0
                    assert len(row['t_midpoints'])==len(row['t_penalties'])==12
                    state['new_rows_checked']+=1
                regions.append(region)
            leaf_regions.append(indices)
        assert len(regions)==len(regions_saved)==1278 and state['new_rows_checked']==20448
        # Precompute fixed half restrictions through the independent POWER derivation above.
        W=[restriction_matrix(Q(0),Q(1,2)),restriction_matrix(Q(1,2),Q(1))]
        assert restriction_matrix(Q(0),Q(1))==[[Q(i==j) for j in range(4)] for i in range(4)]
        for di,(d,item) in enumerate(zip(directions,manifest['directions'])):
            guard();state['in_flight']={'direction':d['name'],'stage':'model_and_region_streams'};save()
            assert item['name']==d['name'] and item['old_column']==d['old_column']
            assert item['new_models']==591 and item['new_region_records']==1278 and item['legacy_models']==1017
            assert item['legacy_models_file']=='../round14/'+d['models_file'] and item['legacy_models_sha256']==d['models_sha256']
            specpath=ROOT/item['spec_file'];assert digest(specpath)==item['spec_sha256']
            spec=read(specpath)
            assert spec==dict(alpha='47/40',R='33/20',h='1/32',N=80,c='0',a=d['coefficients'])
            mass=Q(0)
            for key,value in spec['a'].items():
                i,j=map(int,key.split(','));v=Q(value)
                assert 0<=j<=i<=80 and (not v or h*(i+j+4)<=Q(3,4))
                mult=1 if i==j==0 else 4 if j==0 or i==j else 8
                mass+=h*h*mult*v/2
            assert mass==Q(d['half_mass'])
            mp,rp=ROOT/item['new_models_file'],ROOT/item['new_regions_file']
            assert digest(mp)==item['new_models_sha256'] and digest(rp)==item['new_regions_sha256']
            assert mp.stat().st_size==item['new_models_bytes'] and rp.stat().st_size==item['new_regions_bytes'] and specpath.stat().st_size==item['spec_bytes']
            with mp.open() as ms,rp.open() as rs:
                for li,leaf in enumerate(cover['leaves']):
                    guard();raw=ms.readline();assert raw
                    model=json.loads(raw);line_sha=sha256(raw.encode()).hexdigest()
                    assert model['direction']==d['name'] and model['leaf_node_id']==leaf['node_id']
                    assert model['parent']==leaf['parent']==[52,1] and model['local']==leaf['local'] and model['physical']==leaf['physical']
                    assert model['source_cover_sha256']==COVER_SHA and model['spec_sha256']==item['spec_sha256']
                    elo,E=Q(model['remainder_lower']),Q(model['remainder_upper']);assert 0<=elo<=E
                    assert len(model['bernstein'])==16
                    lower=[[Q(0) for _ in range(4)] for _ in range(4)];upper=[[Q(0) for _ in range(4)] for _ in range(4)]
                    for ell,bb in enumerate(model['bernstein']):
                        i,j=divmod(ell,4);assert bb['index']==[i,j]
                        lo,hi=Q(bb['lower']),Q(bb['upper']);assert lo<=hi
                        lower[i][j]=lo;upper[i][j]=hi
                    for ri in leaf_regions[li]:
                        guard();record=json.loads(rs.readline());region=regions[ri];pid=region['proof_piece_id']
                        assert record['region_index']==ri and record['region_key']==region['region_key'] and record['direction']==d['name']
                        assert record['leaf_node_id']==leaf['node_id'] and record['source_model_line_sha256']==line_sha
                        assert record['whole_leaf_E_unchanged'] is True and Q(record['remainder_upper'])==E
                        if pid=='whole':ll,uu=lower,upper
                        else:
                            iu,iv=((0,0),(1,0),(0,1),(1,1))[pid]
                            ll,uu=transform(lower,W[iu],W[iv]),transform(upper,W[iu],W[iv])
                        assert len(record['bernstein'])==16
                        for ell,bb in enumerate(record['bernstein']):
                            guard();i,j=divmod(ell,4);assert bb['index']==[i,j]
                            lo,hi=ll[i][j],uu[i][j];mid=(lo+hi)/2;rad=(hi-lo)/2;pen=rad+E
                            assert lo<=hi and pen>=0
                            assert tuple(Q(bb[k]) for k in ('lower','upper','midpoint','radius','penalty'))==(lo,hi,mid,rad,pen)
                            row=rows[ri*16+ell]
                            assert row['t_midpoints'][di]==bb['midpoint'] and row['t_penalties'][di]==bb['penalty']
                            state['region_coefficients_checked']+=1
                        state['region_records_checked']+=1
                    state['new_models_checked']+=1
                assert not ms.readline() and not rs.readline()
            state['directions_checked']+=1;save()
        counts=result['counts']
        for key,value in dict(models_started=7092,models_completed=7092,AS_started=7092,AS_returned=7092,
            region_records_started=15336,region_records_completed=15336,baseline_rows_reused=20448,
            new_rows_written=20448,directions_completed=12,specs_started=12,specs_completed=12).items():assert counts[key]==value
        assert all(counts[k]==0 for k in ('models_failed','models_interrupted','AS_failed','region_records_failed','region_records_interrupted'))
        assert state['new_models_checked']==7092 and state['region_records_checked']==15336 and state['region_coefficients_checked']==245376
        state.update(status='passed_independent_saved_T_columns',complete=True,in_flight=None,
            scope='Only new saved T enclosures/restrictions/midpoints/penalties/rows and immutable baseline region binding. No new model or primitive evaluation; legacy 16272 rows hash-pinned only.',
            restriction_algorithm='Fixed cubic matrices derived by Bernstein-to-power, affine substitution, power-to-Bernstein; no de Casteljau calls.',
            legacy_rows_hash_pinned=16272,result_sha256=digest(ROOT/'build_result.json'),rows_sha256=result['new_rows_sha256'],
            model_manifest_sha256=result['model_manifest_sha256'],verifier_sha256=digest(Path(__file__)))
        guard()
    except TimeoutError as exc:
        state.update(status='stopped_90_second_budget',complete=False,reason=str(exc))
    except Exception as exc:
        state.update(status='failed_saved_T_audit',complete=False,reason=f'{type(exc).__name__}: {exc}')
    save();print(json.dumps(state))
    if not state['complete']:raise SystemExit(1)


if __name__=='__main__':main()
