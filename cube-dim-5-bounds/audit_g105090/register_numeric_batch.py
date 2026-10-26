#!/usr/bin/env python3
"""Register a frozen new-only catalogue; preserve every earlier success set."""
import argparse,json
from pathlib import Path
from interfaces import REG,HERE,common as C

def main():
    p=argparse.ArgumentParser()
    for name in ['batch','points','point_sha','overlay_sha','journal_sha','registration_sha']:p.add_argument('--'+name.replace('_','-'),required=True,type=int if name in ['batch','points'] else str)
    a=p.parse_args();root=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089')/('batch%04d'%a.batch);mf=root/'point_manifest.json'
    C.require(C.file_sha(mf)==a.point_sha,'frozen point manifest mismatch');rows=json.loads(mf.read_text())['points'];C.require(len(rows)==a.points,'finite point count mismatch')
    prior=set()
    for path in HERE.glob('numeric_batch*/input_registration.json'):
        reg=json.loads(path.read_text())
        if reg['batch']>=a.batch:continue
        C.require(C.file_sha(reg['point_manifest'])==reg['point_manifest_sha256'],'prior source catalogue changed')
        prior.update(','.join(r['point']) for r in json.loads(Path(reg['point_manifest']).read_text())['points'])
    C.require(len({tuple(r['point']) for r in rows})==len(rows) and all(','.join(r['point']) not in prior for r in rows),'old successful/new anchor repeats')
    reg={'batch':a.batch,'new_points':a.points,'cache_hits_not_recomputed':0,'block_size':64,'workers':12,'target':REG['target'],'epsilon':REG['epsilon'],'source086_sha256':REG['source086_sha256'],'source_cache_count_not_recomputed':180040,
         'point_manifest':str(mf),'point_manifest_sha256':a.point_sha,'overlay':str(root/'state.json'),'overlay_sha256':a.overlay_sha,'journal':str(root/'journal.jsonl'),'journal_sha256':a.journal_sha,'production_registration':str(root/'registration.json'),'production_registration_sha256':a.registration_sha,
         'prior_new_success_points_not_recomputed':len(prior),'point_sets_disjoint_preflight':True,'wall_deadline':None,'checkpoint_policy':'fixed64point blocks;coordinator fsync+perpointACK;never repeat successful prefix','global_claim':False}
    output=HERE/('numeric_batch%04d'%a.batch);output.mkdir(exist_ok=True);path=output/'input_registration.json'
    if path.exists():C.require(json.loads(path.read_text())==reg,'existing finite registration differs')
    else:path.write_text(json.dumps(reg,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'batch':a.batch,'registered_new_points':a.points,'prior_distinct_new_points':len(prior),'input_registration_sha256':C.file_sha(path)}))

if __name__=='__main__':main()
