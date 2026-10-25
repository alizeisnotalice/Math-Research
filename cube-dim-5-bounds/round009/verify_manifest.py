"""Check the explicit round009 payload; optionally read a committed Git snapshot."""
import argparse, hashlib, json, subprocess
from pathlib import Path
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--revision',help='Git revision of the round009 commit, for mutable ledger/STATE snapshots')
    args=parser.parse_args();project=Path(__file__).resolve().parents[1];repo=project.parent
    def read(rel):
        if args.revision:
            return subprocess.check_output(['git','show',args.revision+':cube-dim-5-bounds/'+rel],cwd=repo)
        return (project/rel).read_bytes()
    manifest=json.loads(read('round009/manifest.json'));checked=0
    for row in manifest['white_list']:
        rel=Path(row['path']);assert not rel.is_absolute() and '..' not in rel.parts
        data=read(str(rel));assert len(data)==row['bytes'],str(rel)
        assert hashlib.sha256(data).hexdigest()==row['sha256'],str(rel);checked+=1
    print(json.dumps({'status':'PASS_explicit_round009_payload','files':checked,
                      'mode':'git_snapshot' if args.revision else 'current_filesystem',
                      'revision':args.revision,'mathematical_experiments_run':0}))
if __name__=='__main__':main()
