"""Verify explicit manifests from round009 onward; optionally a Git snapshot."""
from pathlib import Path
import argparse,hashlib,json,subprocess
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('round',type=int);ap.add_argument('--revision');a=ap.parse_args()
    P=Path(__file__).resolve().parents[1];repo=P.parent
    def read(rel):
        if a.revision:return subprocess.check_output(['git','show',a.revision+':cube-dim-5-bounds/'+rel],cwd=repo)
        return (P/rel).read_bytes()
    manifest=json.loads(read(f'round{a.round:03}/manifest.json'));rows=manifest['white_list']
    for row in rows:
        r=Path(row['path']);assert not r.is_absolute() and '..' not in r.parts
        b=read(str(r));assert len(b)==row['bytes'] and hashlib.sha256(b).hexdigest()==row['sha256'],str(r)
    print(json.dumps({'status':'PASS_explicit_round_payload','round':a.round,'files':len(rows),'revision':a.revision,'mathematical_experiments_run':0}))
if __name__=='__main__':main()
