#!/usr/bin/env python3
"""Materialize a repository-relative packet in a new reader directory.

Existing shared library files may replace byte-identical metadata copies in
the archive. SHA is checked before hardlink/copy. No mathematical evaluation.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('packet',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args();p=a.packet.resolve()
    done=json.loads((p/'completion.json').read_bytes());refs=json.loads((p/'repository_reuse.json').read_bytes())['files']
    a.output.mkdir(exist_ok=False)
    for n,h in done['files'].items():
        src=p/n
        if not src.exists():
            assert n in refs and refs[n]['sha256']==h
            src=(p/refs[n]['relative_from_repository_packet']).resolve()
        assert sha(src)==h;dest=a.output/n;dest.parent.mkdir(parents=True,exist_ok=True)
        try:os.link(src,dest)
        except OSError:shutil.copyfile(src,dest)
        assert sha(dest)==h
    shutil.copyfile(p/'completion.json',a.output/'completion.json')
    print(json.dumps({'status':'byte_identical_packet_materialized','mathematical_calls':0,'output':str(a.output)}))
if __name__=='__main__':main()
