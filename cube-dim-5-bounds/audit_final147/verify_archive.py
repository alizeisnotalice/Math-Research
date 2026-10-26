#!/usr/bin/env python3
"""Verify exact byte reconstruction of compact audit evidence; no arithmetic."""
import gzip,hashlib,json
from pathlib import Path

def verify(path):
    data=json.loads(gzip.decompress(Path(path).read_bytes()))
    assert data['format']=='audit147-json-shared-locks-v1'
    count=0
    for batch in data['batches'].values():
        locks=batch['input_locks']
        for name,item in batch['files'].items():
            value=item['json']
            assert 'input_locks' not in value;value={**value,'input_locks':locks}
            raw=(json.dumps(value,sort_keys=True,indent=2)+'\n').encode()
            assert hashlib.sha256(raw).hexdigest()==item['original_sha256'],name
            assert len(raw)==item['original_bytes'],name
            count+=1
    return count

if __name__=='__main__':
    path=Path(__file__).resolve().parent/'evidence/numeric_batches.json.gz'
    print(json.dumps(dict(exact_original_json_files=verify(path),new_H=0,new_log=0)))
