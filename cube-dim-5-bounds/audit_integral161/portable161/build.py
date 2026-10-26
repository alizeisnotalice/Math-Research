"""Pure byte archive: union the frozen 95-file index and 178-member archive."""
from pathlib import Path
from collections import Counter
import gzip, hashlib, io, json, tarfile

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent

def sha(data):
    return hashlib.sha256(data).hexdigest()

def put(name, obj):
    path = HERE/name
    with path.open('x') as out:
        json.dump(obj,out,indent=2,sort_keys=True);out.write('\n')

def main():
    first = SOURCE/'manifest_v2.json'
    second = SOURCE/'execution002_minimal.tar.gz'
    assert sha(first.read_bytes()) == '39781754bba8dcdb87fa390e04098176d611244ecf5acc4e52993e337b8a229c'
    assert sha(second.read_bytes()) == '23a0ec53b181c1b381d596ae43381d3bdbf7d0be838fe92b8a770b2bdeefd8ec'
    paths, blobs = {}, {}
    def add(path, data, origin):
        h = sha(data)
        if h in blobs: assert blobs[h] == data
        else: blobs[h] = data
        if path in paths:
            assert paths[path]['sha256'] == h and paths[path]['bytes'] == len(data)
            paths[path]['origins'].append(origin)
        else:
            paths[path] = {'sha256':h,'bytes':len(data),'origins':[origin]}
    files = json.loads(first.read_text())['files']; assert len(files) == 95
    for path, h in files.items():
        data = (SOURCE/path).read_bytes(); assert sha(data) == h
        add('audit_integral161/'+path,data,'manifest_v2.json')
    with tarfile.open(second,'r:gz') as archive:
        members = archive.getmembers(); assert len(members) == 178
        for member in members:
            assert member.isfile() and member.name.startswith('audit_integral161/')
            data = archive.extractfile(member).read()
            assert data == (SOURCE/member.name.removeprefix('audit_integral161/')).read_bytes()
            add(member.name,data,'execution002_minimal.tar.gz')
    assert len(paths) == 210 and len(blobs) == 165
    put('paths.json',{'format':'portable161-content-addressed-v1','paths':paths,
                      'blobs':{h:len(data) for h,data in blobs.items()}})
    failure = json.loads((SOURCE/'run001_failure_manifest.json').read_text())['files']
    omitted = {'audit_integral161/run001/'+name:row for name,row in failure.items()
               if 'audit_integral161/run001/'+name not in paths}
    put('omitted_failure_files.json',{'files':omitted,
        'reason':'Redundant progress checkpoints and START control marker are not mathematical records.',
        'failure_control_flow_replay_claim':False})
    archive = HERE/'blobs.tar.gz'
    with archive.open('xb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as compressed:
            with tarfile.open(fileobj=compressed,mode='w|',format=tarfile.USTAR_FORMAT) as tar:
                for h,data in sorted(blobs.items()):
                    entry = tarfile.TarInfo('blobs/'+h)
                    entry.size = len(data);entry.mode = 0o644;entry.mtime = 0
                    tar.addfile(entry,io.BytesIO(data))
    assert archive.stat().st_size < 100_000_000
    counts = Counter(row['sha256'] for row in paths.values())
    metadata = ['paths.json','omitted_failure_files.json','build.py','portable.py','README.md']
    put('manifest.json',{'task':'161_pure_archival_appendix','no_numerical_recomputation':True,
        'input_indexes':{'manifest_v2.json':sha(first.read_bytes()),'execution002_minimal.tar.gz':sha(second.read_bytes())},
        'union_relative_paths':len(paths),'unique_content_blobs':len(blobs),
        'mapped_original_bytes':sum(row['bytes'] for row in paths.values()),
        'unique_original_bytes':sum(map(len,blobs.values())),
        'duplicate_mappings_removed':sum(n-1 for n in counts.values()),
        'blob_archive':{'path':'blobs.tar.gz','sha256':sha(archive.read_bytes()),'bytes':archive.stat().st_size},
        'metadata':{name:{'sha256':sha((HERE/name).read_bytes()),'bytes':(HERE/name).stat().st_size} for name in metadata},
        'omitted_failure_files':len(omitted),'failure_control_flow_replay_claim':False,
        'external_raw_H_cache_and_source_packets_not_duplicated':True,
        'all_retained_file_contents_and_relative_paths_preserved':True})
    print('ARCHIVE ONLY:',len(paths),'paths;',len(blobs),'blobs;',archive.stat().st_size,'archive bytes')

if __name__ == '__main__':
    main()
