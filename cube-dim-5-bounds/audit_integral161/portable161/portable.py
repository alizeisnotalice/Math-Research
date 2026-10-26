"""Verify bytes or restore to a fresh directory. Never execute archived files."""
from pathlib import Path, PurePosixPath
from collections import defaultdict
import argparse, hashlib, json, re, tarfile

HERE = Path(__file__).resolve().parent
HEX = re.compile(r'[0-9a-f]{64}')

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as inp:
        for block in iter(lambda:inp.read(1<<20),b''): h.update(block)
    return h.hexdigest()

def inventory():
    manifest = json.loads((HERE/'manifest.json').read_text())
    for name,row in manifest['metadata'].items():
        assert name in {'paths.json','omitted_failure_files.json','build.py','portable.py','README.md'}
        assert digest(HERE/name) == row['sha256'] and (HERE/name).stat().st_size == row['bytes']
    arc = manifest['blob_archive']; assert arc['path'] == 'blobs.tar.gz'
    assert digest(HERE/arc['path']) == arc['sha256']
    assert (HERE/arc['path']).stat().st_size == arc['bytes']
    index = json.loads((HERE/'paths.json').read_text()); paths = index['paths']; blobs = index['blobs']
    assert len(paths) == manifest['union_relative_paths'] and len(blobs) == manifest['unique_content_blobs']
    for name,row in paths.items():
        parts = PurePosixPath(name).parts
        assert name == PurePosixPath(name).as_posix() and not name.startswith('/')
        assert parts and parts[0] == 'audit_integral161' and all(p not in {'','.','..'} for p in parts)
        assert '\\' not in name and '\x00' not in name
        assert HEX.fullmatch(row['sha256']) and row['bytes'] == blobs[row['sha256']]
    assert set(blobs) == {row['sha256'] for row in paths.values()}
    assert all(HEX.fullmatch(h) and isinstance(n,int) and n >= 0 for h,n in blobs.items())
    return manifest,paths,blobs

def verified_members(blobs):
    seen = set()
    with tarfile.open(HERE/'blobs.tar.gz','r:gz') as archive:
        for member in archive:
            assert member.isfile() and member.name.startswith('blobs/')
            h = member.name.removeprefix('blobs/')
            assert HEX.fullmatch(h) and h in blobs and h not in seen
            assert member.size == blobs[h]
            data = archive.extractfile(member).read()
            assert len(data) == member.size and hashlib.sha256(data).hexdigest() == h
            seen.add(h);yield h,data
    assert seen == set(blobs)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['verify','restore'])
    parser.add_argument('--root',type=Path)
    a = parser.parse_args(); manifest,paths,blobs = inventory()
    # Validate every byte before creating any restored output.
    for _ in verified_members(blobs): pass
    if a.mode == 'restore':
        assert a.root is not None, '--root is required'
        root = a.root.expanduser().absolute()
        assert not root.exists() and not root.is_symlink(), 'Use a fresh destination; overwrites are refused'
        root.mkdir(parents=False,mode=0o700)
        byhash = defaultdict(list)
        for name,row in paths.items(): byhash[row['sha256']].append(name)
        for h,data in verified_members(blobs):
            for name in byhash[h]:
                target = root.joinpath(*PurePosixPath(name).parts)
                target.parent.mkdir(parents=True,exist_ok=True)
                with target.open('xb') as out: out.write(data)
        for name,row in paths.items():
            target = root.joinpath(*PurePosixPath(name).parts)
            assert digest(target) == row['sha256'] and target.stat().st_size == row['bytes']
    print('PASS:',a.mode,len(paths),'original paths;',len(blobs),'verified unique blobs; no archived code executed')

if __name__ == '__main__':
    main()
