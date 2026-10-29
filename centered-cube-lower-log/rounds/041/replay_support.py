"""Prepare an explicit isolated data view; never edits archived evidence."""
from pathlib import Path
import hashlib,json,shutil
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare(root):
    root=Path(root)
    for name in ['variance','round041_mcd','round041_audit','round041_post','round041_replay']:
        shutil.copytree(HERE/name,root/name)
    for name in ['round040_stream','round040_source']:
        shutil.copytree(HERE.parent/'040'/name,root/name)
    (root/'round038_highdim').mkdir()
    shutil.copy2(HERE.parent/'038/round038_highdim/source_tables.py',root/'round038_highdim/source_tables.py')
    for name in ['variance_probe.py','CANDIDATE_VARIANCE.md','CANDIDATE_VARIANCE_RESULTS.json']:
        shutil.copy2(root/'variance'/name,root/'round040_source'/name)
    # Actual replay results record exactly the same NPZ SHA as round40. Reuse
    # those shared bytes in this temporary view; this is not a new model run.
    for R in [4,8,16,32]:
        meta=json.loads((root/f'round041_replay/base_replay/R{R}_RESULTS.json').read_text())
        source=root/f'round040_stream/R{R}_paired.npz'
        assert meta['paired_payload_SHA256']==sha(source)
        shutil.copy2(source,root/f'round041_replay/base_replay/R{R}_paired.npz')
    return root
