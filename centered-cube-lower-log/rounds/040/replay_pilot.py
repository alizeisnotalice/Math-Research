"""Replay the single frozen four-case protocol only in a new directory."""
from pathlib import Path
import argparse,os,shutil,subprocess,sys
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True);args=parser.parse_args()
root=Path(args.output).resolve()
if root.exists():raise RuntimeError('Output directory must not exist')
root.mkdir(parents=True)
for folder,names in [('round040_source',['source_tables.py','build_cache.py']),
                     ('round040_stream',['probe.py','helpers.py'])]:
    (root/folder).mkdir()
    for name in names:shutil.copy2(HERE/folder/name,root/folder/name)
subprocess.run([sys.executable,str(root/'round040_stream/probe.py')],cwd=root,check=True,
               env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
