"""Optional full same-protocol reconstruction, only in a new directory."""
from pathlib import Path
import argparse,os,shutil,subprocess,sys
from replay_support import prepare
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True);args=parser.parse_args()
root=Path(args.output).resolve()
if root.exists():raise RuntimeError('Output directory must not exist')
root.mkdir(parents=True);prepare(root)
rep=root/'round041_replay'
shutil.rmtree(rep/'base_replay')
for name in ['PLAN.json','RECONSTRUCTION_RECEIPT.json']+[f'R{r}_upper_inputs.npz'for r in [4,8,16,32]]:
    (rep/name).unlink()
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
subprocess.run([sys.executable,str(rep/'reconstruct.py')],cwd=root,env=env,check=True)
for name in ['FINAL_PLAN.json','FINAL_RESULTS.json']:(root/'round041_post'/name).unlink()
subprocess.run([sys.executable,str(root/'round041_post/aggregate.py')],cwd=root,env=env,check=True)
