import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from stable_profile import run
if __name__=='__main__':
    jobs=[(n,o,18) for n in (262144,1048576,4194304) for o in (256,512)]
    with ProcessPoolExecutor(6) as pool:rows=list(pool.map(run,jobs))
    Path(__file__).with_name('stable.json').write_text(json.dumps(rows,indent=2)+'\n')
    for row in rows:print(row,flush=True)
