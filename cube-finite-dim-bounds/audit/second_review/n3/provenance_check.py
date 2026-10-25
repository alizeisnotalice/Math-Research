"""Optional local historical stdout/raw crosscheck; does not establish exit code."""
from pathlib import Path
import argparse,json,hashlib
ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);args=ap.parse_args()
base=args.directory;checks=[]
for k in range(12):
    p=base/f'part_{k:02}.log';lines=p.read_text().splitlines()
    rows=[json.loads(x) for x in lines[:688]]
    assert len({tuple(r['patch']) for r in rows})==688
    summary=json.loads('\n'.join(lines[688:]));assert summary['part']==k
    raw=json.loads((base/f'bridge-raw-archive/cb11_integral_cheb_part_{k:02}.json').read_text())
    assert all(all(r[key]==raw['patches'][i][key] for key in r) for i,r in enumerate(rows))
    assert all(summary[key]==raw[key] for key in summary)
    checks.append({'part':k,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'progress_rows':688,'raw_fields_match':True})
launch=json.loads((base/'launch.json').read_text());assert launch['parts']==12 and len(set(launch['pids']))==12
out={'status':'PASS historical execution logs/progress/final summary equal archived raw results','logs':checks,'launch_sha256':hashlib.sha256((base/'launch.json').read_bytes()).hexdigest(),'exit_codes':'Inherited original parent execution receipt all_exit_codes=0; stdout files alone do not prove exit status.','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('historical_log_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print('12 logs / 8256 rows match raw data exactly')
