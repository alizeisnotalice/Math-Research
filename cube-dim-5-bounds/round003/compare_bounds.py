"""Compare locked round-003 C against complete rational endpoint certificates."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[1]
lockpath=root/'lower_round003/lock_c_a2_equals_a3.json'
ledgerpath=root/'results/bounds_ledger.json'
oldpath=root/'lower/rectangular_certificate.json'
C=json.loads(lockpath.read_text());ledger=json.loads(ledgerpath.read_text());old=json.loads(oldpath.read_text())
lo=F(C['exact_interval'][0])
thresholds={'checked_public_Aldaz_A5_upper':F(ledger['public_best_verified']['lower']['interval'][1]),'opening_internal_R4_upper':F(ledger['opening_thresholds']['lower']['exact_interval'][1]),'previous_certified_round001_upper':max(F(r['ratio_interval'][1]) for r in old['candidates']),'new_simple_rational':F(19943,10000)}
diffs={name:str(lo-value) for name,value in thresholds.items()}
assert all(F(value)>0 for value in diffs.values())
out={'status':'PASS_complete_exact_endpoint_comparisons','differences':diffs,'qualification':'Strict comparisons of full endpoints; display decimals unused. No new upper bound.','input_hashes':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (lockpath,ledgerpath,oldpath)},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: C above checked public, opening internal, round001, and 19943/10000')
