from pathlib import Path
import json,sys
from upper_certificate import certify
root=Path(__file__).resolve().parents[1]
row=json.loads((root/'results'/f'round{int(sys.argv[1])}_upper.json').read_text())
actual=certify(row['R'],row['epsilon'],row['r0'])
assert all(row[k]==v for k,v in actual.items())
print('PASS',actual['cost'],actual['improvement_vs_public_exact'])
