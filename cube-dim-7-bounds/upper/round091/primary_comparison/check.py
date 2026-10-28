from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,os
p=Path(__file__).resolve().parents[3]
reg=json.loads(Path(__file__).with_name('REGISTRATION.json').read_text())
for f,h in reg['input_sha256'].items():assert hashlib.sha256((p/f).read_bytes()).hexdigest()==h
old89=json.loads((p/'upper/round089/primary_comparison/receipt.json').read_text())
old88=json.loads((p/'upper/round090/primary_comparison/receipt.json').read_text())
c89=Q(old89['c89_upper']);c88=Q(old88['c88_upper']);L=2**33
checks={'H89_gap':Q(11,10)-c89-Q(20000,4987*L)>Q(1989,10000),'H88_gap':Q(11,10)-c88-Q(40,11*L)>Q(1969,10000)}
assert all(checks.values())
Path(__file__).with_name('receipt.json').write_text(json.dumps({'status':'PASS','pid':os.getpid(),'checks':checks,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':reg['input_sha256'],'scope':'new complete-cost comparisons on true common domains; no old baseline reexecution'},indent=2))
print('PASS: 2 new exact branch comparisons from archived certificates')
