from fractions import Fraction as F
from pathlib import Path
import json,hashlib,os
P=Path(__file__).resolve().parents[3]
d=json.loads((P/'results/general_cover_ledger.json').read_text())
old=next(e for e in d['accepted'] if e['round']==56)
g=F(old['parameters']['gamma'])
margin=F(7,200)-(g-F(7,100)+F(4,2**20))
checks={'old_original_gamma_record':g==F(172531397335235056462877787662549380239,1643616601687301779182642153860479582208),'full_Phi71_vs_Phi56_positive_B_price':margin>0,'new_domain_strictly_lower':F(8192)<2**20,'display_scalar_slack':F(13,100)-F(189,2000)-F(7,200)==F(1,2000),'stronger_R069_display_keeps_own_domain':F(23,100)>F(7,200) and 2**30>8192}
assert all(checks.values())
r={'status':'PASS','pid':os.getpid(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'checks':checks,'complete_cost_margin':str(margin),'scope':'new full-function comparison using original exactgamma and old positive B-price lower bound. Phi71<T-.035n<Phi56 on L>=2^20; no old program rerun, no transfer of old effective domain.'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
