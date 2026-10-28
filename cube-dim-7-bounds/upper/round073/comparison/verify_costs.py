from fractions import Fraction as F
from pathlib import Path
import json,os,hashlib
P=Path(__file__).resolve().parents[3]; d=json.loads((P/'results/general_cover_ledger.json').read_text());o=next(v for v in d['accepted'] if v['round']==71)
g=F(o['parameters']['gamma']);s=F(o['parameters']['sigma']);rec=json.loads((P/'upper/round073/cantelli_frontier/verification.json').read_text())['key_rational_bounds'];p=F(rec['price_upper'])
checks={'original_new_price_record':p<F(257,2500),'new_fullPhi_strictly_below_Phi71_entire_domain':g-s+F(4,8192)<F(47,250),'same_domain':o['explicit_threshold']=='natural log n>=8192','full_scalar_price_slack':F(291,1000)-F(257,2500)-F(47,250)==F(1,5000),'R069_fullPhi_still_stronger_in_own_domain':F(291,1000)-F(25,256)+F(4,2**30)<F(23,100)}
assert all(checks.values());r={'status':'PASS','pid':os.getpid(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'checks':checks,'Phi73_vs_Phi71_margin':str(F(47,250)-(g-s+F(4,8192))),'R069_domain_preserved':'natural log n>=2^30','scope':'five new record/full-cost comparisons. Old B positive price lower bounds retained; no old gate/log/producer rerun.'};Path(__file__).with_name('receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
