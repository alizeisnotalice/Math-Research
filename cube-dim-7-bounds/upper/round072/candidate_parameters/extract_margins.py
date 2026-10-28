"""Post-run receipt extraction only; does not replay the 12-row table."""
from pathlib import Path
from fractions import Fraction as F
import datetime,hashlib,json,os
root=Path(__file__).resolve().parent
r=json.loads((root/'verification.json').read_text());v=r['constants'];base=r['base'];best=r['best_registered_row']
x=F(v['x0']);N=F(4096**8,40320);gap=F(9,200)*x-F(v['mean_upper_assumed_profile']);varhi=F(v['variance_upper_assumed_H2'])
common_margins={
 'target_mass_lower':1-F(311,1000)/F(v['target_mass_lower_Taylor20']),
 'new_price_upper':1-F(v['price_upper_Taylor7'])/F(999,10000),
 'standardized_z_squared':1-F(9,25)*varhi/(gap*gap),
 'BE_error_squared_conditional':1-F(183,1000)**2/(19*F(21,500)**2),
 'normal_density_pi_bound':1-F(199,500)**2/F(7,44),
 'retention_lower_conditional':1-F(2,3)/F(v['normal_minus_BE']),
 'source_variance_ge19_conditional':1-19/F(v['variance_lower_assumed_H2']),
 'source_exponent_at_actual_N':1-F(1,225)/(F(9,182)-F(9,200)-F(9,182)/N),
}
base_all={**common_margins,**{'base_payment_'+k:F(q) for k,q in base['normalized_margins'].items()}}
best_all={**common_margins,**{'best_payment_'+k:F(q) for k,q in best['normalized_margins'].items()}}
assert all(q>0 for q in base_all.values()) and all(q>0 for q in best_all.values())
assert min(base_all.values())>F(47,100000)
assert min(best_all.values())>F(9,1000000)
hashof=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
record={'scope':'post-run extraction of new receipt; normalized margin=1-required/available; positive critical catalogue only, not minimum over unregistered constraints or a theorem','actual_pid':os.getpid(),'actual_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'code_sha256':hashof(Path(__file__)),
 'base_catalogue':{k:str(q) for k,q in base_all.items()},'base_minimum_key':min(base_all,key=base_all.get),'base_minimum_value':str(min(base_all.values())),'base_safe_margin_floor':'47/100000',
 'best_catalogue':{k:str(q) for k,q in best_all.items()},'best_minimum_key':min(best_all,key=best_all.get),'best_minimum_value':str(min(best_all.values())),'best_safe_margin_floor':'9/1000000',
 'receipt_sha256':hashof(root/'verification.json'),'no_table_replay':True}
(root/'MARGINS.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'status':'PASS: post-run normalized-margin extraction','actual_pid':os.getpid(),'base_minimum_key':record['base_minimum_key'],'best_minimum_key':record['best_minimum_key']}))
