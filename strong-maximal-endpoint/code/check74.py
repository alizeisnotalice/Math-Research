import json,sys,re,math,hashlib
from pathlib import Path
from decimal import Decimal
p=Path(sys.argv[1]);out=Path(sys.argv[2])
def num(s):return Decimal(re.sub(r'`[^*]*','',s).replace('*^','E'))
py=json.load(open(p/'python/profile74.json'));wl=json.load(open(p/'wolfram/verify74.json'));ps=json.load(open(p/'python/support74.json'));ws=json.load(open(p/'wolfram/support74.json'))
assert hashlib.sha256((p/'python/profile74.json').read_bytes()).hexdigest()==wl['candidate_SHA256']
assert hashlib.sha256((p/'python/support74.json').read_bytes()).hexdigest()==ws['candidate_SHA256']
base={(r['N'],r['D_type']):r for r in py['rows'] if r['refinement']==1};refine=[abs(r['reference_mean']-base[r['N'],r['D_type']]['reference_mean']) for r in py['rows'] if r['refinement']!=1];reld=[abs(r['deletion_reference_mean_loss']-base[r['N'],r['D_type']]['deletion_reference_mean_loss'])/r['deletion_reference_mean_loss'] for r in py['rows'] if r['refinement']!=1]
assert max(refine)<1e-6 and max(reld)<1e-4
assert all(abs(r['mean_balance_residual'])<1e-10 and r['actual_mean']>1 for r in py['rows'])
low={};prec=[];cross=[];nativeagg={};bd={(r['N'],r['D_type']):r for r in py['rows'] if r['refinement']==.25}
for r in wl['rows']:
 assert r['numeric']
 if r['precision']==40:low[r['N'],r['D_type'],r['S']]=r
for r in wl['rows']:
 if r['precision']!=80:continue
 key=(r['N'],r['D_type'],r['S']);lo=low[key];j=bd[r['N'],r['D_type']];node=min(j['quadrature_nodes'],key=lambda x:abs(x['S']-r['S']));assert abs(node['S']-r['S'])<1e-9
 for fld in ['actual_removed','reference_removed','reference_log_density','h']:
  a=num(r[fld]);b=num(lo[fld]);prec.append(abs(a-b))
  if fld in node:cross.append(abs(float(a)-node[fld])/max(abs(node[fld]),1e-100))
 d=j['delta'];nativeagg.setdefault(r['D_type'],0.);nativeagg[r['D_type']]+=r['weight']*math.sqrt(j['N']*1.2784645427610738)*math.exp(float(num(r['reference_log_density'])))*float(num(r['h']))*(float(num(r['reference_removed']))-d)/(1-d)
aggerr=max(abs(v-bd[1024,k]['deletion_reference_mean_loss']) for k,v in nativeagg.items())
assert max(prec)<Decimal('1e-20') and max(cross)<1e-9 and aggerr<1e-10
sb={};spp=[];spcross=[];sd={}
for j in ps['rows']:
 for n in j['rows']:sd[j['N'],j['D_type'],n['radial_node']]=(j,n)
for r in ws['rows']:
 if r['precision']==60:sb[r['N'],r['D_type'],r['radial_node']]=r
for r in ws['rows']:
 if r['precision']!=120:continue
 key=(r['N'],r['D_type'],r['radial_node']);lo=sb[key];j,n=sd[key]
 for fld in ['log_H_ratio','log_density_ratio','log_volume_ratio']:
  a=num(r[fld]);b=num(lo[fld]);spp.append(abs(a-b));v=j[fld] if fld in j else n[fld];spcross.append(abs(float(a)-v))
assert max(spp)<Decimal('1e-20') and max(spcross)<2e-9
modes=json.load(open(p/'python/modes_summary74.json'));assert modes['all_certified'] and modes['count']==1728 and modes['candidate_SHA256']==wl['candidate_SHA256']
res={'exact_binary_mode_certificates':modes['count'],'maximum_mode_logarithm_series_terms':modes['maximum_terms'],'pass':True,'profile_Python_rows':len(py['rows']),'profile_refinement_pairs':len(refine),'profile_refinement_max_mean_absolute':max(refine),'profile_refinement_max_deletion_relative':max(reld),'Wolfram_profile_node_checks':len(wl['rows']),'profile_precision_field_pairs':len(prec),'profile_precision_max_abs':str(max(prec)),'posterior_node_cross_fields':len(cross),'posterior_node_cross_max_relative':max(cross),'native_reconstructed_deletion_max_abs':aggerr,'support_Python_rows':len(ps['rows']),'support_Wolfram_rows':len(ws['rows']),'support_precision_field_pairs':len(spp),'support_precision_max_abs':str(max(spp)),'support_cross_fields':len(spcross),'support_cross_max_abs':max(spcross),'scope':'Numerical consistency, not an interval certificate; nested-rectangle identities and asymptotic proof are separate.'}
out.write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
