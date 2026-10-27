"""Future exact172 bands + same-source closed corners. No calls at import."""
from fractions import Fraction as F
from itertools import product
from util173 import sha,digest,read
class IdentityFault(Exception):
 def __init__(self,reason,details):super().__init__(reason);self.details={'reason':reason,**details}
def rational(x):assert isinstance(x,str);return F(x)
def samepoint(a,b):return len(a)==len(b)==5 and all(rational(x)==rational(y) for x,y in zip(a,b))
def donor_endpoints(source,times):
 md=source['source171_metadata']['new_mass_donor'];binding=source['donor_binding'];model=binding['model_binding'];assert model['bodytable_sha256']=='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9';assert sha(model['registration_path'])==model['registration_sha256'] and sha(model['source_manifest_path'])==model['source_manifest_sha256'];reg=read(model['registration_path']);assert reg['model']==model['model'] and reg['table_sha256']==model['bodytable_sha256'];assert sha(md['packet_path'])==md['packet_sha256'];packet=read(md['packet_path']);assert packet['key']==binding['packet_key'] and samepoint(packet['point'],source['old155']['point']) and samepoint(packet['point'],md['point']);ends=packet['endpoints'];assert len(ends)==73 and [e['t'] for e in ends]==times
 if binding['kind']=='legacy092_compact':
  assert packet['original_raw_sha256']==binding['original_raw_sha256'] and sha(binding['original_raw_path'])==binding['original_raw_sha256'];raw=read(binding['original_raw_path']);p=raw['point'] if 'point' in raw else raw['anchor'];assert samepoint(p,packet['point']) and raw['endpoints']==ends
 else:
  assert binding['kind']=='full_adaptive_evidence' and packet['fingerprint']==model['expected_packet_fingerprint'] and reg['time_endpoints']==times
 for e in ends:
  lo,hi=rational(e['H']['lower']),rational(e['H']['upper']);assert 0<lo<=hi
 return ends,{'packet_path':md['packet_path'],'packet_sha256':md['packet_sha256'],'endpoint_array_sha256':digest(ends),'samepoint_checked':True,'all73_timekeys_checked':True,'producer_model_binding':model,'binding_kind':binding['kind'],'legacy_original_raw_sha256':binding.get('original_raw_sha256'),'no_fixedcutoff_scaling_or_donor_gradient':True}
def band_source_certificate(source,jet,times,old_public_log_upper):
 cid=source['canonical_id'];old=source['old155'];ref=source['original_g151_jet_ref'];assert old['epsilon']==ref['epsilon']==jet['epsilon']=='1/64';assert samepoint(old['point'],ref['point']) and samepoint(old['point'],jet['point']);assert jet['model']==ref['model'] and jet['raw_record_sha256']==ref['raw_record_sha256']==old['raw_record_sha256'];assert jet['jet_id']==ref['jet_id'];gL,gU=map(rational,jet['values'][0]);assert 0<gL<=gU;assert rational(old['Delta115_lower'])+old_public_log_upper==rational(old['Delta110_lower']);ends,identity=donor_endpoints(source,times);bytime={e['t']:e for e in ends};dL,dU=map(rational,(bytime['1/64']['H']['lower'],bytime['1/64']['H']['upper']))
 if max(gL,dL)>min(gU,dU):raise IdentityFault('same_cutoff_H64_intervals_disjoint_stop_no_clip',{'canonical_id':cid,'original_g64_interval':[str(gL),str(gU)],'donor_H64_interval':[str(dL),str(dU)],'identity':identity})
 L64=max(gL,dL);L32=rational(bytime['1/32']['H']['lower'])
 if not 0<=L32<gU:raise IdentityFault('invalid_real_H32_lower_vs_original_H64_upper',{'canonical_id':cid,'L32':str(L32),'original_g64_upper':str(gU),'identity':identity})
 used=[e for e in ends if rational(e['t'])>=F(1,64)];assert len(used)==49 and used[0]['t']=='1/64' and used[-1]['t']=='1';K=F(0);terms=[]
 for i,e in enumerate(used):
  a=rational(e['t']);b=rational(used[i+1]['t']) if i+1<len(used) else None;HUa=rational(e['H']['upper']);HLb=rational(used[i+1]['H']['lower']) if b is not None else F(0);mass=HUa-HLb;assert mass>=0
  weight=F(2,3)*(a**-6-b**-6)/(a**-4-b**-4) if b is not None else F(2,3)/a**2;assert weight>0;term=mass*weight;K+=term;terms.append({'a':str(a),'b':str(b) if b is not None else 'infinity','H_upper_a':str(HUa),'H_lower_b':str(HLb),'mass_upper':str(mass),'weight':str(weight),'weighted_upper':str(term)})
 Qband=K/(16*L64);oldq=rational(old['oldqualified155q']);uniform=F(512,3);split=F(256)*(F(7,10)-F(8,15)*L32/gU);qnew=min(oldq,Qband,uniform,split);assert 0<qnew<=oldq;assert len(terms)==49 and terms[-1]['b']=='infinity'
 return {'canonical_id':cid,'epsilon':'1/64','point':old['point'],'target110':'11/100000','target115':'23/200000','rule':'172_endpoint_bands_with153_real_tailmass_min','oldqualified155q':str(oldq),'q_new':str(qnew),'q_strictly_decreased':qnew<oldq,'Q_band':str(Qband),'Q_153_uniform':str(uniform),'Q_153_realH32_split':str(split),'original_g64_interval':[str(gL),str(gU)],'donor_H64_interval':[str(dL),str(dU)],'same_cutoff_exact_nonempty_overlap':True,'original_g64_U_unchanged':str(gU),'L64_denominator':str(L64),'denominator_uses_lower_not_upper':True,'real_H32_lower':str(L32),'K_upper':str(K),'bands':terms,'band_count':49,'infinite_tail_included':True,'old_Delta_gamma_unchanged':True,'source155_body_sha256':old['source155_body_sha256'],'original_g_raw_record_sha256':old['raw_record_sha256'],'original_jet_projection_sha256':jet['minimal_projection_sha256'],'identity':identity,'new_H_log_gradient_mean_calls':0}
def compile_curve(source,mass_certificate,mass_row_hash,C):
 o=source['old155'];m=mass_certificate;assert m['canonical_id']==source['canonical_id'] and m['old_Delta_gamma_unchanged'] and m['source155_body_sha256']==o['source155_body_sha256'];q=rational(m['q_new']);assert 0<q<=rational(o['oldqualified155q']);assert rational(o['Delta115_lower'])+C==rational(o['Delta110_lower']);return {'source':source,'mass':m,'mass_row_hash':mass_row_hash,'point':list(map(rational,o['point'])),'delta115':rational(o['Delta115_lower']),'gamma':[tuple(map(rational,g)) for g in o['loggradient_directed']],'q':q,'C':C}
def closedbox_certificate(curve,receiver):
 p=curve['point'];delta=curve['delta115'];gamma=curve['gamma'];q=curve['q'];lo=list(map(rational,receiver['box_lower']));hi=list(map(rational,receiver['box_upper']));assert p[0]==1 and len(gamma)==4 and all(a<=b for a,b in gamma);values=[]
 for bits in product((0,1),repeat=4):
  v=[(hi[j] if bits[j] else lo[j])-p[j+1] for j in range(4)];values.append((delta+sum(min(a*x,b*x) for x,(a,b) in zip(v,gamma))-q*sum(x*x for x in v),bits))
 margin,bits=min(values);source=curve['source'];return {'receiver_index':receiver['receiver_index'],'leaf_index':receiver['leaf_index'],'node':receiver['node'],'canonical_id':source['canonical_id'],'epsilon':'1/64','point':source['old155']['point'],'q_new':str(q),'rule':'172_endpoint_bands_with153_real_tailmass_min','source_mass_row_hash':curve['mass_row_hash'],'margin115':str(margin),'margin110':str(margin+curve['C']),'worst_corner_bits':list(bits),'passed115':margin>0,'single_anchor_all16':True,'checked_points':16,'box_external_anchor_allowed':True,'negative_Delta_allowed':True,'new_H_log_gradient_q_calls_in_corner':0}
def initial_best(receiver):
 import json
 b=json.loads(json.dumps(receiver['base_best158']));assert b['first_pass115'] is None and rational(b['best_margin110'])>0 and rational(b['best_margin115'])<=0;b['new173_tested_pairs']=0;b['new173_status']='base158_preserved';return b
def update_best(best,cert):
 assert best['first_pass115'] is None;best['new173_tested_pairs']+=1
 if rational(cert['margin110'])>rational(best['best_margin110']):best.update({'best_margin110':cert['margin110'],'best_margin115':cert['margin115'],'best_source_kind':'new173_band172_single_anchor','best_canonical_id':cert['canonical_id'],'best_certificate_ref':cert})
 if cert['passed115']:best['first_pass115']=cert;best['new173_status']='first115PASS_stop_receiver'
 else:best['new173_status']='best_base158_and_all_tested_newwholebox_including_failed'
 assert rational(best['best_margin110'])>0
