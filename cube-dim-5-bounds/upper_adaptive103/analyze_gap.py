"""Zero-H weighted gap and fixed-current-mean optimistic family diagnostic."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
from math import isqrt,floor,ceil
import gzip,json,hashlib,time,os
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive103');N=1<<160;D=1<<256

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
def root(q,k):
 n=q.numerator*D**k//q.denominator
 if k==4:v=isqrt(isqrt(n))
 else:
  v=1<<((n.bit_length()+4)//5)
  while True:
   z=(4*v+n//v**4)//5
   if z>=v:break
   v=z
  while (v+1)**5<=n:v+=1
  while v**5>n:v-=1
 a=F(v,D);b=a if a**k==q else F(v+1,D);assert a**k<=q<=b**k;return a,b

def main():
 t=time.monotonic();reg=read(H/'gap_registration.json');assert sha(__file__)==reg['source_sha256'];state=read(S/'state.json.gz');assert sha(S/'state.json.gz')==reg['state_sha256'];nodes=state['nodes'];bounds={};inputs=[]
 for p in sorted((H.parent/'upper_range092').glob('block*_cells.json.gz')):
  for r in read(p):bounds[r['key']]={'L_i':r['rangeL_i'],'U_i':r['rangeU_i'],'M_i':r['mean_i'],'integral_lower_i':r['new_power_lower_i'],'integral_upper_i':r['new_chord_upper_i']}
  inputs.append({'path':str(p.resolve()),'sha256':sha(p)})
 for folder in (H.parent/'upper_adaptive095',H):
  p=folder/'children_certificates.json.gz';bounds.update({r['key']:r['bounds'] for r in read(p)});inputs.append({'path':str(p.resolve()),'sha256':sha(p)})
 assert set(bounds)==set(nodes);front=set(state['frontier']);records=[];levels=defaultdict(lambda:Counter());shapes=defaultdict(lambda:Counter());tot=Counter();optimisticLo={};optimisticHi={}
 for k,n in nodes.items():
  b=bounds[k];l,u,m=F(b['L_i'],N),F(b['U_i'],N),F(b['M_i'],N);assert 0<l<=m<=u;v=F(n['weight'],2**(4*n['level']));ml,mh=root(m**5,4);optimisticLo[k]=min(n['selfU'],floor(v*ml*N));optimisticHi[k]=min(n['selfU'],ceil(v*mh*N))
  if k not in front:continue
  ll,lh=root(l**5,4);ul,uh=root(u**5,4);cl,ch=(v*ll,v*lh) if l==u else (v*((u-m)*ll+(m-l)*ul)/(u-l),v*((u-m)*lh+(m-l)*uh)/(u-l));r={'key':k,'level':n['level'],'weight':n['weight'],'gap_units':n['selfU']-n['selfL'],'selfL':n['selfL'],'selfU':n['selfU'],'range_power_lower_units':floor(v*ll*N),'range_power_upper_units':ceil(v*uh*N),'chord_lower_units':floor(cl*N),'chord_upper_units':ceil(ch*N),'mean_power_lower_units':floor(v*ml*N),'mean_power_upper_units':ceil(v*mh*N),'range_to_chord_gain_lower_units':floor((v*ul-ch)*N),'range_to_chord_gain_upper_units':ceil((v*uh-cl)*N),'chord_to_mean_power_surplus_lower_units':floor((cl-v*mh)*N),'chord_to_mean_power_surplus_upper_units':ceil((ch-v*ml)*N),'mean_power_to_range_lower_spread_lower_units':floor(v*(ml-lh)*N),'mean_power_to_range_lower_spread_upper_units':ceil(v*(mh-ll)*N),'cert_chord_upper_units':n['weight']*b['integral_upper_i'],'cert_power_lower_units':n['weight']*b['integral_lower_i'],'cert_aggregate_upper_retention_units':n['weight']*b['integral_upper_i']-n['selfU'],'cert_aggregate_lower_retention_units':n['selfL']-n['weight']*b['integral_lower_i']};assert r['cert_aggregate_upper_retention_units']>=0 and r['cert_aggregate_lower_retention_units']>=0
  # Root enclosure widths can move an exact zero component by at most one grid unit.
  for name in ('range_to_chord_gain','chord_to_mean_power_surplus','mean_power_to_range_lower_spread'):assert r[name+'_upper_units']>=0
  records.append(r);fields={a:z for a,z in r.items() if a.endswith('_units') or a in ('selfL','selfU')};fields['cells']=1
  for cc in (tot,levels[str(n['level'])],shapes[str(sorted(Counter(n['indices']).values(),reverse=True))]):cc.update(fields)
 for k in sorted(nodes,key=lambda k:nodes[k]['level'],reverse=True):
  n=nodes[k]
  if n['children']:optimisticLo[k]=min(optimisticLo[k],sum(optimisticLo[c] for c in n['children']));optimisticHi[k]=min(optimisticHi[k],sum(optimisticHi[c] for c in n['children']))
 lo=sum(nodes[k]['effL'] for k in state['roots']);up=sum(nodes[k]['effU'] for k in state['roots']);aU=tot['selfU']-up;aL=lo-tot['selfL'];assert aU>=0 and aL>=0;assert up-lo==tot['gap_units']-aU-aL;assert tot['cert_chord_upper_units']-tot['cert_power_lower_units']==tot['gap_units']+tot['cert_aggregate_upper_retention_units']+tot['cert_aggregate_lower_retention_units'];oL=sum(optimisticLo[k] for k in state['roots']);oU=sum(optimisticHi[k] for k in state['roots']);assert oL<=oU<=up;g=F(21,200000);p=read(H/'registration.json');E=F(p['E076_upper']);IL=F(lo,N);IU=F(read(H/'run_receipt.json')['after']['I_upper']);TL,TU=root(32**4*(g/5)**5,4);B0L=root(F(oL,N)**4,5)[0]+E;B0U=root(F(oU,N)**4,5)[1]+E;I0L=root(B0L**5,4)[0];I0U=root(B0U**5,4)[1];rank=sorted(records,key=lambda r:(-r['gap_units'],r['level'],nodes[r['key']]['indices']));concentration={}
 for count in (1,64,256,1024,2048,8192):
  got=sum(r['gap_units'] for r in rank[:count]);concentration[str(count)]={'sum_gap_units':got,'fraction_of_physical_selfgap':str(F(got,tot['gap_units']))}
 out={'status':'pass_zero_H_diagnostic','source_guidance_sha256':reg['guidance_sha256'],'all_nodes':len(nodes),'physical_leaves':len(front),'grid_denominator':str(N),'weighted_volume':'1','global_LH_units':lo,'global_UH_units':up,'global_hierarchy_gap_units':up-lo,'physical_selfgap_units':tot['gap_units'],'ancestor_upper_retention_gain_units':aU,'ancestor_lower_retention_gain_units':aL,'exact_gap_decomposition_verified':True,'fullG_I_lower':str(IL),'I_lower_bridge':'p_t>=0 implies G>=H512 pointwise; monotone q5/4 gives full integral I>=integral H512^q>=LH. This is not a pointwise global G or g lower bound.','fullG_I_upper':str(IU),'threshold_T_g105_lower':str(TL),'threshold_T_g105_upper':str(TU),'certificate_S_lower':str(TL-IU),'certificate_S_upper':str(TU-IU),'fixed_g105_trueI_S_lower':str(TL-IU),'fixed_g105_trueI_S_upper':str(TU-IL),'certificate_S_exact_fourth_power_margin':str(32**4*g**5-5**5*IU**4),'certificate_S_positive':32**4*g**5>5**5*IU**4,'all_totals_units':dict(tot),'physical_levels':{k:dict(v) for k,v in levels.items()},'index_multiplicity_classes':{k:dict(v) for k,v in shapes.items()},'weighted_selfgap_concentration':concentration,'top64_physical_selfgaps':rank[:64],'optimistic_fixed_current_M_zero_variance_family':{'is_actual_I_upper_bound':False,'assumptions':'Current saved M_U at every old/new node fixed; same E_upper; variance term hypothetically zero. Any legitimate mean-variance self upper with V_U>=0 is >= volume*M_U^q before oldself intersection. Min/sum hierarchy is monotone.','hypothetical_UH_lower':str(F(oL,N)),'hypothetical_UH_upper':str(F(oU,N)),'maximum_UH_reduction_upper':str(F(up-oL,N)),'hypothetical_budget_I_lower':str(I0L),'hypothetical_budget_I_upper':str(I0U),'hypothetical_S_lower':str(TL-I0U),'hypothetical_S_upper':str(TU-I0L),'strict32_route_excluded_even_at_zero_variance':5**5*B0L**5>=32**4*g**5,'scope':'Only fixed-current-M,same-fixed-E mean-variance certificate family; not actual model/cost/C5 exclusion.'},'time_envelope_loss_separately_identified':False,'time_loss_reason':'73 bands are exact H endpoint differences at each anchor; RU applied to each mass then summed. Current certificate has no independent true box band means/variances, so bandwise worstcase sum and spatial envelope/mean slack cannot be uniquely separated. The global E_upper/Minkowski gap is an upper-budget allowance, not measured true shorttime error.','source_inputs':inputs,'new_H_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-t,'source_sha256':sha(__file__)};raw=S/'gap_components.json.gz';raw.write_bytes(gzip.compress(json.dumps(records,sort_keys=True,separators=(',',':')).encode(),mtime=0));out['external_components']={'path':str(raw),'sha256':sha(raw),'rows':len(records)};(H/'gap_analysis.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ('status','all_nodes','physical_leaves','certificate_S_positive','time_envelope_loss_separately_identified','pid','seconds')},indent=2));print('S',float(TL-IU),float(TU-IU),'fixedI_Supper',float(TU-IL));print('zeroVar','I',float(I0L),float(I0U),'UHrecoverUpper',float(F(up-oL,N)),'excluded',out['optimistic_fixed_current_M_zero_variance_family']['strict32_route_excluded_even_at_zero_variance'])
if __name__=='__main__':main()
