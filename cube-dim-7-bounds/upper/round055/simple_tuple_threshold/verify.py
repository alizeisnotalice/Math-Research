"""Standalone exact certificate for one specified grid member; not a search."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math

root=Path(__file__).resolve().parent
L0=2**27;Lr=2**32
A=F(9);S=F(19);d=F(1,1024);e=F(1,1024);eA=F(1,64);e0=F(1,2**20)
P=F(3579,4096);qmin=1-3*d/4
ul=P/(2*(1+e));uh=1/(P*(1-d/2)*(1-e));Mstar=max(F(4),2*uh/qmin)
m=(1-eA)*d/(2*(1+3*d/4));bracket=(1-3/A)*ul-2*e
Gamma=m*bracket/Mstar;gamma=(1-e)*Gamma;sigma=S*d*d;price=sigma/(1-sigma)
craw=gamma-price;c=F(1,2**16)
assert m==F(63,131168) and bracket==F(151679,524800)
assert Gamma==F(9555777,275347865600)
assert gamma==F(9775559871,281956214374400)
assert craw==F(4893063658542547,295647162275777740800)
h0=e0/(L0*d);nmin=F(L0**2,2);tmin=F(L0,2)
rmin=1+d/4;rmax=1+3*d/4;qb=1-3*d/4-h0;Rmin=rmin-2*h0
rhomin=tmin*qb;zt=F(2,L0*nmin*nmin)
mean=math.factorial(8)*(24*d**3/F(L0**6)+4*d/F(L0**7))
var=48*math.factorial(8)*d**3/F(L0**6);Jerr=4*d*math.factorial(8)/F(L0**7)
eta=d/32-A*d*d;nu=(S/2-A)*d*d
cap=F(24*2**27,2944**4)
tau=F(2**51*math.factorial(64),2**384)
rem=F(1,2**53)+F(math.factorial(16),2**176)+F(2**15*math.factorial(16),2**160)
private=Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds')
rawdir=private/'upper/round053/explicit_candidate_threshold'
rawreceipt=json.loads((rawdir/'verification.json').read_text())
rawc=F(rawreceipt['exact_fields']['c_raw'])
rawrem=F(1,2**63)+F(math.factorial(16),2**256)+F(2**17*math.factorial(16),2**240)
checks={
 'registered_member':A*8==72 and S-2*A==1 and d*2**20==1024,
 'domain_a_s':A>3 and S>2*A,
 'domain_d_sigma':0<d<=F(1,16) and d<F(1,32)/A and 0<sigma<F(1,2),
 'all_error_domains':0<e<1 and 0<eA<1,
 'positive_bracket':bracket>0,
 'exact_M_branch':Mstar==4,
 'craw_above_chosen_c':craw>c,
 'h_over_d':h0/d==F(1,2**27)<F(1,8),
 'zeta_t_over_h':zt/h0<F(1,2),
 'shift_below_half_d':h0+zt<d/2,
 'q_bridge_lower':qb>F(1,2),
 'q0_upper':1-d/4+h0<1,
 'R_bridge':Rmin>1,
 'rho_geometry':rhomin>4*rmax and rmax<2,
 'g_range':0<d-h0 and d-h0+zt<=d,
 'baseline_abs_mean':mean<=F(1,32),
 'tight_baseline_variance':var<=F(3,32)**2*eA/3,
 'Kw_Cheb_probability':nmin*d>=2 and F(8)/(9*nmin*d)<=eA/3,
 'target_cap_probability':nmin>=33 and F(3,4)**32<eA/3,
 'target_cap_geometry':1/rmax-F(2,L0)>F(1,2),
 'individual_J_precision':Jerr<=e,
 'product_exact_P0':1-F(1,8)-1280*e0==P,
 'Ld':L0*d>=1,
 'source_rate':eta==F(23,2**20)>0,
 'source_cap_epsS':cap<e,
 'tau_rate':nu==F(1,2**21),
 'tau_base_argument':nu*L0==64,
 'tau_base_prefactor_squared':(2**20*d)**2*L0**3<2**102,
 'tau_epsT':tau<e,
 'fixed_epsD_margin':Gamma-gamma==e*Gamma>0,
 'budget_positivity':0<gamma<1 and L0-27-F(4,L0)>8,
 'scalar_arguments':sigma*L0==2432>2**11 and sigma*L0/2==1216>2**10,
 'scalar_prefactor_squared':4*L0<2**30,
 'full_scalar_slack':rem<craw-c,
 'raw_full_cost_order':rawc>gamma-sigma+F(4,Lr)+rawrem,
}
assert all(checks.values()),[k for k,v in checks.items() if not v]
source=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md')
sourcehash=hashlib.sha256(source.read_bytes()).hexdigest()
rawhash=hashlib.sha256((rawdir/'THRESHOLD_PROOF.md').read_bytes()).hexdigest()
assert sourcehash=='aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6'
assert rawhash=='26adb4e31f0026637ea17d8f7943bc2933dce97b8655e121f2604b6ddeea2775'
fields={'a':A,'s':S,'d':d,'u_lo':ul,'u_hi':uh,'mass':m,'M':Mstar,
        'bracket':bracket,'Gamma':Gamma,'gamma':gamma,'sigma':sigma,'price':price,'c_raw':craw}
receipt={
 'status':'PASS: standalone exact certificate and finite threshold gates',
 'identity':'one specified registered-list member; no new search or enumeration replay',
 'R051_status':'generic source accepted by coordinator; explicit threshold subject to independent audit',
 'threshold':'log(n)>=2^27','chosen_c':str(c),
 'fields':{k:str(v) for k,v in fields.items()},'checks':checks,
 'sources':{'R051':{'path':str(source),'sha256':sourcehash},'R053':{'path':str(rawdir/'THRESHOLD_PROOF.md'),'sha256':rawhash}},
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('PROOF.md','verify.py')},
 'no_grid_replay':True,'no_production_search':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),
 'proof_sha256':receipt['source_sha256']['PROOF.md'],'c_raw':str(craw),'threshold':'L>=2^27'},ensure_ascii=False,indent=2))
