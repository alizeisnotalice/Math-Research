"""One incremental payment point, from preserved R072 arithmetic; no prefix/table replay."""
from fractions import Fraction as F
from pathlib import Path
import datetime,hashlib,json,os
root=Path(__file__).resolve().parent
hashof=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
oldpath=root/'verification.json';old=json.loads(oldpath.read_text());v=old['constants']
assert hashof(oldpath)=='a653d150ce195cd74212a251d6a7017133ed5b17d1610ccdfb969955a23ee971'
eR=F(1,2**20); eD=eR; epsRet=F(8,25); epsS=F(37,500); epsT=F(1,5000)
ulo=F(v['Ulo']);M=F(v['M']);mL=F(311,1000);gamma=F(39,200);c=F(19,200)
cap=F(v['cap_bound_S64']);tau=F(v['tau_bound_S64']);R=F(v['scalar_remainder_upper']);kept=F(v['normal_minus_BE']);price_upper=F(999,10000)
b=2*(1-epsRet)*(1-eR)*ulo-epsS-epsT;Gamma=mL*b/M;budget=(1-eD)*Gamma
checks={'all_old_common_arithmetic_passed':all(old['common_checks'].values()),'conditional_retention_budget':kept>1-epsRet==F(17,25),'cap_payment':cap<epsS,'tau_payment':tau<epsT,'energy_budget':0<gamma<budget<1,'scalar_slack_identity':gamma-price_upper-c==F(1,10000),'scalar_full_payment':R<F(1,10000),'display_improves_original':c>F(9,100)}
passed=all(checks.values())
source=Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/upper/round070/second_profile/PROOF.md')
expected='4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8';assert hashof(source)==expected
margins={'conditional_retention':(kept-(1-epsRet))/kept,'cap':(epsS-cap)/epsS,'tau':(epsT-tau)/epsT,'energy':(budget-gamma)/budget,'scalar':(F(1,10000)-R)/F(1,10000)}
receipt={'status':('PASS' if passed else 'FAIL')+': conditional incremental arithmetic only; not Cn','actual_pid':os.getpid(),'actual_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'check_count':len(checks),'parameters':{'epsilon_Ret':str(epsRet),'epsilon_S':str(epsS),'epsilon_T':str(epsT),'gamma':str(gamma),'display_c_conditional':str(c)},'Gamma':str(Gamma),'gamma_budget':str(budget),'budget_minus_gamma':str(budget-gamma),'conditional_retention_lower':str(kept),'scalar_slack':'1/10000','normalized_margins':{k:str(q) for k,q in margins.items()},'minimum_normalized_margin':str(min(margins.values())),'R070_local_source':{'path':str(source),'sha256':hashof(source),'scope':'actual new source binding; LOCAL independently being reviewed, not a BE/full-operator theorem'},'old_receipt_sha256':hashof(oldpath),'own_sources':{name:hashof(root/name) for name in ['RETENTION_REGISTRATION.md','verify_retention_point.py']},'no_prefix_replay':True,'no_12row_replay':True,'BE_source':'missing qualified theorem and constant; remains conditional','no_operator_activation':True}
(root/'retention_verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'actual_pid':os.getpid(),'checks':checks,'Gamma':str(Gamma),'budget_minus_gamma':str(budget-gamma),'display_c_conditional':str(c)},indent=2))
