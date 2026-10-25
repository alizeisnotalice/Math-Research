from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,os
p=Path(__file__).resolve().parents[1];began=time.monotonic();source=p/'upper_sip015/exact_witness_results.json';d=json.loads(source.read_text())['witness'];atoms=[(F(x['z']),F(x['d']),F(x['weight'])) for x in d['atoms']];assert len(atoms)==6
assert all(0<z<=1 and dd==z**20 and w>0 for z,dd,w in atoms)
corner=p/'astra_upper006/verification.json';mom=json.loads(corner.read_text())['moments'];targets=[F(1)]+[F(mom[str(k)][0]) for k in range(1,5)]
residuals=[sum((w*dd**k for z,dd,w in atoms),F(0))-targets[k] for k in range(5)];Q=F(926346056693,500000000000);residuals.append(sum((w/z**4 for z,dd,w in atoms),F(0))-Q);assert all(x==0 for x in residuals)
L=sum((w/z**3 for z,dd,w in atoms),F(0));assert L>F(3,2) and L==F(d['objective_exact'])
out={'status':'PASS_primary_exact_new_abstract_moment_obstruction','positive_atoms':6,'moment_residuals':[str(x) for x in residuals],'objective_exact':str(L),'objective_minus_three_halves':str(L-F(3,2)),'Q_exact':str(Q),'proof':'For d=z^20, integrate any valid majorant against these six positive weights. Moment equality gives proposed budget>=L>3/2. No actual determinant distribution or C5 lower bound is asserted.','independence':'No producer LP/solver/verification modules imported; Python Fraction backend shared.','source_sha256':{str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [source,corner,p/'upper_I027/manifest.json',p/'audit_I028/report.md',Path(__file__),p/'round009/check_sip_registration.md']},'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began}}
(p/'round009/check_sip.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: all six exact moment equations, positive weights, objective >3/2; no C5 result')
