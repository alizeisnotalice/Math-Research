"""Exact finite smooth-input existence check; no huge atom generation."""
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path

root = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument('--input', type=Path, default=root/'lower/n3/round006_structural_input.json')
ap.add_argument('--audit', type=Path, default=root/'audit/012_results.json')
ap.add_argument('--target', type=Path, default=root/'audit/007c_results.json')
ap.add_argument('--save', type=Path)
args = ap.parse_args()
inp = json.loads(args.input.read_text())
cert = json.loads(args.audit.read_text())
old = json.loads(args.target.read_text())
assert hashlib.sha256(args.input.read_bytes()).hexdigest() == cert['input_sha256']
p = list(map(F,inp['periods']))
mass = sum(map(F,inp['weights']))
assert len(inp['weights']) == 3 and mass == F(cert['mass_per_index_triple'])
tau = F(1,10**14); eta = tau/4; N = 10**8
assert 2*tau < F(cert['minimum_inner_side'])
assert all(tau < F(r['radius_lower']) for r in cert['radius_records'])
loss = 2*len(inp['keys'])*tau*sum(p[j]*p[k] for j in range(3) for k in range(j+1,3))
ratio = (F(cert['ratio_interval'][0])-loss/mass)*F(2*N+1,2*N+5)**3
gain = ratio-F(old['ratio_interval'][1])
assert cert['b'] == 2 and gain > 0
out = {'status':'PASS','input_sha256':cert['input_sha256'], 'audit_sha256':hashlib.sha256(args.audit.read_bytes()).hexdigest(),
       'N':N,'tau':str(tau),'eta':str(eta),'mass':str(mass),'union_loss_upper':str(loss),
       'finite_ratio_lower':str(ratio),'gain_above_target_full_upper':str(gain),
       'target_sha256':hashlib.sha256(args.target.read_bytes()).hexdigest(),
       'display_only':float(ratio),'gain_display_only':float(gain),
       'scope':'finite smooth-input existence parameters; no atom generation'}
if args.save: args.save.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
