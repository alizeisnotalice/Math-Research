"""New R048 endpoint checks and exact preservation, not theorem numerics."""
from pathlib import Path
from fractions import Fraction as F
import difflib, hashlib, json, os, sys

ROOT = Path('/Users/zhengzhihao/Desktop/T')
WRITER = ROOT/'cube-cover-manuscript/cube-dim-7-bounds'
OUT = Path(__file__).parent
sources = {
 'main': WRITER/'paper/main.tex',
 'accepted_base': WRITER/'upper/round048/manuscript_notes/private/main_94c5.tex',
 'family': ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md',
 'local': ROOT/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round045/one_coordinate/LEMMA.md',
 'comparison': ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round045/comparison/PROOF.md',
 'comparison_audit': ROOT/'cube-r040-star-lp/cube-dim-7-bounds/audit/round045/comparison/AUDIT.md',
}
hashes = {k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in sources.items()}
assert hashes['main']=='ecd7fd7487083b7eb6318a0a94f969f9d291b9d4df443cc3374497a7f103b248'
assert hashes['accepted_base']=='94c5e61d3bf86b16736e328ffe3e412d46df957dd2b8f8760939edc8831e8692'
p2=lambda exponent:F(2)**exponent
d=p2(-24); sigma=32*d*d; delta=d/65536
checks={
 'endpoint_sigma':sigma==p2(-43),
 'endpoint_budget':delta==p2(-40),
 'scalar_quarter_at_closed_endpoint':2*sigma==delta/4,
 'bridge':F(1280,2**20)==F(5,4096)<F(1,8),
 'cap_exponent_margin':12*d<F(1,64),
 'full_cost_tail':p2(-479)+p2(-393)+p2(-270)<p2(-41)-p2(-42),
 'R042_actual_cost_lower_gap':p2(-238)+p2(-111)<p2(-41),
 'old_d_actual_cost_lower_gap':p2(-238)+p2(-46)<p2(-41),
}
assert all(checks.values())
old=sources['accepted_base'].read_text().splitlines()
new=sources['main'].read_text().splitlines()
ops=[{'operation':tag,'old_inclusive':[i+1,j],'new_inclusive':[k+1,l]}
     for tag,i,j,k,l in difflib.SequenceMatcher(None,old,new,autojunk=False).get_opcodes() if tag!='equal']
# Entire retained proof after the soft insertion up to the application insertion.
# Exact line mapping was inspected independently; assert all old R042 and old-energy body.
assert old[539:1282]==new[900:1643]
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,
 'hashes':hashes,'sources':{k:str(p) for k,p in sources.items()},'changes':ops,
 'scope':'Exact new endpoint arithmetic plus byte/line preservation. No author verifier, 13-check rerun, numerical search, cap integration, old baseline, or compilation.'}
(OUT/'microcheck.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'new exact checks; retained old proof body identical; pid',os.getpid())
