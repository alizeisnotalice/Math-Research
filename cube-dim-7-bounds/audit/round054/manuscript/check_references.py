"""R054 source identity/reference check, not author mathematical validation."""
from pathlib import Path
import hashlib,json,re,os,sys
ROOT=Path('/Users/zhengzhihao/Desktop/T')
OUT=Path(__file__).parent
p=ROOT/'cube-cover-manuscript/cube-dim-7-bounds/paper/main.tex'
s=p.read_text();labels=re.findall(r'\\label\{([^}]+)\}',s)
refs=re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',s)
cites=[item for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}',s) for item in group.split(',')]
items=re.findall(r'\\bibitem\{([^}]+)\}',s)
undefined=sorted(set(refs)-set(labels));duplicates=sorted({x for x in labels if labels.count(x)>1})
assert not undefined and not duplicates and not(set(cites)-set(items))
sources={
 'manuscript':p,
 'initial_manuscript':ROOT/'cube-cover-manuscript/cube-dim-7-bounds/upper/round054/manuscript_notes/private/main_a05.tex',
 'definition_delta':ROOT/'cube-cover-manuscript/cube-dim-7-bounds/upper/round054/manuscript_notes/kernel_definition_delta.patch',
 'family':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md',
 'family_audit':ROOT/'cube-r030-soft/cube-dim-7-bounds/audit/round051/flexible_budget/proof_review.md',
 'raw_threshold':ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round053/explicit_candidate_threshold/THRESHOLD_PROOF.md',
 'simple_threshold':ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round055/simple_tuple_threshold/PROOF.md',
 'raw_independent_check':ROOT/'cube-r030-uniform-audit/cube-dim-7-bounds/audit/round053/threshold/microcheck.json',
 'simple_independent_check':ROOT/'cube-r030-uniform-audit/cube-dim-7-bounds/audit/round055/threshold/microcheck.json'}
receipt={'pid':os.getpid(),'python':sys.version,'line_count':len(s.splitlines()),
 'labels':len(labels),'references':len(refs),'undefined':undefined,'duplicate_labels':duplicates,
 'undefined_citations':sorted(set(cites)-set(items)),
 'hashes':{k:hashlib.sha256(v.read_bytes()).hexdigest() for k,v in sources.items()},
 'sources':{k:str(v) for k,v in sources.items()},
 'scope':'Read-only source/reference lock. No 51-check author program, old mathematical rerun, production run or TeX compilation.'}
(OUT/'sourcehash.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS references/citations; current SHA',receipt['hashes']['manuscript'],'pid',os.getpid())
