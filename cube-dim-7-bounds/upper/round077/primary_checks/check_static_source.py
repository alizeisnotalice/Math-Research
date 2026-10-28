from pathlib import Path
from collections import Counter
import re,hashlib,json,os
p=Path('/Users/zhengzhihao/Desktop/T/cube-cover-manuscript/cube-dim-7-bounds/paper/main.tex')
s=p.read_text();h=hashlib.sha256(p.read_bytes()).hexdigest()
assert h=='989d28ff83bd51f2c2dfe1822d9ecb59dd71502e5deb361db845cd87359edc9d'
labels=re.findall(r'\\label\{([^}]+)\}',s); refs=re.findall(r'\\(?:eqref|ref|autoref)\{([^}]+)\}',s)
bib=set(re.findall(r'\\bibitem(?:\[[^]]*\])?\{([^}]+)\}',s));cites=[x.strip() for group in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}',s) for x in group.split(',')]
dups=[x for x,n in Counter(labels).items() if n>1]; missing=sorted(set(refs)-set(labels));missing_cites=sorted(set(cites)-bib)
secs=re.findall(r'\\section\{([^}]+)\}',s)
checks={'no_duplicate_labels':not dups,'all_refs_defined':not missing,'all_citations_defined':not missing_cites,'four_section_structure':secs==['引言','预备知识','主要证明','应用'],'new_lower_threshold_present':'L_*=8192' in s and '47/250' in s,'stronger_large_domain_present':'2^{30}' in s and '23/100' in s,'key_source_parameter':'11569/12800' in s and '25/256' in s,'different_rates_retained':'13627383/3790929920' in s and '1447/327680' in s,'complete_cantelli_certificate':'9717/10000' in s and '291/1000' in s,'no_discouraged_missing_factor_sentence':'不可遗漏' not in s}
assert all(checks.values()),(checks,dups,missing,missing_cites)
rec={'status':'PASS','pid':os.getpid(),'source_sha256':h,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'labels':len(labels),'references':len(refs),'citations':len(cites),'checks':checks,'scope':'Independent static source and cross-reference checks; complements rather than substitutes for whole-text proof review or visual PDF review.'}
Path(__file__).with_name('check_static_source.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2)+'\n');print(json.dumps(rec,ensure_ascii=False))
