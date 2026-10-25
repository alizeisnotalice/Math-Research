from pathlib import Path
from fractions import Fraction as Q
import importlib.util,json,hashlib,time
from flint import arb,ctx
r=Path(__file__).resolve().parents[3];code=r/'upper/round002/anchors/truncated_anchor_batch.py'
spec=importlib.util.spec_from_file_location('old_formula',code);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);ctx.prec=384
point=(Q(1),Q(751,1000),Q(751,1000),Q(749,1000),Q(749,1000),Q(1,1000));start=time.monotonic();v,meta=m.h_epsilon(point,Q(1,100))
src=r/'upper/round004/local_box/arb384_results.json';row=next(x for x in json.loads(src.read_text())['outputs'] if x['index']==2);stored=arb(row['H_epsilon_arb_384']);assert 0 in v-stored
out={'status':'PASS_second_formula_corner_crosscheck','corner_index':2,'value':v.str(95,more=True),'difference':(v-stored).str(95,more=True),'seconds':time.monotonic()-start,'metadata':meta,'formula_source_sha256':hashlib.sha256(code.read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scope':'second explicit-real-derivative DD implementation at one new corner; shared audited mathematical fan'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
