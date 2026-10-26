"""Exactly qualified123 direct EPS128 PF core;read-only,no H/log on import."""
import importlib.util,hashlib
from pathlib import Path
p=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-audit123-sol/cube-dim-5-bounds/audit_g123_sol/pf128.py')
assert hashlib.sha256(p.read_bytes()).hexdigest()=='4e0de0638e44ac15dd239bfc8f5bd33a7e0e7034783f32573eb3c598379f3713'
spec=importlib.util.spec_from_file_location('readonly123_directPF',p)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
values=base.values
dyadic=base.dyadic
logarithm=base.logarithm
S=base.S
