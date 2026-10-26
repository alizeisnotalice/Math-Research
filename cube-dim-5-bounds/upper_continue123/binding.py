"""Read-only binding to the already reviewed098 H128 numerical kernel."""
from pathlib import Path
import hashlib,importlib.util
KERNEL_PATH=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_mixed098/epsilon_kernel.py')
KERNEL_SHA='128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0'
MODEL='n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
assert hashlib.sha256(KERNEL_PATH.read_bytes()).hexdigest()==KERNEL_SHA
spec=importlib.util.spec_from_file_location('readonly098_H128',KERNEL_PATH)
api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)
