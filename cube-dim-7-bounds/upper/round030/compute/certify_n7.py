#!/usr/bin/env python3
"""Reproduce exact separation from the previously certified seven-dimensional point."""
import sys,json,time,resource,hashlib
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
root=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(root/'revision_v3/numeric'))
from verify import certify,finite_certify
start=time.monotonic(); cpu=time.process_time()
old=next(r for r in json.loads((root/'revision_v3/numeric/table.json').read_text()) if r['n']==7)
t=Q(808656059921,125000000000); U=Q(24554960525964662,10**15)
A,bl,bu,fl,fu=certify(7,t)
Aa,delta,Eup,tau,finite=finite_certify(7,t,bu)
oldlo=Q(old['F_lower'])
assert fl<=fu<U<oldlo and finite<U
out={'n':7,'t':str(t),'F_lower':str(fl),'F_upper':str(fu),'U':str(U),'old_F_lower':str(oldlo),'strict_gap_old_lower_minus_new_upper':str(oldlo-fu),'finite_a':'2^256','finite_bound_below_U':True,'comparison':'F_lower <= F_upper < U < old_F_lower and finite_bound < U','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':__import__('os').getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('certificate_n7.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
