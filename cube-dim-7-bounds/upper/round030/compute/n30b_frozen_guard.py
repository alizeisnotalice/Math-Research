#!/usr/bin/env python3
"""Correct frozen parameters: sigma=5/(3t), guard R=7/4>t*sigma."""
import sys,os,json,time,hashlib,resource
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
from n30b import certify,Pbounds,Bbounds,sqrtbounds
from verify import log_bounds
start=time.monotonic();cpu=time.process_time();n=7;t=Q(3391749,524288);s=Q(5,3)/t;R=Q(7,4);assert 1<t*s<R<t
cert=certify(n,t,s);B,Bu=Bbounds(n,t);D=Q(cert['D']);assert D==Q(61741,2187)/t**7
A=(1+1/t)**n;tau=sqrtbounds(A*(Bu-D))[1];pl,_=Pbounds(n,t,R/t);coarse=B-1-n*log_bounds(R)[1]
assert pl-D>tau and coarse-D>tau
out={'n':n,'t':str(t),'sigma':str(s),'t_sigma':str(t*s),'guard_R':str(R),'D':str(D),'coarse_P_lower_minus_D_gt_tau_upper':True,'guard_P_lower_minus_D_gt_tau_upper':True,'coarse_margin_decimal':format(float(coarse-D-tau),'.12g'),'guard_margin_decimal':format(float(pl-D-tau),'.12g'),'explanation':'R>t sigma is a valid guard; the prior mismatch claim was incorrect.','pid':os.getpid(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('n30b_frozen_guard.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
