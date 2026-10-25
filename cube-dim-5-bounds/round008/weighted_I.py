"""Finite type-weight search, then outward ball comparison; I020 theory pending audit."""
from pathlib import Path
from fractions import Fraction as F
from scipy.optimize import minimize
from scipy.special import gammaln
import numpy as np
from flint import arb,ctx
import math,json,hashlib,time,os,sys
P=Path(__file__).resolve().parent;began=time.monotonic();ctx.prec=256
w0=np.array([.5,.375,.125]);rs=np.array([8**.25,1+1/math.sqrt(2),(7+8*math.sqrt(2)+3*math.sqrt(3))/2**3.75])
def fun(z):
 w=np.array([z[0],z[1],1-z[0]-z[1]])
 if np.any(w<=0):return 1e10
 rho=np.dot(w,[.75,.875,1.])
 return math.log(2)+gammaln(1.75)+.15*np.sum(w*np.log(w/w0))+5/rho*(gammaln(1+rho/5)-gammaln(1+.15*rho))-.2*np.dot(w,np.log(rs))
def ar(x):
 x=F(x);return arb(x.numerator)/x.denominator
def calc(w):
 rho=sum((a*b for a,b in zip(w,[F(3,4),F(7,8),F(1)])),F(0));rr=ar(rho)
 bases=list(map(ar,[F(1,2),F(3,8),F(1,8)]));ww=list(map(ar,w));D=sum((x*(x/y).log() for x,y in zip(ww,bases)),arb(0))
 R=[ar(8)**ar(F(1,4)),1+1/ar(2).sqrt(),(7+8*ar(2).sqrt()+3*ar(3).sqrt())/ar(2)**ar(F(15,4))]
 J=2*ar(F(7,4)).gamma()*(ar(F(3,20))*D).exp()*((1+rr/5).gamma()/(1+ar(F(3,20))*rr).gamma())**(5/rr)
 for x,r in zip(ww,R):J*=r**(-x/5)
 c=3/(ar(2)**ar(F(15,2))*arb.pi()**ar(F(5,2)));I=c**ar(F(10,7))*arb.pi()**2/60*J**ar(F(10,7))
 return J,I,rho
opt=minimize(fun,w0[:2],method='SLSQP',bounds=[(1e-6,1-1e-6)]*2,constraints=[{'type':'ineq','fun':lambda z:1-1e-6-z.sum()}],options={'maxiter':40,'ftol':1e-12})
assert np.isfinite(opt.x).all() and time.monotonic()-began<120
candidates=[]
for mode in ('floor','nearest'):
 n=[math.floor(x*10000) if mode=='floor' else round(x*10000) for x in opt.x];w=[F(n[0],10000),F(n[1],10000),1-F(sum(n),10000)]
 if min(w)>0 and w not in candidates:candidates.append(w)
defJ,defI,defrho=calc([F(1,2),F(3,8),F(1,8)]);old=arb(json.loads((P.parent/'round006/joint_cost.json').read_text())['I_ball'])
records=[]
for w in candidates:
 J,I,rho=calc(w)
 records.append({'weights':[str(x) for x in w],'rho':str(rho),'J_ball':J.str(70,more=True),'I_ball':I.str(70,more=True),'new_minus_default_I_ball':(I-defI).str(70,more=True),'strictly_improves_default_I':bool(I<defI),'strictly_improves_round006_I':bool(I<old)})
out={'status':'conditional_complete_I_candidate_pending_I020_theory_audit','formula':'J=2Gamma(7/4) exp((3/20)KL(w||theta0)) [Gamma(1+rho/5)/Gamma(1+3rho/20)]^(5/rho) prod R_i^(-w_i/5)','default':{'weights':['1/2','3/8','1/8'],'rho':str(defrho),'J_ball':defJ.str(70,more=True),'I_ball':defI.str(70,more=True)},'optimizer_diagnostic_only':{'success':bool(opt.success),'message':str(opt.message),'iterations':int(opt.nit),'w12':list(map(float,opt.x)),'logJ':float(opt.fun)},'candidates':records,'not_a_global_optimality_certificate':True,'C5_upper_changed':False,'source_sha256':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [Path(__file__),P/'weighted_I_registration.md']},'runtime':{'pid':os.getpid(),'python':sys.executable,'seconds':time.monotonic()-began,'arb_bits':ctx.prec}}
(P/'weighted_I.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
