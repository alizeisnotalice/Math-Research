#!/usr/bin/env python3
"""New original-symbol diagnostics. Not interval arithmetic or a weak counterexample."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, math, time
import numpy as np

HERE=Path(__file__).resolve().parent
PREFIX="ordered_capped_test_energy"
REG=HERE/(PREFIX+"_registration_20261007.json")
START=time.time()
checks=[]
def check(label, ok):
    checks.append({"label":label,"pass":bool(ok)})
    if not ok: raise AssertionError(label)
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
registration=json.loads(REG.read_text())
check("preregistered",registration["status"]=="registered_before_execution")

# Exact rational constant chain; transcendental bounds are proved in the companion note.
check("exp9 via e>8/3",(F(8,3)**9)>4097)
check("B4096>450",F(4096,9)-1>450)
check("g_bound<1/100",1/(1+F(450,4))<F(1,100))
check("sqrt((1-g)/2)>7/10",F(99,200)>F(7,10)**2)
check("flat error<1/200",F(5,1998)<F(1,200))
check("flat d>1",(F(3,2)*F(7,10)-F(1,200))**2>1)
check("late a<7/10",F(2,3)+F(1,100)<F(7,10))
check("late m8<1/16",F(7,10)**8<F(1,16))
check("late t error<=1/6",F(5,2)*F(1,16)/(1-F(1,16))==F(1,6))
check("gap probability variance>=3/20",F(1,5)*(1-F(1,4))==F(3,20))
check("gap step d>4/5",2*F(99,100)*F(9,4)*F(3,20)>F(4,5)**2)
check("gap perturbed energy>1/3",(F(4,5)-F(1,6))**2>F(1,3))
check("gap M>lambda",1-F(1,16)>F(1,2))
check("paired coefficient<3/2",2*F(100,99)**2<F(3,2)**2)
for n in (8,32,128):
    check(f"n{n}:late power",F(7,10)**n<F(1,16))
    check(f"n{n}:gap time energy n/6",F(n,2)*F(1,3)==F(n,6))
    check(f"n{n}:rational epsilon cap",F(1,{8:1000,32:10000,128:100000}[n])<=F(1,1000))

def B(v):
    v=np.asarray(v,dtype=float)
    out=np.zeros_like(v)
    mask=v>0
    out[mask]=v[mask]/np.log1p(v[mask])-1
    return out
def simpson(a):
    panels=len(a)-1
    assert panels%2==0
    return float((a[0]+a[-1]+4*a[1:-1:2].sum()+2*a[2:-1:2].sum())/(3*panels))
def entropy(Bmean,A,m,lam):
    a=lam+Bmean
    return Bmean-lam*math.log((a+math.sqrt(a*a-A*A*m*m))/(2*lam))

rows=[]
saved=[]
for spec in registration["rounds"]:
    n,K=spec["n"],spec["K"]
    eps=float(F(spec["epsilon"]))
    baseB=float(B(np.array([K*K]))[0])
    g1=1/(1+baseB)
    mZ=(math.exp(-1)+(1-math.exp(-1))*g1)**n
    for variant,A in (("flat",eps),("gap",1.0)):
        previous=None
        for N,panels in zip(registration["phase_grids"],registration["time_panels"]):
            phase=2*np.pi*(np.arange(N)+.5)/N
            co=np.cos(phase)
            f=1+A*co
            M=np.where(co>0,f,1+A*mZ*co)
            selected=np.ones(N,dtype=bool) if variant=="flat" else (co<=-.75)
            active=selected & (co<0)
            lam=.5
            theta=np.where(selected,lam/M,0.0)
            target=float(lam*selected.mean())
            initial=float(np.mean(theta*f))
            expected_flow=target-initial
            z=np.linspace(0,1,panels+1)
            freq=np.arange(N//2+1,dtype=float)
            freqB=B((K*freq)**2)
            weights=np.full(N//2+1,2.0)
            weights[0]=weights[-1]=1.0
            values={k:np.zeros(panels+1) for k in
                    ("test","scaled_test","dissipation","commutator","positive","flow","pair","B_mass")}
            minB=math.inf
            max_t=0.0
            maxres=0.0
            for iz,zz in enumerate(z):
                d=math.exp(-float(zz))
                a=d+(1-d)*g1
                m=a**n
                gd=1/(1+d*baseB)
                mp=-n*(1-gd)*m
                maxres=max(maxres,abs(mp-(-n*d*(1-g1)*a**(n-1))))
                v=1+A*m*co
                # Endpoint z=1 has a zero strict gate; changing one time point
                # must not alter its integral. Integrate the left continuous
                # value at z=1, and disclose this endpoint convention.
                t=np.where(active,(lam+v)/M,0.0)
                max_t=max(max_t,float(t.max()))
                H=np.where(active,lam/M,0.0)
                logv=np.log1p(v/lam)
                tf=np.fft.rfft(t)/N
                lf=np.fft.rfft(logv)/N
                multip=1/(1+d*freqB)
                jlog=n*np.fft.irfft((multip-1)*lf*N,n=N)
                jv=A*mp*co
                b=lam/(lam+v)*jv-lam*jlog
                minB=min(minB,float(b.min()))
                energy=2*n*(float(np.mean(t*t))-float(np.sum(weights*multip*(tf.real**2+tf.imag**2))))
                aa=lam*energy/2
                root=math.sqrt((lam+1)**2-(A*m)**2)
                dd=lam*n*(1-gd)*A*A*m*m/(root*(lam+1+root))
                values["test"][iz]=energy
                values["scaled_test"][iz]=aa
                values["dissipation"][iz]=dd
                values["commutator"][iz]=float(lam*np.mean(t*jlog))
                values["positive"][iz]=float(np.mean(t*b))
                values["flow"][iz]=float(np.mean(H*jv))
                values["pair"][iz]=math.sqrt(max(0,dd*aa))
                values["B_mass"][iz]=float(b.mean())
            integrals={k:simpson(v) for k,v in values.items()}
            closed_drop=entropy(1,A,1,lam)-entropy(1,A,mZ,lam)
            global_pair=math.sqrt(max(0,closed_drop*integrals["scaled_test"]))
            q=target
            lower=n if variant=="flat" else n/6
            check(f"{n}/{variant}/{N}:source nonnegative",float(f.min())>=-1e-13)
            check(f"{n}/{variant}/{N}:M>lambda on selection",bool(np.all(M[selected]>lam)))
            check(f"{n}/{variant}/{N}:t<=2",max_t<=2+1e-13)
            check(f"{n}/{variant}/{N}:energy analytic lower",integrals["test"]>=lower-1e-7)
            check(f"{n}/{variant}/{N}:positive B screen",minB>=-1e-8)
            check(f"{n}/{variant}/{N}:symbol derivative",maxres<1e-10)
            check(f"{n}/{variant}/{N}:paired <=1.5W",integrals["pair"]<1.5+1e-5)
            check(f"{n}/{variant}/{N}:flow decomposition",
                  abs(integrals["flow"]-integrals["positive"]-integrals["commutator"])<1e-9)
            check(f"{n}/{variant}/{N}:comm <=time pair",
                  abs(integrals["commutator"])<=integrals["pair"]+1e-8)
            check(f"{n}/{variant}/{N}:B mass vs dissipation",
                  abs(integrals["B_mass"]-integrals["dissipation"])<1e-7)
            fn=HERE/f"{PREFIX}_n{n}_{variant}_N{N}_20261007.npz"
            np.savez_compressed(fn,z=z,**values)
            saved.append({"file":fn.name,"sha256":sha(fn)})
            row={"n":n,"K":K,"epsilon":spec["epsilon"],"variant":variant,
                 "phase_grid":N,"time_panels":panels,"lambda":lam,"mass_W":1,
                 "selected_volume":float(selected.mean()),"lambda_volume":q,
                 "initial_paid":initial,"exact_flow":expected_flow,"mZ":mZ,
                 "closed_entropy_drop":closed_drop,"integrals":integrals,
                 "global_cauchy":global_pair,"energy_lower_analytic":lower,
                 "min_B_float":minB,"max_t_float":max_t,"symbol_derivative_residual":maxres,
                 "entropy_time_error":abs(integrals["dissipation"]-closed_drop),
                 "winner_flow_time_error":abs(integrals["flow"]-expected_flow)}
            if previous:
                row["nested_differences"]={k:abs(integrals[k]-previous["integrals"][k]) for k in integrals}
                row["volume_nested_difference"]=abs(row["selected_volume"]-previous["selected_volume"])
            rows.append(row)
            previous=row
            print(f"n={n} {variant} phase={N} panels={panels} "
                  f"T={integrals['test']:.9g} pair={integrals['pair']:.9g} "
                  f"comm={integrals['commutator']:.9g}",flush=True)

result={"status":"PASS_DIAGNOSTIC_AND_EXPLICIT_RATIONAL_CHAIN",
        "registration_sha256":sha(REG),"script_sha256":sha(Path(__file__)),
        "checks":checks,"check_count":len(checks),"rows":rows,"saved_arrays":saved,
        "runtime_seconds":time.time()-START,
        "limitations":["No interval certification of FFT/Simpson",
                       "Strict gate at terminal time integrated using left endpoint representative",
                       "Periodic finite Fourier source, plus analytical finite-L1 large-box lifting",
                       "Not a weak endpoint counterexample or actual geom-qualified sample",
                       "Time-local pairing upper shown only for these single-phase sources",
                       "No general source-once paired-budget theorem or cube fee"]}
out=HERE/(PREFIX+"_results_20261007.json")
out.write_text(json.dumps(result,indent=2)+"\n")
print(f"Saved {out.name}: {len(checks)} checks, {len(rows)} rows, runtime {result['runtime_seconds']:.2f}s",flush=True)
