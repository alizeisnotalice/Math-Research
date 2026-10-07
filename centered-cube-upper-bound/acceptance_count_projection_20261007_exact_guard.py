#!/usr/bin/env python3
"""New finite exact count/CI guards; NOT original Rn/FIRST samples."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib, json, time

BASE = Path(__file__).resolve().parent
PREFIX = "acceptance_count_projection_20261007"
CHECKS = 0

def check(value, label):
    global CHECKS
    if not value:
        raise AssertionError(label)
    CHECKS += 1

def add(a, b):
    return [x+y for x,y in zip(a,b)]

def scale(a, t):
    return [t*x for x in a]

def had(a):
    a = list(a)
    width = 1
    while width < len(a):
        for start in range(0,len(a),2*width):
            for i in range(start,start+width):
                x,y = a[i],a[i+width]
                a[i],a[i+width] = x+y,x-y
        width *= 2
    return a

def markov(a, rho):
    """Tensor binary symmetric kernel with Walsh eigenvalue rho per axis."""
    n = len(a)
    f = had(a)
    f = [v*rho**i.bit_count() for i,v in enumerate(f)]
    return [v/n for v in had(f)]

def sensor(a,j,d,varying):
    if not varying:
        return list(a)
    bit = 1 << ((j-1) % d)
    return [(a[x]+a[x^bit])/2 for x in range(len(a))]

def generator(a,d):
    return [(d*a[x]-sum(a[x^(1<<i)] for i in range(d)))/3
            for x in range(len(a))]

def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),F(0))

def binlaw(M,p):
    return [F(comb(M,k))*p**k*(1-p)**(M-k) for k in range(M+1)]

def inversions(M,J):
    delta = F(1,20*J)
    grid = [F(i,16) for i in range(17)]
    laws = {p:binlaw(M,p) for p in grid}
    out = []
    for s in range(M+1):
        lower = F(0) if s==0 else max(
            p for p in grid if sum(laws[p][s:]) <= delta/2)
        upper = F(1) if s==M else min(
            p for p in grid if sum(laws[p][:s+1]) <= delta/2)
        check(lower<=upper, "CI interval ordered")
        out.append((lower,upper))
    return delta,out

def strf(x):
    return str(x)

def case(d,J,varying,tau_factor):
    before = CHECKS
    size = 1 << d
    u = [F(0)]*size
    u[0],u[1] = F(1),F(2)
    Omega = {0,1}
    Su = generator(u,d)
    kappa = max(-Su[x] for x in range(size) if x not in Omega)
    nu = [kappa+Su[x] if x in Omega else F(0) for x in range(size)]
    mu = [nu[x]-Su[x] for x in range(size)]
    W = sum(nu)
    check(W>0 and min(nu)>=0 and min(mu)>=0 and max(mu)<=kappa,
          "finite saturated source positivity")
    check(all(mu[x]==kappa for x in Omega) and sum(mu)==W,
          "finite saturation and mass")
    rho = [F(1)-F(j,3*J) for j in range(J+1)]
    increments = [None]+[rho[j]/rho[j-1] for j in range(1,J+1)]
    tau = tau_factor*kappa
    fields = [None]+[sensor(markov(nu,rho[j]**2),j,d,varying)
                     for j in range(1,J+1)]
    M = [max(fields[j][x] for j in range(1,J+1)) for x in range(size)]
    E = [x for x in range(size) if x not in Omega and M[x]>tau]
    winners = {x:next(j for j in range(1,J+1) if fields[j][x]==M[x]) for x in E}
    gs = [None]+[[tau/M[x] if x in winners and winners[x]==j else F(0)
                  for x in range(size)] for j in range(1,J+1)]
    qs = [None]+[markov(sensor(gs[j],j,d,varying),rho[j]) for j in range(1,J+1)]
    for j in range(1,J+1):
        check(min(qs[j])>=0 and max(qs[j])<=1, "coin range")
        check(sum(fields[j])==W and sum(qs[j])==sum(gs[j]), "receiver mass")
        check(markov(markov(nu,rho[j]),rho[j])==markov(nu,rho[j]**2),
              "P square equals K")
        check(markov(markov(nu,rho[j-1]),increments[j])==markov(nu,rho[j]),
              "positive cocycle")
        check(dot(markov(nu,rho[j]),qs[j])==dot(gs[j],fields[j]),
              "Fubini response direction")
    check(all(sum(gs[j][x] for j in range(1,J+1))<=1 for x in range(size)),
          "disjoint truewinner weights")
    # Complete joint law, state and count, with all prior coin failures preserved.
    forward = [scale(nu,1/W)]
    for j in range(1,J+1):
        advanced = [markov(a,increments[j]) for a in forward]
        nextlaw = [[F(0)]*size for _ in range(j+1)]
        for k,a in enumerate(advanced):
            for x in range(size):
                nextlaw[k][x] += a[x]*(1-qs[j][x])
                nextlaw[k+1][x] += a[x]*qs[j][x]
        forward = nextlaw
    law = [sum(a) for a in forward]
    check(min(law)>=0 and sum(law)==1,"joint law probability")
    # Independent backward polynomial recursion, conditional on initial source.
    future = [[F(1)]*size]
    for j in range(J,0,-1):
        if j<J:
            future = [markov(a,increments[j+1]) for a in future]
        b = [[F(0)]*size for _ in range(len(future)+1)]
        for k,a in enumerate(future):
            for x in range(size):
                b[k][x] += (1-qs[j][x])*a[x]
                b[k+1][x] += qs[j][x]*a[x]
        future = b
    conditional = [markov(a,increments[1]) for a in future]
    check([dot(nu,a)/W for a in conditional]==law, "independent backward law")
    check(all(sum(a[x] for a in conditional)==1 for x in range(size)),
          "conditional source law")
    mean = sum(F(k)*p for k,p in enumerate(law))
    stop = 1-law[0]
    deficit = sum(F(max(k-1,0))*p for k,p in enumerate(law))
    pair = sum(F(comb(k,2))*p for k,p in enumerate(law))
    tails = [sum(law[m:]) for m in range(1,J+1)]
    check(mean==tau*len(E)/W, "mean equals Lebesgue volume")
    check(mean==stop+deficit, "stop deficit decomposition")
    check(mean==sum(tails) and deficit==sum(tails[1:]), "count tail sums")
    check(pair==sum(F(m-1)*tails[m-1] for m in range(2,J+1)), "pair tail sum")
    # Backward union/overlap recurrence, separate from full PGF.
    c,dvec = list(qs[J]),[F(0)]*size
    for j in range(J-1,0,-1):
        vc,vd = markov(c,increments[j+1]),markov(dvec,increments[j+1])
        c = [qs[j][x]+(1-qs[j][x])*vc[x] for x in range(size)]
        dvec = [vd[x]+qs[j][x]*vc[x] for x in range(size)]
    check(dot(nu,markov(c,increments[1]))/W==stop, "stop union recurrence")
    check(dot(nu,markov(dvec,increments[1]))/W==deficit, "deficit recurrence")
    pair_ind = F(0)
    for i in range(1,J+1):
        src = markov(nu,rho[i])
        for j in range(i+1,J+1):
            vq = markov(qs[j],rho[j]/rho[i])
            pair_ind += sum(src[x]*qs[i][x]*vq[x] for x in range(size))/W
    check(pair_ind==pair, "independent joint pair")
    # First-accept settlement against the genuinely future count.
    R = [None]*(J+1)
    R[J] = [F(0)]*size
    for j in range(J-1,0,-1):
        R[j] = markov(add(qs[j+1],R[j+1]),increments[j+1])
    alive = scale(nu,1/W)
    first_mass,future_def = F(0),F(0)
    for j in range(1,J+1):
        ell = markov(alive,increments[j])
        zeta = [ell[x]*qs[j][x] for x in range(size)]
        first_mass += sum(zeta)
        future_def += dot(zeta,R[j])
        alive = [ell[x]-zeta[x] for x in range(size)]
    check(first_mass==stop and first_mass<=1, "first accept source once")
    check(future_def==deficit, "first accept future count identity")
    check(sum(alive)==law[0], "surviving no accept")
    # Receiver occurrence intensity and Palm, without source=Lebesgue.
    for j in range(1,J+1):
        receiver = sensor(markov(markov(nu,rho[j]),rho[j]),j,d,varying)
        intensity = [gs[j][x]*receiver[x]/W for x in range(size)]
        check(intensity==[tau/W if winners.get(x)==j else F(0) for x in range(size)],
              "uniform accepted receiver intensity")
    if mean>0:
        palm = [F(k)*p/mean for k,p in enumerate(law)]
        reciprocal = sum(palm[k]/k for k in range(1,J+1))
        check(sum(palm)==1 and reciprocal==stop/mean, "Palm reciprocal")
        check(sum(palm[k]*F(k-1,k) for k in range(1,J+1))==deficit/mean,
              "Palm overlap")
        check(sum(F(k-1)*palm[k] for k in range(1,J+1))==2*pair/mean,
              "Palm pair")
        for K in (1,max(1,J//2),J):
            small = sum(palm[:K+1])
            check(reciprocal>=small/K, "Palm truncated sufficient bound")
    else:
        check(stop==deficit==pair==0 and E==[], "empty strict event convention")
    # Source conditional tail and exact signed saturation duality.
    total_g = sum(sum(gs[j]) for j in range(1,J+1))
    for m in range(1,J+1):
        cm = [sum(conditional[k][x] for k in range(m,J+1)) for x in range(size)]
        check(min(cm)>=0 and max(cm)<=1, "conditional tail range")
        check(sum(cm)<=total_g/m, "Lebesgue tail trace bound")
        check(dot(nu,cm)==dot(mu,cm)+dot(u,generator(cm,d)),
              "saturated signed tail identity")
        check(dot(nu,cm)/W==tails[m-1], "conditional tail source weighting")
    # Conditional whole-path law: only here the coin variables are independent.
    for path in ([0]*J,[(j*3)%size for j in range(J)]):
        poly = [F(1)]
        probs = [qs[j][path[j-1]] for j in range(1,J+1)]
        for p in probs:
            out = [F(0)]*(len(poly)+1)
            for k,val in enumerate(poly):
                out[k] += (1-p)*val
                out[k+1] += p*val
            poly = out
        check(sum(poly)==1, "conditional path coin normalization")
        check(sum(F(k)*p for k,p in enumerate(poly))==sum(probs),
              "conditional path mean")
        check(sum(F(comb(k,2))*p for k,p in enumerate(poly))==
              sum(probs[i]*probs[j] for i in range(J) for j in range(i+1,J)),
              "conditional path pair")
    # Exact frequentist CI calibration, not a Monte Carlo run.
    delta,intervals = inversions(8,J)
    for p in sorted(set(tails+[F(0),F(1,2),F(1)])):
        sampling = binlaw(8,p)
        coverage = sum(sampling[s] for s,(lo,hi) in enumerate(intervals) if lo<=p<=hi)
        check(coverage>=1-delta, "exact CI coverage")
    return {
        "sensor": "varying symmetric flip mixture" if varying else "identity",
        "tau_factor": strf(tau_factor), "kappa": strf(kappa), "W": strf(W),
        "receiver_count": len(E), "count_law": [strf(p) for p in law],
        "mean": strf(mean), "stop": strf(stop), "deficit": strf(deficit),
        "pair": strf(pair), "checks": CHECKS-before,
        "empty_event": not E
    }

def reflection(J):
    before = CHECKS
    size = 2*J+1
    source = [F(int(x==0)) for x in range(size)]
    responses,gs = [],[]
    for a in range(1,J+1):
        response = [source[(a-x)%size] for x in range(size)]
        g = [F(1,2) if x==a else F(0) for x in range(size)]
        responses.append(response); gs.append(g)
        check(response==[F(int(x==a)) for x in range(size)],"reflection response")
        check(all((a-(a-x))%size==x for x in range(size)), "selfadjoint reflection")
        check(g[a]==F(1,2) and response[a]>F(1,2), "strict truewinner weight")
    law = binlaw(J,F(1,2))
    mean = sum(F(k)*p for k,p in enumerate(law))
    reciprocal = (1-law[0])/mean
    check(mean==F(J,2),"reflection count mean")
    check(reciprocal==2*(1-F(1,2)**J)/J,"reflection Palm harmonic pressure")
    palm = [F(k)*p/mean for k,p in enumerate(law)]
    check(palm[1:]==binlaw(J-1,F(1,2)),"reflection size bias law")
    check(sum(dot(g,f) for g,f in zip(gs,responses))==mean,
          "reflection uniform receiver intensity")
    return {"J":J,"group_size":size,"count_law":[strf(p) for p in law],
            "mean":strf(mean),"Palm_reciprocal":strf(reciprocal),
            "checks":CHECKS-before,
            "scope":"Selfadjoint bistochastic sensor boundary, not centered h or actual obstacle."}

def main():
    start = time.monotonic()
    reg = BASE/(PREFIX+"_registration.json")
    registration = json.loads(reg.read_text())
    check(registration["registered_before_execution"],"registration exists")
    rounds = []
    for d,J in ((3,3),(4,6),(5,9)):
        before = CHECKS
        cases = [case(d,J,varying,tau) for varying in (False,True)
                 for tau in (F(1,16),F(1,8),F(1,4),F(1))]
        ref = reflection(J)
        rounds.append({"d":d,"J":J,"cases":cases,"reflection":ref,
                       "checks":CHECKS-before})
        print(json.dumps({"round":d,"J":J,"checks":CHECKS-before,
                          "empty_cases":sum(c["empty_event"] for c in cases)}),flush=True)
    result = {
        "status":"PASS", "checks":CHECKS, "random_seed":None,
        "scope":registration["scope"], "original_source_MC":"NOT_EXECUTED",
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "registration_sha256":hashlib.sha256(reg.read_bytes()).hexdigest(),
        "elapsed_seconds":time.monotonic()-start,"rounds":rounds
    }
    (BASE/(PREFIX+"_results.json")).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("status","checks","elapsed_seconds","script_sha256")}),flush=True)

if __name__=="__main__":
    main()
