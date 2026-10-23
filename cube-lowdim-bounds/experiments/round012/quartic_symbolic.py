"""Round 12 quartic angular formula and finite diagnostics.

The complete SymPy differentiation was executed once during derivation.
--derive reproduces that full symbolic check; default only verifies the
remaining rational identities and runs one fixed finite diagnostic group.
No LP, search, or file writes.
"""
from __future__ import annotations
import argparse
import json
import sympy as S
import mpmath as mp


def symbolic(derive=False):
    r, d, x = S.symbols("r d x", real=True)
    A = 25*r**3-42*r**2
    B = -46*r**3+84*r**2-12*r+8
    C = 5*r**3+6*r**2-36*r+8
    Q = -(A*x**2+B*x+C)/(2*r**2*(1-x)**2)
    expected = 2*(1-r)**2*((9*r-6)-(r+2)*x)/(r**2*(1-x)**3)
    assert S.cancel(S.diff(Q, x)-expected) == 0
    q0 = -S.Rational(5, 2)*r-3+18/r-4/r**2
    qstar = (165-101*r)/8+1/(2*r**2)
    assert S.cancel(Q.subs(x, 0)-q0) == 0
    assert S.cancel(Q.subs(x, (9*r-6)/(r+2))-qstar) == 0
    assert S.cancel((A+B+C)-16*(1-r)**3) == 0
    assert S.simplify((q0-qstar).subs(r, S.Rational(2, 3))) == 0
    d0 = 12*S.log(r/2)+q0
    assert S.cancel(S.diff(d0, r)+(r-2)**2*(5*r-4)/(2*r**3)) == 0
    out = {"rational_Q": str(Q), "Q_x": str(expected),
           "maximum_low": str(q0), "maximum_high": str(qstar),
           "rational_identities_verified": True}
    if derive:
        P = (88-132*r+78*r**2-17*r**3+(78-27*r)*d**2)/12
        Lp = S.log((r+d)/(2*(1+d)))
        Lm = S.log((r-d)/(2*(1-d)))
        B4 = (P+(1+d)**3*Lp+(1-d)**3*Lm
              +(1-r)**3*S.log((r**2-d**2)/(4*(1-r)**2)))
        expected_dd = 6*((1+d)*Lp+(1-d)*Lm)+Q.subs(x, d**2/r**2)
        assert S.cancel(S.diff(B4, d, 2)-expected_dd) == 0
        out["full_second_derivative_verified"] = True
    return out


def diagnostics():
    mp.mp.dps = 65
    def h(t):
        return 2*t*t+mp.mpf(11)/3+(1+t)**3*mp.log(t/(1+t))+(1-t)**3*mp.log(t/(1-t))
    def direct_B(r, d):
        u=(r+d)/2; v=(r-d)/2; a=1-v; b=1-u
        return a**3*h(u/a)+b**3*h(v/b)
    def closed_B(r, d):
        P=(88-132*r+78*r*r-17*r**3+(78-27*r)*d*d)/12
        return (P+(1+d)**3*mp.log((r+d)/(2*(1+d)))
                +(1-d)**3*mp.log((r-d)/(2*(1-d)))
                +(1-r)**3*mp.log((r*r-d*d)/(4*(1-r)**2)))
    def Q(r, x):
        A=25*r**3-42*r*r; B=-46*r**3+84*r*r-12*r+8
        C=5*r**3+6*r*r-36*r+8
        return -(A*x*x+B*x+C)/(2*r*r*(1-x)**2)
    def dd(r, d):
        return (6*((1+d)*mp.log((r+d)/(2*(1+d)))
                   +(1-d)*mp.log((r-d)/(2*(1-d))))+Q(r,d*d/(r*r)))
    def cap(r):
        if r<=mp.mpf(2)/3:
            q=-mp.mpf(5)*r/2-3+18/r-4/(r*r)
        else:
            q=(165-101*r)/8+1/(2*r*r)
        return 12*mp.log(r/2)+q
    radii=[mp.mpf(p)/q for p,q in [(1,5),(1,2),(3,5),(2,3),(3,4),(9,10),(99,100),(999,1000)]]
    xs=[mp.mpf(0),mp.mpf(1)/4,mp.mpf(3)/4,mp.mpf(99)/100]
    closed_error=mp.mpf(0);dd_error=mp.mpf(0)
    min_curvature=mp.inf;min_taylor=mp.inf;count=0
    for r in radii:
        for x in xs:
            d=r*mp.sqrt(x)
            closed_error=max(closed_error,abs(direct_B(r,d)-closed_B(r,d)))
            actual_dd=mp.diff(lambda z:direct_B(r,z),d,2)
            dd_error=max(dd_error,abs(actual_dd-dd(r,d)))
            min_curvature=min(min_curvature,cap(r)-actual_dd)
            min_taylor=min(min_taylor,4*cap(r)*d*d-8*(direct_B(r,d)-direct_B(r,0)))
            count+=1
    tol=mp.mpf("1e-45")
    assert closed_error<tol and dd_error<tol
    assert min_curvature>=-tol and min_taylor>=-tol
    return {"points":count,"mpmath_digits":mp.mp.dps,
            "max_slice_vs_closed_B4_error":str(closed_error),
            "max_numeric_differentiation_vs_B4dd_error":str(dd_error),
            "minimum_curvature_bound_slack":str(min_curvature),
            "minimum_Taylor_bound_slack":str(min_taylor),
            "all_finite_checks_passed":True,
            "scope":"finite diagnostic only; analytic proof is in quartic_angular.md"}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--derive",action="store_true")
    args=parser.parse_args()
    print(json.dumps({"symbolic":symbolic(args.derive),"diagnostics":diagnostics()},indent=2))
