#!/usr/bin/env python3
"""Exact endpoint integrals and Fraction logarithm sign checks; no grid."""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import importlib.util
import json
from math import comb
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent


def locate_core(explicit):
    if explicit is not None:
        return explicit.resolve()
    choices = [
        ROOT.parent / "round008" / "certify.py",
        ROOT.parent.parent / "Math-Research-lowdim-agent" / "cube-lowdim-bounds"
        / "experiments" / "round008" / "certify.py",
    ]
    for path in choices:
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError("Pass --core <round008/certify.py>.")


def endpoint_coefficients(p):
    # h_p(1)=2/(p-1)+2*integral_1^2 (2-z)^p/z^2 dz.
    rational = Q(2, p-1)
    logarithmic = Q(0)
    for k in range(p+1):
        coefficient = Q(comb(p,k)*2**(p-k)*(-1)**k)
        if k == 1:
            logarithmic += 2*coefficient
        else:
            primitive_difference = (Q(2)**(k-1)-1)/(k-1)
            rational += 2*coefficient*primitive_difference
    return rational, logarithmic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--core", type=Path)
    arguments = parser.parse_args()
    spec = importlib.util.spec_from_file_location("outer_constant_fraction_core",
                                                locate_core(arguments.core))
    core = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = core
    spec.loader.exec_module(core)
    expected = {2:(Q(8),Q(-8)), 3:(Q(18),Q(-24)), 4:(Q(136,3),Q(-64))}
    coefficients = {}
    bounds = {}
    ln2 = core.log_point(2)
    for p in (2,3,4):
        pair = endpoint_coefficients(p)
        if pair != expected[p]:
            raise ValueError("Endpoint integral identity mismatch")
        coefficients[p] = pair
        bounds[p] = pair[0]+pair[1]*ln2
        if bounds[p].lo <= 0:
            raise ValueError("Positivity not certified")
    differences = {}
    for p,q in ((2,3),(2,4),(3,4)):
        a,b = coefficients[p]
        c,d = coefficients[q]
        gap = (a-c)+(b-d)*ln2
        if gap.lo <= 0:
            raise ValueError("Strict ordering not certified")
        differences[f"h{p}_minus_h{q}"] = {
            "interval":gap.pair(),"lower_display":float(gap.lo)}
    print(json.dumps({
        "status":"exact_endpoint_integrals_and_strict_signs_verified",
        "method":"finite exact polynomial integration; reviewed Fraction log enclosure",
        "endpoint_coefficients":{str(p):[str(z) for z in coefficients[p]] for p in (2,3,4)},
        "endpoint_bounds":{str(p):{"interval":bounds[p].pair(),"display":float(bounds[p].lo)}
                           for p in (2,3,4)},
        "strict_differences":differences,
        "ordering":"0<h4(1)<h3(1)<h2(1)",
        "scope":"constant identities/signs only; no angular or radial sampling"
    },indent=2))


if __name__ == "__main__":
    main()
