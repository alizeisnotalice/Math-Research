"""Replay the rational margin and volume bounds for a finite example."""
from fractions import Fraction
from math import comb

m = 8
eps = Fraction(1, m*m)
width = 4 - 2*eps
source_base = sum(comb(m, j) * 2**j for j in range(m//2+1))
hit = Fraction(1, 2) + Fraction(comb(m, m//2), 2**(m+1))
margin = hit * ((1-eps)/(1+eps))**m - Fraction(1, 3)
source = source_base * (width-eps)**m
conflict_lower = 3**m * eps * width**(m-1)
assert margin > 0 and conflict_lower > 0
print("strict margin >=", margin)
print("single-axis conflict/source >=", conflict_lower/source)
