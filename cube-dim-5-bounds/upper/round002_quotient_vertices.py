#!/usr/bin/env python3
"""Exact complete vertex/facet enumeration for the four-dimensional quotient ball."""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations, product
from functools import reduce
from math import gcd
import hashlib, json
from pathlib import Path

def det(matrix):
    a=[[Fraction(x) for x in row] for row in matrix]; n=len(a); previous=Fraction(1); sign=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        d=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=a[i][j]*d-a[i][k]*a[k][j]
                a[i][j]=z/previous
        previous=d
        for i in range(k+1,n):a[i][k]=0
    return sign*a[-1][-1]

def rank(rows):
    for r in range(min(len(rows),4),0,-1):
        for selected in combinations(rows,r):
            for columns in combinations(range(4),r):
                if det(tuple(tuple(row[j] for j in columns) for row in selected)):return r
    return 0

def kernel(rows):
    cof=[]
    for j in range(4):cof.append(int((-1)**j*det(tuple(tuple(x for k,x in enumerate(row) if k!=j) for row in rows))))
    divisor=reduce(gcd,(abs(x) for x in cof if x),0); result=tuple(x//divisor for x in cof)
    return result if next(x for x in result if x)>0 else tuple(-x for x in result)

def affine_rank(points):
    return 0 if len(points)<2 else rank([tuple(x-y for x,y in zip(point,points[0])) for point in points[1:]])

def main():
    forms=tuple(tuple((mask>>j)&1 for j in range(4)) for mask in range(1,16))
    vertices=set(); independent=0
    for triple in combinations(forms,3):
        if rank(triple)!=3:continue
        independent+=1; ray=kernel(triple)
        psi=2*sum(abs(sum(a*b for a,b in zip(form,ray))) for form in forms)
        assert psi>0
        for sign in (1,-1):
            vertex=tuple(Fraction(sign*a,psi) for a in ray)
            zero_forms=[form for form in forms if sum(a*b for a,b in zip(form,vertex))==0]
            assert rank(zero_forms)==3
            vertices.add(vertex)
    support=sorted({tuple(sum(signs[i]*form[j] for i,form in enumerate(forms)) for j in range(4)) for signs in product((-1,1),repeat=15)})
    assert all(2*sum(a*b for a,b in zip(z,vertex))<=1 for z in support for vertex in vertices)
    facets=[]; histogram={}
    for z in support:
        equality=[v for v in vertices if 2*sum(a*b for a,b in zip(z,v))==1]
        key=(len(equality),affine_rank(equality));histogram[str(key)]=histogram.get(str(key),0)+1
        if affine_rank(equality)==3:facets.append((z,equality))
    payload={"quotient_symbol":"2 sum_{nonempty S subset [4]} |sum_{j in S} xi_j|","form_count":15,"independent_triple_count":independent,"complete_vertex_count":len(vertices),"support_vector_count":len(support),"exact_facet_count":len(facets),"facet_vertex_count_histogram":{str(k):sum(len(e)==k for _,e in facets) for k in sorted({len(e) for _,e in facets})},"support_equality_rank_histogram":histogram,"vertex_principle":"A vertex requires and is certified by zero directions of rank 3 plus the nonzero linear gradient; every rank-3 triple was enumerated exactly.","facet_principle":"For each exact support vector, equality vertices were checked by Fraction and a facet is precisely affine rank 3."}
    assert len(vertices)==90 and len(facets)==370
    raw=json.dumps(payload,sort_keys=True,indent=2).encode()+b"\n";Path(__file__).with_name("round002_quotient_vertices_results.json").write_bytes(raw)
    print(raw.decode(),end="");print("sha256",hashlib.sha256(raw).hexdigest())

if __name__=="__main__":main()
