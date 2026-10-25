#!/usr/bin/env python3
"""Exact H-to-V, facet, and circuit-triangulation audit for quotient volume."""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations, product
import hashlib, json
from pathlib import Path
from round002_quotient_vertices import affine_rank, det, kernel, rank

def affine_circuit(ids, vertices):
    """The unique affine dependence of five 3D-affinely-spanning vertices."""
    columns=[(Fraction(1),)+vertices[i] for i in ids]
    rows=tuple(tuple(columns[col][row] for col in range(5)) for row in range(5))
    for chosen_rows in combinations(range(5),4):
        block=tuple(rows[row] for row in chosen_rows)
        coefficients=tuple((-1)**j*det(tuple(tuple(x for k,x in enumerate(row) if k!=j) for row in block)) for j in range(5))
        if any(coefficients):
            assert sum(coefficients)==0
            assert all(sum(coefficients[col]*columns[col][row] for col in range(5))==0 for row in range(5))
            return coefficients
    raise RuntimeError("facet has no affine circuit")

def main():
    forms=tuple(tuple((mask>>j)&1 for j in range(4)) for mask in range(1,16))
    vertices=set(); independent=0
    for triple in combinations(forms,3):
        if rank(triple)!=3:continue
        independent+=1; ray=kernel(triple)
        psi=2*sum(abs(sum(a*b for a,b in zip(form,ray))) for form in forms)
        vertices.add(tuple(Fraction(a,psi) for a in ray));vertices.add(tuple(Fraction(-a,psi) for a in ray))
    vertices=sorted(vertices)
    support=sorted({tuple(sum(signs[i]*form[j] for i,form in enumerate(forms)) for j in range(4)) for signs in product((-1,1),repeat=15)})
    facets=[]
    for z in support:
        equality=frozenset(i for i,v in enumerate(vertices) if 2*sum(a*b for a,b in zip(z,v))==1)
        if affine_rank([vertices[i] for i in equality])==3:facets.append(tuple(sorted(equality)))
    simplices=[]; facet_counts=[]; circuit_hist={}
    for facet in facets:
        if len(facet)==4:
            local=[facet]
        else:
            assert len(facet)==5
            coefficients=affine_circuit(facet,vertices)
            positive=[i for i,value in enumerate(coefficients) if value>0]
            negative=[i for i,value in enumerate(coefficients) if value<0]
            assert {len(positive),len(negative)}=={2,3}
            chosen=positive if min(positive)<min(negative) else negative
            local=[tuple(facet[k] for k in range(5) if k!=omitted) for omitted in chosen]
            circuit_key=str((len(positive),len(negative)))
            circuit_hist[circuit_key]=circuit_hist.get(circuit_key,0)+1
        assert all(abs(det([vertices[i] for i in simplex]))>0 for simplex in local)
        simplices.extend(local);facet_counts.append(len(local))
    volume=sum(abs(det([vertices[i] for i in simplex])) for simplex in simplices)/24
    assert len(vertices)==90 and len(facets)==370 and len(simplices)==402
    assert volume==Fraction(14455,191102976)
    payload={"status":"exact_vertex_facet_circuit_triangulation_audit","rank3_triple_count":independent,"exact_vertex_count":len(vertices),"exact_facet_count":len(facets),"exact_simplex_count":len(simplices),"facet_to_simplex_count_histogram":{str(k):facet_counts.count(k) for k in sorted(set(facet_counts))},"five_vertex_facet_circuit_sign_histogram":circuit_hist,"triangulation_principle":"For a five-vertex affine circuit with Radon parts of sizes 2 and 3, deleting each vertex in either chosen part gives a disjoint tetrahedral triangulation of its convex hull. The deterministic choice uses the part with smaller least local index.","volume_from_exact_determinants":str(volume),"Cperp":"14455/(127401984*pi^4)","Cface":"14455/(31104*pi^5)"}
    raw=json.dumps(payload,sort_keys=True,indent=2).encode()+b"\n";Path(__file__).with_name("round002_quotient_volume_audit_results.json").write_bytes(raw)
    print(raw.decode(),end="");print("sha256",hashlib.sha256(raw).hexdigest())

if __name__=="__main__":main()
