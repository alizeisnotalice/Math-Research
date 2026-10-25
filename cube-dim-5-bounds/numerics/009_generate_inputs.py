#!/usr/bin/env python3
"""Arb time inputs and exact four-region ordered-face partition for 009."""
from __future__ import annotations
import hashlib, json, math, struct
from fractions import Fraction
from itertools import product
from pathlib import Path
from flint import arb, ctx

P, PANELS, D, N = 9, 9, 24, 25
REGIONS = (
 ("interior", ((Fraction(55,64),Fraction(57,64)),(Fraction(39,64),Fraction(41,64)),(Fraction(23,64),Fraction(25,64)),(Fraction(7,64),Fraction(9,64)))),
 ("near_S2", ((Fraction(4079,4096),Fraction(4087,4096)),(Fraction(4055,4096),Fraction(4063,4096)),(Fraction(24,4096),Fraction(32,4096)),(Fraction(8,4096),Fraction(16,4096)))),
 ("boundary", ((Fraction(47,64),Fraction(49,64)),(Fraction(31,64),Fraction(33,64)),(Fraction(15,64),Fraction(17,64)),(Fraction(0),Fraction(1,2048)))),
 ("near_corner", ((Fraction(251,256),Fraction(253,256)),(Fraction(243,256),Fraction(245,256)),(Fraction(235,256),Fraction(237,256)),(Fraction(227,256),Fraction(229,256)))),
)
def bits(x): return f"0x{struct.unpack('>Q',struct.pack('>d',x))[0]:016x}"
def enc(x):
 lo,hi=x.lower(),x.upper(); a=float(x); b=float(x)
 while arb(a)>lo:a=math.nextafter(a,-math.inf)
 while arb(math.nextafter(a,math.inf))<=lo:a=math.nextafter(a,math.inf)
 while arb(b)<hi:b=math.nextafter(b,math.inf)
 while arb(math.nextafter(b,-math.inf))>=hi:b=math.nextafter(b,-math.inf)
 return bits(a),bits(b)
def fracenc(x): return enc(arb(x.numerator)/x.denominator)
def main():
 ctx.prec=256; pi=arb.pi(); rho=arb(49)/10; cp=arb(35509)/(arb(764411904)*pi**5); err=arb(0); rows=[]; smallest=None
 for panel in range(PANELS):
  a=arb(1)/(2**(P-panel)); sig=a*(arb(3)/2-(rho+1/rho)/4); err+=a*4*cp*sig**-5*rho**(-D)/(rho-1)
  for j in range(N):
   theta=arb(2*j+1)*pi/(2*N); t=a*(arb(3)/2+theta.cos()/2); w=a/N
   for k in range(1,D//2+1): w-=(a/N)*2*(2*k*theta).cos()/(4*k*k-1)
   assert w.lower()>0; smallest=w if smallest is None or w.lower()<smallest.lower() else smallest
   rows.append(f"{panel} {j} {' '.join(enc(t)+enc(w))}")
 Path("009_time_nodes.bits").write_text("\n".join(rows)+"\n",encoding="ascii")
 spatial=[]; descriptions=[]
 for rid,(name,bounds) in enumerate(REGIONS):
  assert all(bounds[k][0]>=bounds[k+1][1] for k in range(3)); cells=[]; vol=Fraction(0)
  for cid,sel in enumerate(product((0,1),repeat=4)):
   sub=[]
   for bit,(lo,hi) in zip(sel,bounds):
    mid=(lo+hi)/2; sub.append((lo,mid) if bit==0 else (mid,hi))
   assert all(sub[k][0]>=sub[k+1][1] for k in range(3)); v=math.prod((hi-lo for lo,hi in sub),start=Fraction(1));vol+=v
   parts=[str(rid),str(cid)]
   for lo,hi in ((Fraction(1),Fraction(1)),)+tuple(sub):parts.extend(fracenc(lo)+fracenc(hi))
   spatial.append(" ".join(parts));cells.append({"cell":cid,"volume":str(v)})
  descriptions.append({"id":rid,"name":name,"bounds":[[str(a),str(b)] for a,b in bounds],"volume":str(vol),"cells":cells})
 Path("009_space_boxes.bits").write_text("\n".join(spatial)+"\n",encoding="ascii")
 out={"status":"009_inputs_complete","precision_bits":256,"epsilon":"1/512","T":1,"panels":PANELS,"degree":D,"N":N,"nodes":len(rows),"rho":"49/10","all_weights_certified_positive":True,"minimum_weight_arb":str(smallest),"Cpeak_arb":str(cp),"interpolation_error_arb":str(err),"tail_error_arb":str(cp/4),"time_nodes_sha256":hashlib.sha256(Path("009_time_nodes.bits").read_bytes()).hexdigest(),"space_boxes_sha256":hashlib.sha256(Path("009_space_boxes.bits").read_bytes()).hexdigest(),"regions":descriptions,"generator_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 Path("009_inputs.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8");print(json.dumps({k:out[k] for k in ("nodes","time_nodes_sha256","space_boxes_sha256","interpolation_error_arb","tail_error_arb")},indent=2))
if __name__=="__main__":main()
