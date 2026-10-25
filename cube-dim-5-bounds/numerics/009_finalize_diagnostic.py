#!/usr/bin/env python3
"""Quantify the first 009 box job against the strict H-integral budget."""
from __future__ import annotations
import hashlib,json,math,struct
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from body_fan_arb_anchor import exact_terms
from body_fan_h_epsilon_arb import truncated_green
def fbits(s):return struct.unpack('>d',int(s,16).to_bytes(8,'big'))[0]
def main():
 ctx.prec=256; raw=json.loads(Path('009_region_0_heat.json').read_text());inp=json.loads(Path('009_inputs.json').read_text())
 assert raw['time_nodes_sha256']==inp['time_nodes_sha256']==hashlib.sha256(Path('009_time_nodes.bits').read_bytes()).hexdigest()
 rho=arb(49)/10;cp=arb(35509)/(arb(764411904)*arb.pi()**5);err=arb(0)
 for k in range(9):
  a=arb(1)/(2**(9-k));sig=a*(arb(3)/2-(rho+1/rho)/4);err+=a*4*cp*sig**-5*rho**(-24)/(rho-1)
 tail=cp/4;q=arb(5)/4;region=inp['regions'][0];bycell={x['cell']:Fraction(x['volume']) for x in region['cells']};upper=arb(0);maxup=0.0
 for cell in raw['cells']:
  u=fbits(cell['quadrature']['hi_bits'])+float(err.upper())+float(tail.upper());maxup=max(maxup,u);upper+=arb(bycell[cell['cell']].numerator)/bycell[cell['cell']].denominator*arb(u)**q
 vol=Fraction(region['volume']); density=upper/(arb(vol.numerator)/vol.denominator)
 # Global short-time L^(5/4) norm budget supplied by the upper-theory interface.
 eps=arb(1)/512
 B=arb(16)**(-arb(4)/5)*(arb(5)/6)*(arb(256)/(3*arb.pi()**5))**(arb(1)/5)*(arb(16)/arb.pi())**(arb(4)/5)*eps**(arb(6)/5)
 hnorm=upper**(arb(4)/5); combined=(hnorm+B)**q
 center=(Fraction(1),Fraction(7,8),Fraction(5,8),Fraction(3,8),Fraction(1,8))
 direct=truncated_green(center,Fraction(1,512),tuple(exact_terms()))[0].real
 threshold=arb('3.96333458634898456e-5')
 out={'status':'direct_box_method_stopped_after_first_region','precision_bits':256,'executed_region':region['name'],'executed_volume':str(vol),'executed_cells':16,'unexecuted_preregistered_regions':[x['name'] for x in inp['regions'][1:]],'epsilon':'1/512','interpolation_error_arb':str(err),'tail_error_arb':str(tail),'max_cell_H_upper_float':maxup,'H_integral_upper_for_executed_region':str(upper),'H_integral_upper_per_volume':str(density),'H_integral_budget_arb':str(threshold),'budget_fraction_for_executed_region':str(upper/threshold),'short_time_global_L5_4_budget':str(B),'H_L5_4_upper_for_executed_region':str(hnorm),'Minkowski_G_integral_upper_for_executed_region':str(combined),'direct_center_H_1_512_arb_ball':str(direct),'cancellation_stop_reason':'At 16 subboxes, whole-box Fourier-cone interval cancellation yields H upper endpoints orders of magnitude above the direct center scale; further pre-registered regions would repeat the same non-covering mechanism.','input_hashes':{'time':inp['time_nodes_sha256'],'space':inp['space_boxes_sha256'],'heat_output':hashlib.sha256(Path('009_region_0_heat.json').read_bytes()).hexdigest()},'finalizer_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Finite feasibility diagnosis only; no full-face I certificate or C5 bound.'}
 Path('009_diagnostic.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
