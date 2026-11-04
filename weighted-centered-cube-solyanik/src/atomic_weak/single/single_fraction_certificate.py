#!/usr/bin/env python3
"""Exact root-power certificate: J_{1/8}(17/16)>27/4>old B."""
from fractions import Fraction as Q
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
def main():
 lminus=Q(707,1000);lplus=Q(547,500);upperB=Q(13,12)
 checks={'707/1000 squared below 1/2':lminus*lminus<Q(1,2),'547/500 eighth below 33/16':lplus**8<Q(33,16),'13/12 eighth above 5/3':upperB**8>Q(5,3)}
 assert all(checks.values())
 lower=Q(64,17)*(lminus+lplus);assert lower==Q(14408,2125) and lower-Q(27,4)==Q(257,8500)
 bupper=6*upperB;assert bupper==Q(13,2) and Q(27,4)>bupper
 out={'tau':'1/8','r':'17/16','power_checks':checks,'J_lower':str(lower),'J_minus_27over4':str(lower-Q(27,4)),'B_upper':str(bupper),'conclusion':'J > 27/4 > 13/2 > B'}
 (ROOT/'single_fraction_certificate.json').write_text(json.dumps(out,indent=2));print('single_fraction_certificate PASS: J>14408/2125>27/4>13/2>B')
if __name__=='__main__':main()
