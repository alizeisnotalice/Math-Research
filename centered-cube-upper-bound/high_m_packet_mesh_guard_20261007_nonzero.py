#!/usr/bin/env python3
from fractions import Fraction as F
import json,pathlib,hashlib
P=pathlib.Path(__file__).resolve().parent
PRE="high_m_packet_mesh_guard_20261007_nonzero"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+"\n")
def main():
    regp=P/(PRE+"_registration.json");reg=json.loads(regp.read_text())
    assert sha(pathlib.Path(__file__))==reg["script_sha256"]
    checks=1;rows=[]
    for n in reg["n"]:
        for xs in reg["xi"]:
            xi=F(xs);d=xi/(2*n);R=F(7,5)
            z=(R-1)/d;j=-(-z.numerator//z.denominator);r=1+j*d
            Mt=R**(-n);Mm=r**(-n)
            tests=[1<R<2,1<=r<=2,r>R,r-d<R, j>=z and j-1<z,
                   (r/(r-d))**n<=1+xi, Mm<Mt, Mt<=(1+xi)*Mm,
                   (Mt/Mm)==(r/R)**n, Mt-Mm>0]
            assert all(tests);checks+=len(tests)
            rows.append({"n":n,"xi":xs,"source":{"atoms":[["0"]*n],"weights":["1"]},
                         "receiver":["7/10"]+["0"]*(n-1),"a":"1","b":"2","delta":str(d),
                         "true_R":str(R),"upward_mesh_index":j,"mesh_R":str(r),
                         "true_M":str(Mt),"mesh_M":str(Mm),"positive_error":str(Mt-Mm),
                         "multiplicative_ratio":str(Mt/Mm),"checks":len(tests)})
    rp=P/(PRE+"_results.json")
    dump(rp,{"status":"PASS","checks":checks,"rows":rows,"script_sha256":sha(pathlib.Path(__file__)),
             "registration_sha256":sha(regp),"oracle_invocations":0,
             "scope":"Original hard atom input with analytically unique continuous winner R=7/5; exact nonzero upward-mesh error. Not nearmax, general continuous numerical proof, or actual FIRST."})
    dump(P/(PRE+"_receipt.json"),{"status":"PASS","successful_execution_count":1,"failed_execution_count":0,
                                 "checks":checks,"script_sha256":sha(pathlib.Path(__file__)),"registration_sha256":sha(regp),"results_sha256":sha(rp)})
    print(json.dumps({"status":"PASS","records":len(rows),"checks":checks,"nonzero_error_records":len(rows)}))
if __name__=="__main__":main()

