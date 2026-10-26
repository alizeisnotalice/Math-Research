"""Read-only precheck verification, not a full launch or baseline replay."""
import sys
sys.dont_write_bytecode=True
import json,hashlib
from pathlib import Path
from fractions import Fraction as Q
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 m=json.loads((P/'precheck_manifest.json').read_text())
 for f,h in m['whitelist_sha256'].items():assert sha(P/f)==h,('local',f)
 for f,h in m['external_readonly_sha256'].items():assert sha(Path(f))==h,('external',f)
 s=json.loads((P/'S2_certificate.json').read_text());assert s['actual_new_preflight'] and not s['C053_preflight_inherited']
 assert s['status']=='PASS_STRICT_OFF_TRUE_S2' and s['direction_count']==152 and s['pair_count']==11476 and Q(s['pair_min_squared'])==Q(363,4375)>0 and not s['zero_pairs']
 assert Q(s['line_min_squared'])==Q(2991,10000)>0
 tests=json.loads((P/'structure_checks.json').read_text());assert len(tests['tests'])==6<=8 and all(x['PASS'] for x in tests['tests'])
 data=json.loads((P/'input_certificate.json').read_text());assert data['input_sha256']==sha(P/'native_input.txt')
 assert data['paired_terms']==25399680 and data['logical_terms_each_DD_field']==609592320 and data['joint_DD_field_evaluations']==50799360
 assert json.loads((P/'native_validate.json').read_text())['positive_time_evaluations']==0 and json.loads((P/'abi.json').read_text())['full_calls']==0
 build=json.loads((P/'build_record.json').read_text());assert build['source_sha256']==sha(P/'sparse4_probe.cpp')
 dep=json.loads((P/'dependency_pins.json').read_text());assert dep['status']=='F050_AND_F052_FINAL_PASS_PINNED'
 print('PASS saved precheck: hashes,true152/11476certificate,6 structural receipts,input/core/ABI/dependencies;0 kernel or baseline replay/0 writes')
if __name__=='__main__':main()
