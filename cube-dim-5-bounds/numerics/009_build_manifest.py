#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json
from pathlib import Path
def h(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 files=['009_generate_inputs.py','009_time_nodes.bits','009_space_boxes.bits','009_inputs.json','009_face_box_heat.cpp','009_region_0_heat.json','009_finalize_diagnostic.py','009_diagnostic.json','009_registration.md']
 out={'status':'finite_009_manifest','binary_path':'/tmp/cube5-009_face_box_heat-macos-arm64','binary_sha256':h('/tmp/cube5-009_face_box_heat-macos-arm64'),'files_sha256':{x:h(x) for x in files},'scope':'Only region0 ran; regions1-3 have no heat outputs by declared stop condition.'}
 Path('009_build_manifest.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
