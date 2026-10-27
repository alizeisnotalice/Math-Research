"""Verify frozen n=7 R006 manifests and exact union interval correspondence."""
import hashlib, json
from pathlib import Path

root = Path(__file__).resolve().parents[3]
n7 = root / "cube-dim-7-bounds"
lower = n7 / "lower/round006"
audit = n7 / "audit/round006/lower"

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def manifest_check(folder, manifest_name):
    manifest_path = folder / manifest_name
    manifest = json.loads(manifest_path.read_text())
    entries = manifest["files"]
    observed = {name: sha(folder / name) for name in entries}
    return {"manifest_sha256": sha(manifest_path), "files": {
        name: {"expected": expected, "actual": observed[name], "match": observed[name] == expected}
        for name, expected in entries.items()
    }}

lower_manifest = manifest_check(lower, "manifest.json")
audit_manifest = manifest_check(audit, "manifest.json")
orig_results = json.loads((lower / "results.json").read_text())
audit_results = json.loads((audit / "exactresults.json").read_text())
audit_registration = json.loads((audit / "registration.json").read_text())
interval_match = orig_results["full_construction_comparison_interval"] == [audit_results["union_lower"], audit_results["union_upper"]]
assert all(x["match"] for x in lower_manifest["files"].values())
assert all(x["match"] for x in audit_manifest["files"].values())
assert audit_registration["candidate_manifest_sha256"] == lower_manifest["manifest_sha256"]
assert audit_registration["candidate_results_sha256"] == sha(lower / "results.json")
assert interval_match
print(f"lower_manifest_sha256={lower_manifest['manifest_sha256']}")
print(f"lower_manifest_files={len(lower_manifest['files'])}; all_match=True")
print(f"audit_manifest_sha256={audit_manifest['manifest_sha256']}")
print(f"audit_manifest_files={len(audit_manifest['files'])}; all_match=True")
print(f"audit_exactresults_sha256={sha(audit / 'exactresults.json')}")
print(f"audit_registration_sha256={sha(audit / 'registration.json')}")
print("candidate_manifest_binding=True")
print("candidate_results_binding=True")
print("complete_8192_label_union_interval_match=True")
print("PASS: frozen source hashes and full rational interval correspondence verified")
