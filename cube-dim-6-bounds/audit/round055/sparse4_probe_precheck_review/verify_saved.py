"""Read-only B055 closure verification; no structural rerun or producer execution."""
import sys
sys.dont_write_bytecode = True
import hashlib
import json
import tarfile
from collections import Counter
from pathlib import Path

B = Path(__file__).resolve().parent
EXPECTED_C055_MANIFEST = "e5a61613fa0af7d235c5bd8017e8a9b1171b3577d06509c31ccb75f74088892c"
EXPECTED_C055_SOURCE = "d44781595a04e8be9a4292e724fa7f3b3d63703494c12a2790399c63859c5af2"

def read(path):
    return json.loads(Path(path).read_text())

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify():
    checks = []
    def pinned(group, path, expected):
        actual = sha(path)
        assert actual == expected, (group, str(path), expected, actual)
        checks.append((group, str(path)))
    cfg = read(B / "input_manifest.json")
    C = Path(cfg["producer"])
    receipt = read(B / "execution_receipt.json")
    result = read(B / "results.json")
    header = read(B / "header_binding.json")
    frozen = read(C / "precheck_manifest.json")
    build = read(C / "build_record.json")
    pinned("final_manifest", C / "precheck_manifest.json", EXPECTED_C055_MANIFEST)
    assert build["source_sha256"] == EXPECTED_C055_SOURCE
    for path, digest in cfg["sha256"].items():
        pinned("original_input_pin", path, digest)
    for path, digest in cfg["auditor_source_sha256"].items():
        pinned("auditor_pin", B / path, digest)
    pinned("receipt_input", B / "input_manifest.json", receipt["input_manifest_sha256"])
    pinned("receipt_result", B / "results.json", receipt["results_sha256"])
    for path, digest in frozen["whitelist_sha256"].items():
        pinned("final_whitelist", C / path, digest)
    for path, digest in frozen["external_readonly_sha256"].items():
        pinned("final_external", path, digest)
    for path, digest in header["headers_sha256"].items():
        pinned("compile_header", path, digest)
    assert header["archive_sha256"] == build["headers_sha256"]
    with tarfile.open(build["headers_archive"]) as archive:
        members = {}
        for member in archive.getmembers():
            if member.isfile():
                name = Path(member.name).name
                assert name not in members, ("ambiguous archive basename", name)
                members[name] = member
        for path, digest in header["headers_sha256"].items():
            actual = hashlib.sha256(archive.extractfile(members[Path(path).name]).read()).hexdigest()
            assert actual == digest, ("header_archive_content", path)
            checks.append(("header_archive_content", path))
    assert frozen["binary_sha256"] == build["binary_sha256"] == result["driver"]["binary_sha256"]
    assert result["driver"]["source_sha256"] == EXPECTED_C055_SOURCE
    assert result["status"] == receipt["status"] == "PASS_PRECHECK_STRUCTURE_AND_SOURCE_ONLY"
    assert receipt["exit_code"] == 0 and receipt["all_handles_closed"]
    assert receipt["configurations"] == result["configurations"] == 6
    assert receipt["matrix_operations"] == result["matrix_operations"] == 22952
    assert 0 <= receipt["actual_seconds"] <= 900
    assert receipt["production_native_calls"] == result["production_G_H_p_native"] == 0
    assert receipt["author_writer_runs"] == result["author_writer_runs"] == 0
    assert receipt["old_formal_repeats"] == result["old_formal_repeats"] == 0
    assert result["production_launch_authorized_by_this_audit"] is False
    summary = {
        "status": "PASS_SAVED_EVIDENCE_AND_FINAL_SOURCE_BINDING",
        "hash_checks": len(checks),
        "hash_check_groups": dict(sorted(Counter(g for g, _ in checks).items())),
        "producer_precheck_manifest_sha256": EXPECTED_C055_MANIFEST,
        "producer_source_sha256": EXPECTED_C055_SOURCE,
        "existing_configurations_reexecuted": 0,
        "new_scientific_configurations": 0,
        "new_matrix_operations": 0,
        "native_calls": 0,
        "author_writer_runs": 0,
    }
    manifest_path = B / "manifest.json"
    if manifest_path.exists():
        manifest = read(manifest_path)
        for path, digest in manifest["whitelist_sha256"].items():
            assert sha(B / path) == digest, ("sealed audit file", path)
        summary["sealed_audit_files_verified"] = len(manifest["whitelist_sha256"])
    return summary

if __name__ == "__main__":
    assert sys.argv[1:] == ["--verify"]
    print(json.dumps(verify(), indent=2))
