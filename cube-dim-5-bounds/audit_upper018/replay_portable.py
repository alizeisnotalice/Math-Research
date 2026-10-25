#!/usr/bin/env python3
"""Rebuild locked audit inputs from existing project files without archiving copies.

The default action checks hashes and stages inputs in an automatically removed
directory under this audit folder.  Pass --run to execute the unchanged checker.
No source files or resources outside this audit folder are written.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_relative(name):
    p = Path(name)
    if p.is_absolute() or ".." in p.parts:
        raise ValueError(f"Not a safe relative path: {name}")
    return p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=HERE.parents[1],
                        help="Repository containing the already imported source files")
    parser.add_argument("--run", action="store_true", help="Also run the unchanged checker")
    args = parser.parse_args()
    registration = json.loads((HERE / "registration.json").read_text())
    expected = json.loads((HERE / "check_results.json").read_text())
    assert digest(HERE / "registration.json") == expected["registration_sha256"]
    assert digest(HERE / "check_interfaces.py") == expected["source_sha256"]
    entries = registration["inputs"] + registration.get("extra_inputs", [])
    with tempfile.TemporaryDirectory(prefix="portable-stage-", dir=HERE) as tmp:
        stage = Path(tmp)
        for name in ["check_interfaces.py", "registration.json"]:
            shutil.copy2(HERE / name, stage / name)
        for entry in entries:
            source = args.source_root.resolve() / safe_relative(entry["path"])
            assert digest(source) == entry["sha256"], f"Source hash mismatch: {source}"
            target = stage / safe_relative(entry["snapshot"])
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            assert digest(target) == entry["sha256"]
        if args.run:
            subprocess.run([sys.executable, "-B", str(stage / "check_interfaces.py")],
                           cwd=stage, check=True, stdout=subprocess.DEVNULL)
            assert (stage / "check_results.json").read_bytes() == (HERE / "check_results.json").read_bytes()
    print(json.dumps({"status": "pass", "locked_inputs_staged_and_hash_checked": len(entries),
                      "checker_run": args.run, "temporary_tree_removed": True,
                      "original_files_modified": False}, indent=2))


if __name__ == "__main__":
    main()
