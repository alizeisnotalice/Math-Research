#!/usr/bin/env python3
"""Portable verifier for the frozen audit_lower019 files.

Default: python3 -B portable.py --project-root /checkout/cube-dim-5-bounds
Only --run replays the bounded checker, in a new retained temporary directory.
The project root may also be the repository containing cube-dim-5-bounds/.
Original absolute source paths are never opened or used as fallbacks.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
FROZEN = {
    "audit.md": "0c7e11dce9cadcd043375782a1edc21a9e44c783d669f740d84b4b9732bc640d",
    "check_saved_coefficients.py": "a10fa1c72946171018e5f8d9c803fcefe724860d2c39fdb915e939b02439d305",
    "registration.json": "72f495a4331ea7d7f493c0617bb3351fea37aaffc27181a7e65122761e6cd1fe",
    "results.json": "a580bb8a1018a692d4b3a90f3be027f7b0459084fb9f3b147ff3925e78f79fdd",
    "completion_receipt.json": "a351ba64b58fb8170247dadbd4d735f6bcb9e158bfe76851ff4720486c659b00",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mapped_path(original, project):
    parts = Path(original).parts
    if "cube-dim-5-bounds" in parts:
        suffix = parts[parts.index("cube-dim-5-bounds")+1:]
        if not suffix or suffix[0] not in {"lower_alt014", "lower", "lower_round003"}:
            raise ValueError("unexpected registered project source: " + original)
        mapped = project.joinpath(*suffix).resolve()
        if not mapped.is_relative_to(project):
            raise ValueError("registered source escapes project root")
        return mapped
    if Path(original).name == "AGENTS.md":
        return project.parent / "AGENTS.md"
    raise ValueError("unmapped registered source: " + original)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--run", action="store_true", help="replay only in a new temporary directory")
    args = parser.parse_args()
    project = args.project_root.resolve()
    if (project/"cube-dim-5-bounds").is_dir():
        project = project/"cube-dim-5-bounds"
    for name, expected in FROZEN.items():
        if sha(HERE/name) != expected:
            raise ValueError("frozen audit hash mismatch: " + name)
    registration = json.loads((HERE/"registration.json").read_text())
    if registration["script_sha256"] != FROZEN["check_saved_coefficients.py"]:
        raise ValueError("registered checker hash mismatch")
    registration["author_directory"] = str(project/"lower_alt014")
    registration["standing_C_path"] = str(project/"lower_round003/a2_equals_a3_results.json")
    for source in registration["locked_sources"]:
        path = mapped_path(source["path"], project)
        if sha(path) != source["sha256"]:
            raise ValueError("mapped source hash mismatch: " + str(path))
        source["path"] = str(path)
    # A non-main import creates definitions only; disable bytecode side effects.
    spec = importlib.util.spec_from_file_location("audit_lower019_frozen_checker", HERE/"check_saved_coefficients.py")
    module = importlib.util.module_from_spec(spec)
    previous = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    summary = {"status": "PORTABLE_HASH_AND_SAFE_IMPORT_PASS", "run": False,
               "project_root": str(project), "frozen_files": len(FROZEN),
               "mapped_sources": len(registration["locked_sources"]),
               "original_absolute_paths_used": False}
    if args.run:
        output = Path(tempfile.mkdtemp(prefix="audit_lower019_replay_"))
        checker = output/"check_saved_coefficients.py"
        shutil.copyfile(HERE/checker.name, checker)
        if sha(checker) != registration["script_sha256"]:
            raise ValueError("temporary checker hash mismatch")
        (output/"registration.json").write_text(json.dumps(registration, indent=2)+"\n")
        run = subprocess.run([sys.executable, "-B", str(checker)], cwd=output,
                             check=True, capture_output=True, text=True)
        summary.update(status="PORTABLE_TEMPORARY_REPLAY_PASS", run=True,
                       output_directory=str(output), checker_stdout=run.stdout.strip())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
