#!/usr/bin/env python3
"""Validate frozen round011 files without importing or executing their script.

Default (stdlib only; no Gamma calculation and no files written):
  python round011_replay_portable.py --project-root /path/to/repository

The project root may also be /path/to/repository/cube-dim-5-bounds.
Explicit reproduction, requiring the original python-flint dependency:
  python round011_replay_portable.py --project-root /path/to/repository --run
  python round011_replay_portable.py --project-root /path/to/repository \
      --run --output-dir /path/to/new-empty-output-directory

--run copies the script and registration into a separate output directory,
rewrites only the three upstream provenance paths in that copy, and executes
the copy. Its unchanged `root = Path(__file__).parent` redirects all output
to the copy's directory. Frozen originals are never imported or overwritten.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


FROZEN = {
    "round011_shape_integral.md":
        "eda3dd844c83c8a333505f0023923fe5fccd1fc24edec9f81b98016dc6b6b96b",
    "round011_shape_integral.py":
        "ba1fcfe90e273515a7afb14ce64aef5200273cec480604dd6af29a65aa8efe9d",
    "round011_shape_integral_registration.json":
        "3b2ef0a65a29d4a00a89699d99a36abe5596f966f1975d31ef6ae64444fd397e",
    "round011_shape_integral_results.json":
        "5b77e6a098f314e8d3cac1fbbc1625ea846e6baeec562f839e9d9d298846211d",
}
UPSTREAM = {
    "report.md": "99d1ad30fa9b1df3f59bc88ae0dfa782c4928280efb3550ecfd11a0cee58a462",
    "verify.py": "53b7bba30a503c38fb9e516680f45d6e2d1de7dbcb7a2fee2a88a26212088d70",
    "verification.json": "d36533519b150137d9a87f0e85c3a8fc793627ea79874cfebc4c39f033cbfad5",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def project_directory(given: Path | None) -> Path:
    candidate = (given if given is not None else Path(__file__).resolve().parent.parent).resolve()
    if (candidate / "cube-dim-5-bounds").is_dir():
        candidate /= "cube-dim-5-bounds"
    require((candidate / "upper").is_dir(), f"No upper/ project directory under {candidate}")
    return candidate


def validate(project: Path) -> tuple[dict, str]:
    upper = project / "upper"
    frozen_hashes = {name: sha256(upper / name) for name in FROZEN}
    require(frozen_hashes == FROZEN, "A frozen round011 file has changed; refusing execution")
    saved = json.loads((upper / "round011_shape_integral_results.json").read_text())
    source_path = upper / "round011_shape_integral.py"
    source = source_path.read_text()
    tree = ast.parse(source, filename=str(source_path))
    root_assignments = [node for node in ast.walk(tree)
                        if isinstance(node, ast.Assign)
                        and any(isinstance(target, ast.Name) and target.id == "root"
                                for target in node.targets)]
    require(len(root_assignments) == 1, "Expected exactly one production output-root assignment")
    require(ast.unparse(root_assignments[0].value) == "Path(__file__).parent",
            "Unexpected production output-root expression")
    mappings = []
    rewritten = source
    for name, expected in UPSTREAM.items():
        candidates = [old for old, digest in saved["source_sha256"].items()
                      if old.endswith("/astra_upper001/" + name) and digest == expected]
        require(len(candidates) == 1, f"Missing locked upstream provenance record for {name}")
        old = candidates[0]
        mapped = project / "astra_upper001" / name
        require(mapped.is_file(), f"Missing mapped upstream file: {mapped}")
        require(sha256(mapped) == expected, f"Mapped upstream hash mismatch: {mapped}")
        literal_nodes = [node for node in ast.walk(tree)
                         if isinstance(node, ast.Constant) and node.value == old]
        require(len(literal_nodes) == 1, f"Expected one source literal for {old}")
        literal = ast.get_source_segment(source, literal_nodes[0])
        require(literal is not None and rewritten.count(literal) == 1,
                f"Ambiguous provenance path replacement for {old}")
        rewritten = rewritten.replace(literal, repr(str(mapped)), 1)
        mappings.append({"original": old, "mapped": str(mapped), "sha256": expected})
    ast.parse(rewritten, filename="<portable-round011-copy>")
    receipt = {
        "status": "validated_without_execution",
        "project_root": str(project),
        "frozen_sha256": frozen_hashes,
        "upstream_path_mappings": mappings,
        "production_AST_valid": True,
        "rewritten_copy_AST_valid": True,
        "production_imported": False,
        "production_executed": False,
        "Gamma_baseline_recomputed": False,
        "wrapper_sha256": sha256(Path(__file__).resolve()),
    }
    return receipt, rewritten


def replay(project: Path, receipt: dict, rewritten: str, requested: Path | None) -> None:
    if requested is None:
        output = Path(tempfile.mkdtemp(prefix="round011-portable-"))
    else:
        output = requested.expanduser().resolve()
        require(output != project and project not in output.parents,
                "Replay output must be outside the archived project")
        if output.exists():
            require(output.is_dir() and not any(output.iterdir()), "Replay output directory must be empty")
        else:
            output.mkdir(parents=True)
    require(output != project and project not in output.parents,
            "Replay output must be outside the archived project")
    upper = project / "upper"
    copied_script = output / "round011_shape_integral.py"
    copied_script.write_text(rewritten)
    registration_name = "round011_shape_integral_registration.json"
    (output / registration_name).write_bytes((upper / registration_name).read_bytes())
    receipt.update({"status": "explicit_replay_requested", "output_dir": str(output),
                    "temporary_script_sha256": sha256(copied_script),
                    "output_root_redirection": "Path(__file__).parent is the temporary copy directory"})
    completed = subprocess.run([sys.executable, "-B", str(copied_script)], cwd=output,
                               text=True, capture_output=True, timeout=60)
    receipt.update({"production_executed": True, "Gamma_baseline_recomputed": True,
                    "returncode": completed.returncode,
                    "replay_stdout": completed.stdout, "replay_stderr": completed.stderr})
    # Revalidate originals even if the optional child process reported failure.
    after, _ = validate(project)
    receipt["frozen_hashes_unchanged_after_replay"] = after["frozen_sha256"] == FROZEN
    if completed.returncode == 0:
        result = output / "round011_shape_integral_results.json"
        require(result.is_file(), "Replay did not produce its result JSON")
        receipt.update({"status": "explicit_replay_completed", "replay_result": str(result),
                        "replay_result_sha256": sha256(result)})
    else:
        receipt["status"] = "explicit_replay_failed"
    (output / "round011_portable_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    require(completed.returncode == 0, f"Replay failed; see {output}/round011_portable_receipt.json")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project-root", type=Path, help="Repository root or cube-dim-5-bounds directory")
    parser.add_argument("--run", action="store_true", help="Explicitly recompute only in a separate temporary copy")
    parser.add_argument("--output-dir", type=Path, help="New/empty replay directory outside the project; requires --run")
    args = parser.parse_args()
    if args.output_dir is not None and not args.run:
        parser.error("--output-dir requires explicit --run; default validation writes no files")
    try:
        project = project_directory(args.project_root)
        receipt, rewritten = validate(project)
        if args.run:
            replay(project, receipt, rewritten, args.output_dir)
    except (OSError, ValueError, SyntaxError, subprocess.SubprocessError) as exc:
        parser.exit(1, f"round011 portable validation failed: {exc}\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
