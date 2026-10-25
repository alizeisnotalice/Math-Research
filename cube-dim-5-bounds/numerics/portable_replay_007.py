#!/usr/bin/env python3
"""Clean, finite replay of the locked 007 single-point certificate.

This launcher deliberately does not copy body_fan_check.py: its temporary
Arb helpers are rewritten to import a clean geometry-only module instead.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


REPLAY_SOURCES = (
    "005b_interval_heat.cpp", "007_interval_fejer.cpp",
    "007_generate_fejer_inputs.py", "007_final_verify.py",
)
DEPENDENCIES = (
    "body_fan_geometry.py", "body_fan_arb_anchor.py",
    "body_fan_h_epsilon_arb.py", "body_fan_green_arb.py",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], cwd: Path, timeout: float, log: list[dict]) -> subprocess.CompletedProcess[str]:
    began = time.monotonic()
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True, timeout=timeout, check=False)
    log.append({"command": command, "returncode": completed.returncode,
                "elapsed_seconds": time.monotonic() - began,
                "stdout": completed.stdout[-4000:], "stderr": completed.stderr[-4000:]})
    if completed.returncode:
        raise RuntimeError(f"command failed: {command[0]} (exit {completed.returncode})")
    return completed


def install_dependency(source: Path, destination: Path) -> str:
    text = source.read_text(encoding="utf-8")
    if source.name == "body_fan_arb_anchor.py":
        old = "from body_fan_check import BASE_CONES, directions"
        new = "from body_fan_geometry import BASE_CONES, directions"
        if old in text:
            text = text.replace(old, new)
        elif "body_fan_geometry" not in text:
            raise RuntimeError("anchor has neither a known legacy nor clean geometry import")
    if source.name == "body_fan_green_arb.py":
        old = "from body_fan_check import directions"
        new = "from body_fan_geometry import directions"
        if old in text:
            text = text.replace(old, new)
        elif "body_fan_geometry" not in text:
            raise RuntimeError("Green helper has neither a known legacy nor clean geometry import")
    destination.write_text(text, encoding="utf-8")
    return digest(destination)


def write_receipt(path: Path, receipt: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True, help="new receipt JSON path")
    parser.add_argument("--dependency-dir", type=Path, default=Path(__file__).resolve().parent,
                        help="clean geometry + Arb helper directory; defaults to script directory")
    parser.add_argument("--force", action="store_true", help="allow replacing an existing receipt")
    parser.add_argument("--timeout", type=float, default=120.0)
    args = parser.parse_args()
    source_dir = Path(__file__).resolve().parent
    output = args.output.resolve()
    if output.exists() and not args.force:
        raise SystemExit(f"refusing to overwrite existing output: {output}")
    receipt: dict = {
        "status": "in_progress", "replay": "007_single_point_clean_portable_replay",
        "source_dir": str(source_dir), "dependency_dir": str(args.dependency_dir.resolve()),
        "python": sys.version, "platform": platform.platform(), "commands": [],
        "forbidden_not_copied": ["body_fan_check.py"],
    }
    began = time.monotonic()
    try:
        needed = [source_dir / item for item in REPLAY_SOURCES] + [args.dependency_dir / item for item in DEPENDENCIES]
        missing = [str(path) for path in needed if not path.is_file()]
        if missing:
            raise FileNotFoundError("missing replay dependency: " + ", ".join(missing))
        receipt["original_sha256"] = {str(path): digest(path) for path in needed}
        with tempfile.TemporaryDirectory(prefix="cube5-replay-007-") as temporary:
            work = Path(temporary)
            for name in REPLAY_SOURCES:
                shutil.copy2(source_dir / name, work / name)
            transformed_hashes = {}
            for name in DEPENDENCIES:
                transformed_hashes[name] = install_dependency(args.dependency_dir / name, work / name)
            receipt["temporary_dependency_sha256"] = transformed_hashes
            receipt["temporary_replay_source_sha256"] = {name: digest(work / name) for name in REPLAY_SOURCES}
            compiler = run(["c++", "--version"], work, 15, receipt["commands"])
            receipt["compiler_version"] = compiler.stdout
            source_hash = digest(work / "007_interval_fejer.cpp")
            # Inputs are generated in the temporary directory, never copied from 007 outputs.
            run([sys.executable, "007_generate_fejer_inputs.py"], work, args.timeout, receipt["commands"])
            input_hash = digest(work / "007_fejer_nodes.bits")
            binary = work / "cube5-007-replay"
            build_label = (
                "c++ -O2 -std=c++20 -fno-fast-math -ffp-contract=off "
                f"-DCUBE7_SOURCE_SHA256={source_hash} -DCUBE7_INPUT_SHA256={input_hash} "
                "-DCUBE7_BUILD_COMMAND=<self-recorded> 007_interval_fejer.cpp -o cube5-007-replay"
            )
            compile_command = [
                "c++", "-O2", "-std=c++20", "-fno-fast-math", "-ffp-contract=off",
                f'-DCUBE7_BUILD_COMMAND="{build_label}"',
                f'-DCUBE7_SOURCE_SHA256="{source_hash}"',
                f'-DCUBE7_INPUT_SHA256="{input_hash}"',
                "007_interval_fejer.cpp", "-o", str(binary),
            ]
            run(compile_command, work, args.timeout, receipt["commands"])
            receipt["binary_sha256"] = digest(binary)
            run([str(binary)], work, args.timeout, receipt["commands"])
            run([sys.executable, "007_final_verify.py"], work, args.timeout, receipt["commands"])
            final = json.loads((work / "007_final_verify.json").read_text(encoding="utf-8"))
            heat = json.loads((work / "007_heat_quadrature.json").read_text(encoding="utf-8"))
            inputs = json.loads((work / "007_fejer_inputs.json").read_text(encoding="utf-8"))
            if not final.get("direct_arb_ball_contained"):
                raise RuntimeError("replayed final interval failed its direct Arb containment check")
            receipt.update({
                "status": "success", "elapsed_seconds": time.monotonic() - began,
                "generated_input_sha256": input_hash,
                "generated_heat_output_sha256": digest(work / "007_heat_quadrature.json"),
                "generated_final_output_sha256": digest(work / "007_final_verify.json"),
                "result_summary": {
                    "input_nodes_sha256": inputs["nodes_bits_sha256"],
                    "quadrature_bits": heat["quadrature"],
                    "H_interval_bits": final["H_interval_bits"],
                    "direct_arb_ball_contained": final["direct_arb_ball_contained"],
                    "H_interval_decimal_display": final["H_interval_decimal_display"],
                },
                "scope": "One clean finite replay only. It does not replay 005b throughput, prove spatial uniformity, or establish g/I/C5.",
            })
    except Exception as error:  # Receipt is useful evidence of a finite failed replay.
        receipt.update({"status": "failed", "elapsed_seconds": time.monotonic() - began,
                        "error_type": type(error).__name__, "error": str(error)})
        write_receipt(output, receipt)
        return 1
    write_receipt(output, receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
