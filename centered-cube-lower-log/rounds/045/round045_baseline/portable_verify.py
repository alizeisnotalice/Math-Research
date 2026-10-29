#!/usr/bin/env python3
"""Replay the unchanged frozen parameter probe in a temporary dependency layout.

Only --output is written persistently. Original script, result and inputs are read
only. Exact JSON equality is required; there is no floating tolerance fallback.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
FROZEN = {
    "script": "305dc70c0bef886e8ae78f893dc3e6a56211d0f6184d2cccf4561a1d72edc9c0",
    "result": "411bf37572183b587c69e71b6cbddcffd7b8345a2dbacb1f47f237c9ce3e5e1d",
    "paper": "23f4b15c1e21815721c7839eabe5cd66bf6760122c20359c41dbdfd2ac977222",
    "labels": "fa741ce7f68780023c870b22248dc44e029b10343ab1614d09043282820278b8",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(a, b, path="$", differences=None):
    if differences is None:
        differences = []
    if type(a) is not type(b):
        differences.append(dict(path=path, expected_type=type(a).__name__, actual_type=type(b).__name__))
    elif isinstance(a, dict):
        if set(a) != set(b):
            differences.append(dict(path=path, expected_keys=sorted(a), actual_keys=sorted(b)))
        for key in a.keys() & b.keys():
            compare(a[key], b[key], path + "." + key, differences)
    elif isinstance(a, list):
        if len(a) != len(b):
            differences.append(dict(path=path, expected_length=len(a), actual_length=len(b)))
        for i, (aa, bb) in enumerate(zip(a, b)):
            compare(aa, bb, path + "[" + str(i) + "]", differences)
    elif a != b:
        differences.append(dict(path=path, expected=a, actual=b))
    return differences


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paper", type=Path, required=True, help="Unmodified accepted paper.tex")
    ap.add_argument("--labels", type=Path, required=True, help="Unmodified source_d10.npz")
    ap.add_argument("--output", type=Path, required=True, help="New verification JSON; must not exist")
    ap.add_argument("--script", type=Path, default=HERE / "parameter_probe.py")
    ap.add_argument("--expected", type=Path, default=HERE / "PARAMETER_PROBE.json")
    args = ap.parse_args()
    output = args.output.resolve()
    if output.exists():
        ap.error("--output already exists; choose a new verification path")
    inputs = dict(script=args.script.resolve(), result=args.expected.resolve(),
                  paper=args.paper.resolve(), labels=args.labels.resolve())
    for name, path in inputs.items():
        if not path.is_file():
            ap.error(name + " is not a regular file: " + str(path))
        if output == path:
            ap.error("--output must differ from every input")
    before = {name: sha(path) for name, path in inputs.items()}
    if before != FROZEN:
        ap.error("input identity differs from frozen round045 bytes: " +
                 json.dumps({k: dict(expected=FROZEN[k], actual=v)
                             for k, v in before.items() if v != FROZEN[k]}))
    expected = json.loads(inputs["result"].read_text(encoding="utf-8"))
    assert expected["source_paper_SHA256"] == before["paper"]
    assert expected["source_sampler_probe"]["label_input_SHA256"] == before["labels"]
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="hl-round045-portable-") as tmp:
        root = Path(tmp)
        work = root / "hl_cube_log_iteration_20261029" / "round045_baseline"
        paper = root / "hl_cube_paper_20261013" / "paper.tex"
        labels = paper.parent / "numerics" / "inputs" / "source_d10.npz"
        work.mkdir(parents=True)
        labels.parent.mkdir(parents=True)
        for src, dst in [(inputs["script"], work / "parameter_probe.py"),
                         (inputs["paper"], paper), (inputs["labels"], labels)]:
            shutil.copyfile(src, dst)
            dst.chmod(0o444)
        # The frozen script reads HERE.parent.parent/...; this temporary layout
        # satisfies both external paths without editing any script byte.
        child_start = time.perf_counter()
        child = subprocess.run([sys.executable, str(work / "parameter_probe.py")],
                               cwd=work, capture_output=True, text=True, timeout=120)
        child_seconds = time.perf_counter() - child_start
        if child.returncode:
            raise RuntimeError("frozen probe failed with exit " + str(child.returncode) +
                               ":\n" + child.stderr[-12000:])
        replay_file = work / "PARAMETER_PROBE.json"
        replay = json.loads(replay_file.read_text(encoding="utf-8"))
        assert json.loads(child.stdout) == replay, "stdout/file payload mismatch"
        replay_sha = sha(replay_file)
        differences = compare(expected, replay)
    after = {name: sha(path) for name, path in inputs.items()}
    assert after == before, "original input bytes changed during replay"
    sampler = replay["source_sampler_probe"]
    try:
        import numpy
        numpy_version = numpy.__version__
    except ImportError:
        numpy_version = "unavailable in verifier process"
    report = dict(
        status="PASS" if not differences else "MISMATCH",
        exact_json_equal=not differences,
        byte_identical_result=replay_sha == before["result"],
        expected_result_SHA256=before["result"], replay_result_SHA256=replay_sha,
        frozen_input_SHA256=before, original_inputs_unchanged=after == before,
        ignored_payload_fields=[], float_tolerance_used=False,
        differences=differences,
        execution=dict(python=sys.version, executable=sys.executable,
                       numpy=numpy_version, platform=platform.platform(),
                       child_wall_seconds=child_seconds,
                       wrapper_wall_seconds=time.perf_counter() - started,
                       timeout_seconds=120),
        sampler_replayed=dict(N=sampler["N"], seed=sampler["seed"],
                              peak_draws=sampler["peak_draws"],
                              background_proposals=sampler["background_proposals"],
                              max_empirical_CDF_grid_error=sampler["max_empirical_CDF_grid_error"]),
        parameter_rows_replayed=len(replay["rows"]),
        scope="Replays parameter arithmetic and finite unmatched source/CDF diagnostic only. Necessary construction gates are not full source qualification, an explicit n0, or a finite C_n certificate.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(json.dumps(report, indent=2))
    return 0 if not differences else 1


if __name__ == "__main__":
    raise SystemExit(main())
