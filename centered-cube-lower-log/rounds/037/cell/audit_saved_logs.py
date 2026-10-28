"""Read saved logs only; no sampling, source evaluation, or file replacement."""
from pathlib import Path
import hashlib
import json
import numpy as np
from replay_stable_ess import stable_weak_scan

HERE = Path(__file__).resolve().parent


def main():
    final = json.loads((HERE / "RESULTS_stable_ESS.json").read_text())
    revision = json.loads((HERE / "ESS_revision.json").read_text())
    hashes = {
        "probe.py": "original_source_SHA256",
        "registration.json": "registration_SHA256",
        "RESULTS.json": "original_RESULTS_SHA256",
        "replay_stable_ess.py": "fix_script_SHA256",
        "ESS_FIX.patch": "minimal_patch_SHA256",
    }
    for name, key in hashes.items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == revision[key]
    checks = []
    for case, saved in zip(final["cases"], final["same_data_protocol_replay_checks"]):
        path = HERE / saved["local_raw_logdata"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == saved["raw_logdata_SHA256"]
        with np.load(path) as data:
            m, d = data["logM"], data["logD"]
        stream = hashlib.sha256()
        for first in range(0, len(d), 64):
            stream.update(np.ascontiguousarray(m[first:first+64]).tobytes())
            stream.update(np.ascontiguousarray(d[first:first+64]).tobytes())
        assert stream.hexdigest() == saved["logM_logD_batch_stream_SHA256"]
        underflow_prefixes = {}
        for i, name in enumerate(case["models"]):
            order = np.argsort(m[:, i])[::-1]
            with np.errstate(under="ignore", invalid="ignore", divide="ignore"):
                w = np.exp(-d[order])
                legacy_ess = np.cumsum(w)**2 / np.cumsum(w*w)
            underflow_prefixes[name] = int(np.count_nonzero(~np.isfinite(legacy_ess)))
            for key, minimum in [("empirical_unrestricted", 1), ("empirical_ESS_ge64", 64)]:
                recomputed = stable_weak_scan(m[:, i], d, minimum)
                assert recomputed == case["models"][name][key]
        checks.append({"n": case["metadata"]["n"], "R": case["metadata"]["R"],
                       "all_saved_scans_exactly_reproduced": True,
                       "legacy_nonfinite_ESS_prefix_count": underflow_prefixes})
    print(json.dumps({"mode": "read-only saved-log audit; zero new samples",
                      "all_original_hashes_match": True, "cases": checks}, indent=2))


if __name__ == "__main__":
    main()
