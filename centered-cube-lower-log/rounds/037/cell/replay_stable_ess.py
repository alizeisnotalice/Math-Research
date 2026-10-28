"""Postprocessing correction, same frozen data protocol; NOT a new validation."""
from pathlib import Path
import hashlib
import json
import inspect
import difflib
import time
import math
import sys
import numpy as np
import probe as original

HERE=Path(__file__).resolve().parent


def stable_weak_scan(log_values,log_d,ess_min):
    order=np.argsort(log_values)[::-1]
    v=log_values[order];log_w=-log_d[order]
    log_sum=np.logaddexp.accumulate(log_w)
    log_square_sum=np.logaddexp.accumulate(2*log_w)
    log_ess=2*log_sum-log_square_sum
    score=v+log_sum-math.log(len(v))
    selected=log_ess>=math.log(ess_min)-1e-13
    ix=int(np.argmax(np.where(selected,score,-np.inf)))
    return {"coefficient":float(np.exp(score[ix])),
            "log_threshold":float(v[ix]),
            "ESS":float(np.exp(log_ess[ix])),"min_ESS":ess_min}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target=HERE/"RESULTS_stable_ESS.json"
    if target.exists():raise RuntimeError("Refusing a second correction replay.")
    initial=json.loads((HERE/"RESULTS.json").read_text())
    registration=json.loads((HERE/"registration.json").read_text())
    assert sha(HERE/"probe.py")==registration["script_SHA256"]
    assert initial["numpy"]==np.__version__
    assert initial["python_version"]==sys.version
    old_function=inspect.getsource(original.weak_scan)
    new_function=inspect.getsource(stable_weak_scan).replace(
        "def stable_weak_scan(","def weak_scan(",1)
    patch="".join(difflib.unified_diff(old_function.splitlines(True),
        new_function.splitlines(True),fromfile="original/probe.py:weak_scan",
        tofile="stable_log_ESS/probe.py:weak_scan"))
    (HERE/"ESS_FIX.patch").write_text(patch)
    revision={"scope":"postprocessing correction, not an independent validation",
        "original_source_SHA256":sha(HERE/"probe.py"),
        "registration_SHA256":sha(HERE/"registration.json"),
        "original_RESULTS_SHA256":sha(HERE/"RESULTS.json"),
        "fix_script_SHA256":sha(Path(__file__)),
        "minimal_patch_SHA256":sha(HERE/"ESS_FIX.patch"),
        "numpy":np.__version__,"python":sys.executable,
        "source_RNG_queries_unchanged":True,
        "replayed_same_registered_seeds":[c["seed"]for c in registration["cases"]],
        "reason":"linear weight squares underflow at n63/R8; stable log cumulative "
                 "sums correct both unrestricted and supported ESS/coefficients",
        "original_raw_logM_logD_saved":False,
        "replay_raw_logM_logD_saved_locally":True}
    (HERE/"ESS_revision.json").write_text(json.dumps(revision,indent=2))
    original.weak_scan=stable_weak_scan
    evaluator=original.evaluate_queries
    state={}
    def capturing_evaluator(src,x):
        m,d=evaluator(src,x)
        state["query_hash"].update(np.ascontiguousarray(x).tobytes())
        state["data_hash"].update(np.ascontiguousarray(m).tobytes())
        state["data_hash"].update(np.ascontiguousarray(d).tobytes())
        state["m"].append(m.copy());state["d"].append(d.copy())
        return m,d
    original.evaluate_queries=capturing_evaluator
    started=time.time();cases=[];checks=[]
    for config,old in zip(registration["cases"],initial["cases"]):
        state={"query_hash":hashlib.sha256(),"data_hash":hashlib.sha256(),
               "m":[],"d":[]}
        new=original.run_case(config,registration["alpha_each"])
        invariant_keys=["metadata","mass_weight_mean","mass_SE_descriptive",
            "weight_max","uniform_rows","source_b_counts_by_model",
            "max_actual_coefficient_contribution","contribution_cap_all_thresholds",
            "max_PB_probability_mass_residual",
            "positive_subnormal_probability_entries",
            "subnormal_loss_magnitude_budget_only"]
        for key in invariant_keys:
            assert new[key]==old[key],(config["n"],config["R"],key)
        unaffected=not(config["n"]==63 and config["R"]==8)
        max_relative_change=0.
        if unaffected:
            for name in old["models"]:
                for scan in ["empirical_unrestricted","empirical_ESS_ge64"]:
                    a=old["models"][name][scan]["coefficient"]
                    b=new["models"][name][scan]["coefficient"]
                    max_relative_change=max(max_relative_change,abs(a-b)/max(1,abs(a)))
            assert max_relative_change<1e-10
        raw_name=f"logdata_n{config['n']}_R{config['R']}.npz"
        np.savez_compressed(HERE/raw_name,logM=np.concatenate(state["m"]),
                            logD=np.concatenate(state["d"]))
        checks.append({"n":config["n"],"R":config["R"],
            "all_source_query_weight_invariants_identical":True,
            "unaffected_case_best_curve_relative_change":max_relative_change
                if unaffected else None,
            "query_coordinates_SHA256":state["query_hash"].hexdigest(),
            "logM_logD_batch_stream_SHA256":state["data_hash"].hexdigest(),
            "local_raw_logdata":raw_name,"raw_logdata_SHA256":sha(HERE/raw_name)})
        cases.append(new)
    corrected={**initial,"cases":cases,
        "runtime_seconds":time.time()-started,
        "original_run_seconds":initial["runtime_seconds"],
        "postprocessing_correction":revision,
        "same_data_protocol_replay_checks":checks,
        "macroscopic_diagnostics_reused_without_replay":True}
    assert sha(HERE/"probe.py")==revision["original_source_SHA256"]
    assert sha(HERE/"RESULTS.json")==revision["original_RESULTS_SHA256"]
    target.write_text(json.dumps(corrected,separators=(",",":")))
    print("CORRECTION_COMPLETE",corrected["runtime_seconds"],flush=True)


if __name__=="__main__":
    main()
