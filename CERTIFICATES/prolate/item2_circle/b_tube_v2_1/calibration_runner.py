"""Calibration runner with an explicit nonbinding diagnostic profile."""
from calibration_context import *
from calibration_candidate import *
from calibration_config import *
from calibration_numeric import *
from calibration_security import *


def run_calibration(out_dir: Path, *, diagnostic: bool = False) -> int:
    assert_no_stale_inputs(out_dir)
    assert_clean_source_tree()
    assert_workflow_security()
    config, config_raw = load_config()
    if diagnostic:
        start = require_diagnostic_mode(config)
    else:
        require_blocal_dependency(config)
        start = Rational.from_json(
            config["blocal_dependency"]["lambda_start"], "blocal_dependency.lambda_start"
        )
    kernel, kernel_path = load_production_kernel()
    from flint import arb, ctx
    ctx.dps = config["dps"]
    out_dir.mkdir(parents=True)
    (out_dir / "config.calibration.json").write_bytes(config_raw)
    records = []
    previous = chain_genesis(CHAIN_DOMAIN)
    first_passing = None
    pairs = _candidate_pairs(config)
    for candidate_index, (width, radius) in enumerate(pairs):
        passed, previous, candidate = _candidate_run(
            config=config, kernel=kernel, arb_type=arb, start=start,
            width=width, radius=radius, candidate_index=candidate_index,
            records=records, previous=previous,
        )
        if passed and first_passing is None:
            first_passing = candidate

    if diagnostic:
        recommendation = None
        state = "CALIBRATION_INCOMPLETE"
        coverage_claim = False
    else:
        recommendation = first_passing
        state = "CALIBRATION_COMPLETE" if recommendation is not None else "CALIBRATION_INCOMPLETE"
        coverage_claim = recommendation is not None
    summary = {
        "binding_to_final_lambda_start": config["binding_to_final_lambda_start"],
        "candidate_count": len(pairs),
        "chain_tip": previous,
        "coverage_claim": coverage_claim,
        "machine_conclusion": {"real_analytic": False},
        "mode": config["mode"],
        "recommendation": recommendation,
        "record_count": len(records),
        "schema": "btube-calibration-summary-v1",
        "state": state,
    }
    assert_result_namespace(summary)
    (out_dir / "calibration_records.jsonl").write_bytes(canonical_jsonl(records))
    (out_dir / "CALIBRATION_SUMMARY.json").write_bytes(canonical_json_bytes(summary))
    source_manifest = {
        "audited_source_commit": AUDITED_SOURCE_COMMIT,
        "binding_to_final_lambda_start": config["binding_to_final_lambda_start"],
        "design_commit": DESIGN_COMMIT,
        "kernel_path": kernel_path.relative_to(REPO_ROOT).as_posix(),
        "kernel_sha256": sha256_hex(kernel_path.read_bytes()),
        "mode": config["mode"],
        "schema": "btube-calibration-source-manifest-v1",
    }
    (out_dir / "SOURCE_MANIFEST.json").write_bytes(canonical_json_bytes(source_manifest))
    return 0


__all__ = [name for name in globals() if not name.startswith("__")]
