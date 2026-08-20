"""Independent record and pre-delivery verification."""
from calibration_context import *
from calibration_config import *
from calibration_security import *


def _verifier_candidate_pairs(config: dict[str, Any]) -> list[tuple[Dyadic, Dyadic]]:
    """Reconstruct configured candidate order without the runner helper."""
    width_items = config.get("candidate_lambda_widths")
    radius_items = config.get("candidate_tube_radii")
    if not isinstance(width_items, list) or not width_items:
        raise CalibrationError("verifier: width candidates missing")
    if not isinstance(radius_items, list) or not radius_items:
        raise CalibrationError("verifier: radius candidates missing")
    widths = [
        Dyadic.from_json(item, f"candidate_lambda_widths[{index}]")
        for index, item in enumerate(width_items)
    ]
    radii = [
        Dyadic.from_json(item, f"candidate_tube_radii[{index}]")
        for index, item in enumerate(radius_items)
    ]
    for name, values in (("width", widths), ("radius", radii)):
        if any(value <= D_ZERO for value in values):
            raise CalibrationError(f"verifier: {name} candidates must be positive")
        if len(set(values)) != len(values):
            raise CalibrationError(f"verifier: duplicate {name} candidate")
        if any(not values[index + 1] < values[index] for index in range(len(values) - 1)):
            raise CalibrationError(f"verifier: {name} candidates not strictly decreasing")
    pairs = []
    for width in widths:
        for radius in radii:
            pairs.append((width, radius))
    return pairs


def _verify_records(out_dir: Path):
    config, config_raw = load_config(out_dir / "config.calibration.json")
    parsed = parse_canonical_jsonl((out_dir / "calibration_records.jsonl").read_bytes())
    previous = chain_genesis(CHAIN_DOMAIN)
    for record, raw in parsed:
        if record.get("previous_record_sha256") != previous:
            raise CalibrationError("record chain mismatch")
        assert_result_namespace(record)
        previous = sha256_hex(raw)
    summary = parse_canonical_json_bytes(
        (out_dir / "CALIBRATION_SUMMARY.json").read_bytes(), allow_display=False,
    )
    _require_exact_keys(summary, {
        "binding_to_final_lambda_start", "candidate_count", "chain_tip", "coverage_claim",
        "machine_conclusion", "mode", "recommendation", "record_count", "schema", "state",
    }, "summary")
    assert_result_namespace(summary)
    if summary["schema"] != "btube-calibration-summary-v1":
        raise CalibrationError("summary schema mismatch")
    if summary["machine_conclusion"] != {"real_analytic": False}:
        raise CalibrationError("machine_conclusion must be exactly present-and-false")
    if summary["state"] not in TERMINAL_STATES:
        raise CalibrationError("invalid terminal state")
    if summary["chain_tip"] != previous or summary["record_count"] != len(parsed):
        raise CalibrationError("summary chain/count mismatch")
    if summary["mode"] != config["mode"]:
        raise CalibrationError("summary mode mismatch")
    if summary["binding_to_final_lambda_start"] is not config["binding_to_final_lambda_start"]:
        raise CalibrationError("summary binding flag mismatch")

    ends = [record for record, _ in parsed if record.get("record_type") == "candidate_end"]
    pairs = _verifier_candidate_pairs(config)
    if len(ends) != summary["candidate_count"] or len(ends) != len(pairs):
        raise CalibrationError("candidate completeness mismatch")
    if [record.get("candidate_index") for record in ends] != list(range(len(pairs))):
        raise CalibrationError("candidate order/index mismatch")
    passing = [record["candidate_index"] for record in ends if record.get("passed") is True]
    first_passing = None
    if passing:
        first = passing[0]
        width, radius = pairs[first]
        first_passing = {
            "candidate_index": first,
            "lambda_width": width.to_json(),
            "tube_radius": radius.to_json(),
        }

    if config["mode"] == CALIBRATION_MODE:
        expected = None
        expected_state = "CALIBRATION_INCOMPLETE"
        expected_coverage = False
        if config["binding_to_final_lambda_start"] is not False:
            raise CalibrationError("diagnostic mode cannot bind final lambda_start")
    else:
        expected = first_passing
        expected_state = "CALIBRATION_COMPLETE" if expected is not None else "CALIBRATION_INCOMPLETE"
        expected_coverage = expected is not None
    if (summary["recommendation"] != expected
            or summary["state"] != expected_state
            or summary["coverage_claim"] is not expected_coverage):
        raise CalibrationError("deterministic recommendation/state policy mismatch")
    return config, summary, config_raw


def verify_pre(out_dir: Path, source_head: str) -> int:
    assert_clean_source_tree()
    assert_workflow_security()
    config, summary, config_raw = _verify_records(out_dir)
    require_blocal_dependency(config)
    load_production_kernel()
    report = {
        "config_sha256": sha256_hex(config_raw),
        "kernel_sha256": KERNEL_SHA256,
        "record_chain_tip": summary["chain_tip"],
        "schema": "btube-calibration-checker-report-v1",
        "source_head": source_head,
        "state": summary["state"],
        "verifier": "PASS",
    }
    assert_result_namespace(report)
    (out_dir / "CHECKER_REPORT.json").write_bytes(canonical_json_bytes(report))
    return 0


__all__ = [name for name in globals() if not name.startswith("__")]
