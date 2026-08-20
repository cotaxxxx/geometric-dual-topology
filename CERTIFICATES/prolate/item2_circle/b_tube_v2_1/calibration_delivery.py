"""Byte-closed deterministic payload delivery."""
from calibration_context import *
from calibration_config import *
from calibration_security import *
from calibration_verify import *

def _payload_files(run_dir: Path) -> dict[str, Path]:
    files = {
        "CALIBRATION_SUMMARY.json": run_dir / "CALIBRATION_SUMMARY.json",
        "CHECKER_REPORT.json": run_dir / "CHECKER_REPORT.json",
        "SOURCE_MANIFEST.json": run_dir / "SOURCE_MANIFEST.json",
        "calibration_records.jsonl": run_dir / "calibration_records.jsonl",
        "config.calibration.json": run_dir / "config.calibration.json",
        "source/.github/workflows/prolate-item2-btube-v2-1-calibration.yml": WORKFLOW_PATH,
        f"source/{KERNEL_RELATIVE.as_posix()}": REPO_ROOT / KERNEL_RELATIVE,
    }
    for relative in SOURCE_FILE_LIST:
        files[f"source/CERTIFICATES/prolate/item2_circle/b_tube_v2_1/{relative}"] = BTUBE_ROOT / relative
    return files

def _build_deterministic_zip(payload_dir: Path, archive_path: Path) -> None:
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(item for item in payload_dir.rglob("*") if item.is_file()):
            relative = path.relative_to(payload_dir).as_posix()
            info = zipfile.ZipInfo(relative, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, path.read_bytes())

def deliver(run_dir: Path, delivery_dir: Path, source_head: str) -> int:
    if delivery_dir.exists():
        raise CalibrationError("delivery directory must not exist")
    config, summary, config_raw = _verify_records(run_dir)
    require_blocal_dependency(config)
    checker = parse_canonical_json_bytes((run_dir / "CHECKER_REPORT.json").read_bytes())
    if checker.get("verifier") != "PASS" or checker.get("source_head") != source_head:
        raise CalibrationError("pre-verifier report mismatch")
    delivery_dir.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix="btube-calibration-payload-") as temporary:
        payload_dir = Path(temporary) / "payload"
        payload_dir.mkdir()
        for relative, source in sorted(_payload_files(run_dir).items()):
            target = payload_dir / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())
        manifest = {
            "files": {
                path.relative_to(payload_dir).as_posix(): sha256_hex(path.read_bytes())
                for path in sorted(item for item in payload_dir.rglob("*") if item.is_file())
            },
            "schema": "btube-calibration-payload-manifest-v1",
        }
        manifest_raw = canonical_json_bytes(manifest)
        manifest_path = payload_dir / "PAYLOAD_SHA256SUMS.json"
        manifest_path.write_bytes(manifest_raw)
        for relative, digest in manifest["files"].items():
            if sha256_hex((payload_dir / relative).read_bytes()) != digest:
                raise CalibrationError("payload changed after manifest creation")
        if manifest_path.read_bytes() != manifest_raw:
            raise CalibrationError("payload manifest byte mismatch")
        archive_path = delivery_dir / "btube-v2-1-calibration.zip"
        _build_deterministic_zip(payload_dir, archive_path)
        receipt = {
            "archive_name": archive_path.name,
            "archive_sha256": sha256_hex(archive_path.read_bytes()),
            "configuration_sha256": sha256_hex(config_raw),
            "kernel_file_sha256": KERNEL_SHA256,
            "payload_manifest_sha256": sha256_hex(manifest_raw),
            "schema": "btube-calibration-delivery-receipt-v1",
            "source_head": source_head,
            "state": summary["state"],
            "workflow_source_sha256": sha256_hex(WORKFLOW_PATH.read_bytes()),
        }
        assert_result_namespace(receipt)
        receipt_raw = canonical_json_bytes(receipt)
        receipt_path = delivery_dir / "DELIVERY_RECEIPT.json"
        receipt_path.write_bytes(receipt_raw)
        if receipt_path.read_bytes() != canonical_json_bytes(receipt):
            raise CalibrationError("receipt canonical-byte mismatch")
        if sha256_hex(archive_path.read_bytes()) != receipt["archive_sha256"]:
            raise CalibrationError("archive changed after receipt creation")
    return 0

__all__ = [name for name in globals() if not name.startswith("__")]
