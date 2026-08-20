from __future__ import annotations
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import calibration
from numeric_schema import parse_canonical_json_bytes, sha256_hex


class CalibrationGuardTests(unittest.TestCase):
    def test_all_repository_python_sources_self_scan_clean(self):
        calibration.assert_clean_source_tree()

    def test_all_python_self_scan_detects_forbidden_path(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            token = "flo" + "at("
            (root / "bad.py").write_text("x = " + token + "1)")
            with self.assertRaises(calibration.CalibrationError):
                calibration.assert_clean_source_tree(root)

    def test_result_namespace_rejects_production_prefix(self):
        with self.assertRaises(calibration.CalibrationError):
            calibration.assert_result_namespace({"state": "CERT" + "IFIED_X"})

    def test_result_namespace_rejects_production_key(self):
        with self.assertRaises(calibration.CalibrationError):
            calibration.assert_result_namespace({"verdict": "x"})

    def test_stale_output_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "existing"
            path.mkdir()
            with self.assertRaises(calibration.CalibrationError):
                calibration.assert_no_stale_inputs(path)

    def test_symlink_dependency_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target.py"
            target.write_text("pass\n")
            link = root / "link.py"
            link.symlink_to(target)
            with self.assertRaises(calibration.CalibrationError):
                calibration._assert_repo_regular_file(link, root)

    def test_binding_calibration_run_blocked_until_blocal_is_pinned(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "run"
            with self.assertRaisesRegex(
                calibration.CalibrationError, "B-LOCAL/B-ENTRY dependency is not pinned"
            ):
                calibration.run_calibration(out)
            self.assertFalse(out.exists())

    def test_diagnostic_mode_is_explicit_and_nonbinding(self):
        config = calibration.load_config()[0]
        self.assertEqual(calibration.require_diagnostic_mode(config), calibration.Rational(21, 10))
        self.assertIs(config["binding_to_final_lambda_start"], False)

    def test_affine_rule_is_frozen(self):
        config = calibration.load_config()[0]
        self.assertEqual(config["q_evaluation_rule"], "exact_endpoint_convex_hull_v1")

    def test_workflow_has_tag_head_guard_and_no_dispatch(self):
        calibration.assert_workflow_security()

    def test_workflow_has_independent_unpinned_binding_gate(self):
        text = calibration.WORKFLOW_PATH.read_text(encoding="utf-8")
        marker = "B-LOCAL/B-ENTRY unpinned: workflow binding run prohibited"
        self.assertIn(marker, text)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "workflow.yml"
            path.write_text(text.replace(marker, "binding gate removed"), encoding="utf-8")
            with self.assertRaisesRegex(
                calibration.CalibrationError, "workflow security/authorization guard missing"
            ):
                calibration.assert_workflow_security(path)

    def test_workflow_does_not_authorize_diagnostic_execution(self):
        text = calibration.WORKFLOW_PATH.read_text(encoding="utf-8")
        self.assertNotIn("--diagnostic", text)

    def test_result_merge_rejects_surviving_workflow(self):
        with self.assertRaises(calibration.CalibrationError):
            calibration.assert_no_workflow_in_result_merge(
                [".github/workflows/prolate-item2-btube-v2-1-calibration.yml"]
            )

    def test_receipt_noncanonical_bytes_rejected(self):
        with self.assertRaises(ValueError):
            parse_canonical_json_bytes(b'{"a":1, "b":2}')

    def test_payload_digest_detects_change(self):
        original = b"alpha"
        recorded = sha256_hex(original)
        self.assertNotEqual(sha256_hex(original + b"x"), recorded)


if __name__ == "__main__":
    unittest.main()
