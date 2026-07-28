from __future__ import annotations
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import calibration
from numeric_schema import Rational, canonical_json_bytes


class CalibrationConfigTests(unittest.TestCase):
    def _config(self):
        return calibration.load_config()[0]

    def _write(self, directory: Path, obj) -> Path:
        path = directory / "config.json"
        path.write_bytes(canonical_json_bytes(obj))
        return path

    def test_valid_diagnostic_profile_and_precision_equality(self):
        config, raw = calibration.load_config()
        self.assertEqual(config["checker_dps"], config["dps"])
        self.assertEqual(config["mode"], "DIAGNOSTIC_ONLY")
        self.assertIs(config["binding_to_final_lambda_start"], False)
        self.assertNotIn("lambda_start", config)
        self.assertEqual(raw, canonical_json_bytes(config))

    def test_diagnostic_start_is_exact_and_above_stage1_upper(self):
        config = self._config()
        start = Rational.from_json(config["diagnostic_lambda_start"])
        self.assertEqual(start, Rational(21, 10))
        self.assertLess(calibration.BLOCAL_STAGE1_UPPER, start)
        self.assertLess(start, Rational.from_json(config["lambda_end"]))

    def test_unpinned_blocal_tuple_is_explicit_and_null(self):
        dependency = self._config()["blocal_dependency"]
        self.assertEqual(dependency["status"], "UNPINNED")
        for key in (
            "artifact_zip_sha256", "certificate_sha256", "config_sha256",
            "lambda_start", "machine_conclusion", "source_head",
        ):
            self.assertIsNone(dependency[key])

    def test_false_blocal_tuple_promotion_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["blocal_dependency"]["status"] = "PINNED"
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_nonnull_unpinned_blocal_field_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["blocal_dependency"]["source_head"] = "0" * 40
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_binding_flag_must_remain_false(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["binding_to_final_lambda_start"] = True
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_diagnostic_start_not_safely_above_boundary_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["diagnostic_lambda_start"] = {"p": "206539", "q": "100000"}
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_blocal_dependency_gate_is_fail_closed(self):
        with self.assertRaisesRegex(
            calibration.CalibrationError, "B-LOCAL/B-ENTRY dependency is not pinned"
        ):
            calibration.require_blocal_dependency(self._config())

    def test_diagnostic_gate_returns_nonbinding_start(self):
        start = calibration.require_diagnostic_mode(self._config())
        self.assertEqual(start, Rational(21, 10))

    def test_kernel_pin_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["production_kernel_sha256"] = "0" * 64
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_cg_tuple_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["cg_match_dependency"]["source_head"] = "0" * 40
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_lambda_end_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["lambda_end"] = {"p": "19", "q": "4"}
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_checker_precision_not_below_runner(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["checker_dps"] = config["dps"] - 1
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_duplicate_candidate_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["candidate_lambda_widths"][1] = config["candidate_lambda_widths"][0]
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_unordered_candidate_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["candidate_tube_radii"] = list(reversed(config["candidate_tube_radii"]))
            with self.assertRaises(calibration.CalibrationError):
                calibration.load_config(self._write(Path(temporary), config))

    def test_floating_json_number_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = self._config()
            config["evaluation_budget"] = 1.5
            path = Path(temporary) / "config.json"
            path.write_text(json.dumps(config, sort_keys=True, separators=(",", ":")))
            with self.assertRaises(ValueError):
                calibration.load_config(path)


if __name__ == "__main__":
    unittest.main()
