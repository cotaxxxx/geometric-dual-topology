from __future__ import annotations
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import calibration
from numeric_schema import Dyadic, Rational, canonical_json_bytes, canonical_jsonl, chain_genesis, parse_canonical_json_bytes, parse_canonical_jsonl
from record_layout_contract import _partition
from record_layout_verifier import verify_record_layout


def _rechain(records: list[dict]) -> tuple[list[dict], str]:
    chained = []
    previous = chain_genesis(calibration.CHAIN_DOMAIN)
    for original in records:
        body = dict(original)
        body.pop("previous_record_sha256", None)
        previous = calibration._append_record(chained, previous, body)
    return chained, previous


class CalibrationRecordTests(unittest.TestCase):
    def _write_case(self, directory: Path, passes: list[bool]):
        config, config_raw = calibration.load_config()
        pairs = calibration._candidate_pairs(config)
        self.assertEqual(len(passes), len(pairs))
        start = calibration.require_diagnostic_mode(config).as_fraction()
        end = Rational.from_json(config["lambda_end"]).as_fraction()
        records = []
        previous = chain_genesis(calibration.CHAIN_DOMAIN)
        for index, (passed, pair) in enumerate(zip(passes, pairs)):
            width, radius = pair
            previous = calibration._append_record(records, previous, {
                "candidate_index": index,
                "lambda_width": width.to_json(),
                "record_type": "candidate_start",
                "tube_radius": radius.to_json(),
            })
            cells = _partition(start, end, width.as_fraction())
            for cell_index, (left, right) in enumerate(cells):
                previous = calibration._append_record(records, previous, {
                    "candidate_index": index,
                    "cell_index": cell_index,
                    "evaluation_count": 3 * (cell_index + 1),
                    "lambda_interval": {
                        "lo": Rational.from_fraction(left).to_json(),
                        "hi": Rational.from_fraction(right).to_json(),
                    },
                    "passed": passed,
                    "record_type": "cell",
                })
            for join_index in range(max(len(cells) - 1, 0)):
                previous = calibration._append_record(records, previous, {
                    "candidate_index": index,
                    "failure_reason": None,
                    "join_index": join_index,
                    "record_type": "join",
                    "width": Dyadic(1, 20).to_json(),
                })
            previous = calibration._append_record(records, previous, {
                "candidate_index": index,
                "cells_attempted": len(cells),
                "cells_passed": len(cells) if passed else 0,
                "evaluation_count": 3 * len(cells),
                "joins_passed": True,
                "passed": passed,
                "record_type": "candidate_end",
            })
        summary = {
            "binding_to_final_lambda_start": False,
            "candidate_count": len(pairs),
            "chain_tip": previous,
            "coverage_claim": False,
            "machine_conclusion": {"real_analytic": False},
            "mode": "DIAGNOSTIC_ONLY",
            "recommendation": None,
            "record_count": len(records),
            "schema": "btube-calibration-summary-v1",
            "state": "CALIBRATION_INCOMPLETE",
        }
        directory.mkdir()
        (directory / "config.calibration.json").write_bytes(config_raw)
        (directory / "calibration_records.jsonl").write_bytes(canonical_jsonl(records))
        (directory / "CALIBRATION_SUMMARY.json").write_bytes(canonical_json_bytes(summary))

    def _count(self) -> int:
        return len(calibration._candidate_pairs(calibration.load_config()[0]))

    def test_diagnostic_mode_suppresses_passing_candidate_recommendation(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "case"
            passes = [False] * self._count()
            passes[3] = True
            passes[5] = True
            self._write_case(out, passes)
            _, summary, _ = calibration._verify_records(out)
            self.assertEqual(summary["state"], "CALIBRATION_INCOMPLETE")
            self.assertIsNone(summary["recommendation"])
            self.assertIs(summary["coverage_claim"], False)
            layout = verify_record_layout(out, allow_unbound_fixture=True)
            self.assertIsNone(layout["recommendation"])
            self.assertIs(layout["coverage_claim"], False)

    def test_false_diagnostic_recommendation_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "case"
            passes = [False] * self._count()
            passes[0] = True
            self._write_case(out, passes)
            summary = parse_canonical_json_bytes((out / "CALIBRATION_SUMMARY.json").read_bytes())
            width, radius = calibration._candidate_pairs(calibration.load_config()[0])[0]
            summary["recommendation"] = {
                "candidate_index": 0,
                "lambda_width": width.to_json(),
                "tube_radius": radius.to_json(),
            }
            summary["state"] = "CALIBRATION_COMPLETE"
            summary["coverage_claim"] = True
            (out / "CALIBRATION_SUMMARY.json").write_bytes(canonical_json_bytes(summary))
            with self.assertRaises(calibration.CalibrationError):
                calibration._verify_records(out)

    def test_verifiers_rederive_candidate_pairs_without_runner_helper(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "case"
            count = self._count()
            self._write_case(out, [False] * count)
            with patch.object(
                calibration, "_candidate_pairs", side_effect=AssertionError("runner helper used")
            ):
                _, summary, _ = calibration._verify_records(out)
                layout = verify_record_layout(out, allow_unbound_fixture=True)
            self.assertEqual(summary["candidate_count"], count)
            self.assertEqual(layout["candidate_count"], count)

    def test_layout_verifier_blocks_unbound_production_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "case"
            self._write_case(out, [False] * self._count())
            with self.assertRaisesRegex(
                calibration.CalibrationError, "B-LOCAL/B-ENTRY dependency is not pinned"
            ):
                verify_record_layout(out)

    def test_missing_attempted_cell_rejected_after_rechain(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "case"
            self._write_case(out, [False] * self._count())
            parsed = parse_canonical_jsonl((out / "calibration_records.jsonl").read_bytes())
            records = [record for record, _ in parsed]
            missing = next(
                index for index, record in enumerate(records)
                if record.get("record_type") == "cell" and record.get("candidate_index") == 0
            )
            del records[missing]
            chained, tip = _rechain(records)
            (out / "calibration_records.jsonl").write_bytes(canonical_jsonl(chained))
            summary = parse_canonical_json_bytes((out / "CALIBRATION_SUMMARY.json").read_bytes())
            summary["record_count"] = len(chained)
            summary["chain_tip"] = tip
            (out / "CALIBRATION_SUMMARY.json").write_bytes(canonical_json_bytes(summary))
            with self.assertRaises(calibration.CalibrationError):
                verify_record_layout(out, allow_unbound_fixture=True)

    def test_record_chain_tamper_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "case"
            self._write_case(out, [False] * self._count())
            records = (out / "calibration_records.jsonl").read_bytes().split(b"\n")
            record = parse_canonical_json_bytes(records[1])
            record["previous_record_sha256"] = "0" * 64
            records[1] = canonical_json_bytes(record)
            (out / "calibration_records.jsonl").write_bytes(b"\n".join(records))
            with self.assertRaises(calibration.CalibrationError):
                calibration._verify_records(out)

    def test_machine_conclusion_present_and_false(self):
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "case"
            self._write_case(out, [False] * self._count())
            summary = parse_canonical_json_bytes((out / "CALIBRATION_SUMMARY.json").read_bytes())
            summary["machine_conclusion"]["real_analytic"] = True
            (out / "CALIBRATION_SUMMARY.json").write_bytes(canonical_json_bytes(summary))
            with self.assertRaises(calibration.CalibrationError):
                calibration._verify_records(out)


if __name__ == "__main__":
    unittest.main()
