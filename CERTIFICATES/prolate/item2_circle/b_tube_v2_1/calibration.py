#!/usr/bin/env python3
"""Process-separated entry points for B-TUBE v2.1 calibration."""
from calibration_context import *
from calibration_config import *
from calibration_security import *
from calibration_numeric import *
from calibration_candidate import *
from calibration_runner import *
from calibration_verify import *
from calibration_delivery import *
from calibration_receipt import *


def _main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--out", type=Path, required=True)
    run_parser.add_argument("--diagnostic", action="store_true")
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--phase", choices=("pre", "final"), required=True)
    verify_parser.add_argument("--source-head", required=True)
    verify_parser.add_argument("--out", type=Path)
    verify_parser.add_argument("--archive", type=Path)
    verify_parser.add_argument("--receipt", type=Path)
    deliver_parser = subparsers.add_parser("deliver")
    deliver_parser.add_argument("--out", type=Path, required=True)
    deliver_parser.add_argument("--delivery", type=Path, required=True)
    deliver_parser.add_argument("--source-head", required=True)
    subparsers.add_parser("verify-config")
    args = parser.parse_args()
    if args.command == "run":
        return run_calibration(args.out, diagnostic=args.diagnostic)
    if args.command == "deliver":
        return deliver(args.out, args.delivery, args.source_head)
    if args.command == "verify-config":
        return verify_config_only()
    if args.phase == "pre":
        if args.out is None:
            raise CalibrationError("--out is required for pre verification")
        return verify_pre(args.out, args.source_head)
    if args.archive is None or args.receipt is None:
        raise CalibrationError("--archive and --receipt are required for final verification")
    return verify_final(args.archive, args.receipt, args.source_head)


if __name__ == "__main__":
    try:
        raise SystemExit(_main())
    except (CalibrationError, SchemaError, OSError, ValueError, KeyError, zipfile.BadZipFile) as exc:
        print(f"CALIBRATION ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
