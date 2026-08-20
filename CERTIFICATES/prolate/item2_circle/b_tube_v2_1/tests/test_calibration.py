#!/usr/bin/env python3
"""Compatibility entry point for the independent layout verifier."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from record_layout_verifier import _main
if __name__ == "__main__":
    raise SystemExit(_main())
