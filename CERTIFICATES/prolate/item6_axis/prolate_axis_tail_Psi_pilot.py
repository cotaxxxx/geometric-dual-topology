#!/usr/bin/env python3
"""Extended-trace wrapper for staged direct-Psi compact-tail pilots.

This wrapper changes only the number of component-level evaluation records kept
in the diagnostic JSON.  The certified arithmetic and subdivision logic remain
those of ``prolate_axis_tail_H_scaled_block_arb.py``.
"""
from __future__ import annotations

import prolate_axis_tail_H_scaled_block_arb as driver


driver.TRACE_LIMIT = 512


if __name__ == "__main__":
    driver.main()
