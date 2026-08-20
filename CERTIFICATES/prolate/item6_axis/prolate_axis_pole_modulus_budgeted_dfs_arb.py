#!/usr/bin/env python3
"""Run the budget-adaptive pole certificate with depth-first compact covers."""
from __future__ import annotations

import prolate_axis_pole_modulus_budgeted_arb as budgeted
import prolate_axis_pole_modulus_compact_dfs_arb as compact_dfs


def main() -> None:
    original_audit = budgeted.endpoint.endpoint_factorization_audit

    def json_safe_audit() -> dict:
        result = original_audit()
        result["checks"] = {
            key: bool(value) for key, value in result.get("checks", {}).items()
        }
        result["status"] = (
            "PASSED" if all(result["checks"].values()) else "FAILED"
        )
        return result

    budgeted.compact = compact_dfs
    budgeted.endpoint.endpoint_factorization_audit = json_safe_audit
    budgeted.main()


if __name__ == "__main__":
    main()
