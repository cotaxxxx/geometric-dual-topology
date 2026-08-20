"""Deterministic evaluation of one calibration candidate."""
from calibration_context import *
from calibration_numeric import *


def _candidate_run(*, config, kernel, arb_type, start, width, radius, candidate_index,
                   records, previous):
    start_fraction = start.as_fraction()
    end = Rational.from_json(config["lambda_end"]).as_fraction()
    cells = _cell_partition(start_fraction, end, width.as_fraction(), config["max_cells"])
    tol = "1e-20"
    depth = config["max_subdivisions"]
    limit = config["evaluation_budget"]
    y_box = DyadicInterval(-radius, radius)
    anchor = _nearest_dyadic((CG_ROOT[0].as_fraction() + CG_ROOT[1].as_fraction()) / 2)
    predictors_reversed = []
    seed = anchor
    refresh = config["predictor_refresh"]
    for reverse_index, (left, right) in enumerate(reversed(cells)):
        right_iterations = 4 if reverse_index % refresh == 0 else 1
        q_right = _newton_predictor(
            kernel, arb_type, right, seed, iterations=right_iterations,
            tol=tol, depth=depth, limit=limit,
        )
        q_left = _newton_predictor(
            kernel, arb_type, left, q_right, iterations=1,
            tol=tol, depth=depth, limit=limit,
        )
        predictor = AffinePredictor(
            Rational.from_fraction(left), Rational.from_fraction(right), q_left, q_right,
        )
        predictors_reversed.append((left, right, predictor))
        seed = q_left
    predictors = list(reversed(predictors_reversed))

    previous = _append_record(records, previous, {
        "candidate_index": candidate_index, "lambda_width": width.to_json(),
        "record_type": "candidate_start", "tube_radius": radius.to_json(),
    })
    cell_passes = []
    joins_pass = True
    evaluation_count = 0
    sections = []
    for cell_index, (left, right, predictor) in enumerate(predictors):
        domain = physical_tube(predictor.range_hull(), y_box)
        reason = None
        image = DyadicInterval.point(domain.midpoint())
        residual = DyadicInterval.point(D_ZERO)
        slope = DyadicInterval.point(D_ZERO)
        preconditioner = D_ZERO
        left_margin = D_ZERO
        right_margin = D_ZERO
        passed = False
        if domain.lo <= D_ZERO or not domain.hi < Dyadic(1, 0):
            reason = "physical_tube_outside_open_unit_interval"
        else:
            lam_box = _fraction_box(left, right, arb_type)
            domain_box = _dyadic_box(domain, arb_type)
            midpoint = domain.midpoint()
            midpoint_lam = (left + right) / 2
            residual = arb_ball_to_exact_interval(kernel.F_arb(
                _dyadic_arb(midpoint, arb_type), lam_box,
                tol=tol, depth=depth, limit=limit,
            ))
            slope = arb_ball_to_exact_interval(kernel.dFdr_arb(
                domain_box, lam_box, tol=tol, depth=depth, limit=limit,
            ))
            center_slope = arb_ball_to_exact_interval(kernel.dFdr_arb(
                _dyadic_arb(midpoint, arb_type), _rational_arb(midpoint_lam, arb_type),
                tol=tol, depth=depth, limit=limit,
            ))
            evaluation_count += 3
            preconditioner = center_slope.midpoint()
            if preconditioner == D_ZERO:
                reason = "preconditioner_zero"
            else:
                image = krawczyk_image(
                    m=midpoint, residual=residual, slope=slope,
                    preconditioner=preconditioner, domain=domain,
                )
                left_margin = image.lo - domain.lo
                right_margin = domain.hi - image.hi
                if not domain.strictly_contains(image):
                    reason = "krawczyk_not_strict"
                elif not slope.hi < D_ZERO:
                    reason = "slope_not_strictly_negative"
                else:
                    passed = True
        cell_passes.append(passed)
        sections.append((predictor, y_box))
        previous = _append_record(records, previous, {
            "candidate_index": candidate_index,
            "cell_index": cell_index,
            "evaluation_count": evaluation_count,
            "failure_reason": reason,
            "krawczyk_image": image.to_json(),
            "lambda_interval": {
                "lo": Rational.from_fraction(left).to_json(),
                "hi": Rational.from_fraction(right).to_json(),
            },
            "left_margin": left_margin.to_json(),
            "passed": passed,
            "predictor": {
                "q_left": predictor.q_left.to_json(),
                "q_right": predictor.q_right.to_json(),
                "rule": Q_RULE,
            },
            "preconditioner": preconditioner.to_json(),
            "record_type": "cell",
            "residual": residual.to_json(),
            "right_margin": right_margin.to_json(),
            "slope": slope.to_json(),
            "subdivision_count": 0,
            "tube_interval": domain.to_json(),
        })

    for join_index in range(len(sections) - 1):
        left_predictor, left_y = sections[join_index]
        right_predictor, right_y = sections[join_index + 1]
        failure = None
        width_value = D_ZERO
        try:
            intersection = exact_join_intersection(
                left_predictor.q_right, left_y, right_predictor.q_left, right_y,
            )
            width_value = intersection.hi - intersection.lo
        except SchemaError:
            intersection = DyadicInterval.point(D_ZERO)
            failure = "join_empty_or_zero_width"
            joins_pass = False
        previous = _append_record(records, previous, {
            "candidate_index": candidate_index,
            "failure_reason": failure,
            "intersection": intersection.to_json(),
            "join_index": join_index,
            "record_type": "join",
            "width": width_value.to_json(),
        })

    passed = all(cell_passes) and joins_pass and evaluation_count <= limit
    previous = _append_record(records, previous, {
        "candidate_index": candidate_index,
        "cells_attempted": len(cells),
        "cells_passed": sum(cell_passes),
        "evaluation_count": evaluation_count,
        "joins_passed": joins_pass,
        "passed": passed,
        "record_type": "candidate_end",
    })
    return passed, previous, {
        "candidate_index": candidate_index,
        "lambda_width": width.to_json(),
        "tube_radius": radius.to_json(),
    }


__all__ = [name for name in globals() if not name.startswith("__")]
