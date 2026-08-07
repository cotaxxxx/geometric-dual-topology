# Run the item 6 audits

```bash
python -m pip install sympy==1.14.0 mpmath==1.3.0
cd CERTIFICATES/prolate/item6_axis
python prolate_axis_symbolic_audit.py
python prolate_axis_center_symbolic_audit.py
python prolate_axis_tail_symbolic_audit.py
python prolate_axis_reference.py --dps 50
python prolate_axis_center_reference.py --dps 40
python prolate_axis_tail_reference.py --dps 40
python -m py_compile \
  prolate_axis_symbolic_audit.py \
  prolate_axis_center_symbolic_audit.py \
  prolate_axis_tail_symbolic_audit.py \
  prolate_axis_reference.py \
  prolate_axis_center_reference.py \
  prolate_axis_tail_reference.py
```

Expected status:

- general symbolic audit: `PASSED`
- center Hessian symbolic audit: `PASSED`
- tail coefficient symbolic audit: `PASSED`
- reference scout: `NON_CERTIFIED_REFERENCE`
- center Hessian samples: positive, but non-certified
- seven tail slopes: within the recorded floating-point tolerance
- all sampled `psi` values: positive

The mpmath results are regression evidence only. The symbolic tail audit fixes the logarithmic coefficient but does not bound the remainder. None of these outputs is an Arb positivity certificate.
