# Do not interpret the current CSV as a ternary prediction result

`ob1a12_ternary_comparison.csv` was produced as a **software smoke test** using
the fixed-template scaffold files, not fresh ColabFold ternary predictions.
Its `/tmp/` input paths make that provenance visible. It will be overwritten by
the evaluator after the three real rank-1 ColabFold outputs are supplied.

Run:

```text
python3 src/evaluate_ob1a12_ternary_models.py /path/to/ternary_colabfold_results
```
