# Release validation

The release was built and checked on arm64 macOS on 21 October 2026 with
Gurobi Optimizer 13.0.2 and IBM ILOG CPLEX Optimization Studio 22.1.1.

The following commands completed successfully from the clean release directory.

```bash
make -j4 all cp
make -j4 check check-cp
python3 scripts/verify_archive.py
python3 scripts/summarize_results.py
```

The solver tests include exhaustive small-instance comparisons, packing-bound
checks, insertion-cut validity checks, HIG screening and target checks,
deterministic and robust scheduling comparisons, CP master checks, and agreement
among the six end-to-end method–initialization combinations. Driver tests cover
input freezing, paired sensitivity data, imported instances, and resume guards.

The archive-selection wrapper was also exercised on a short benchmark run and
both comparison models at Gamma=0, including completed-run resumption. Its dry
runs select 800 algorithmic-comparison runs, 160 larger-instance LBBD runs,
180 additional comparison-model runs, and 490 sensitivity runs under the
documented settings. The 90 two-stage reference solves for the model comparison
are retained separately in the recorded data.

These checks do not replace a rerun of the formal campaigns with the final HIG
implementation. That distinction is described in REPRODUCIBILITY.md.
