# Reproducibility scope

## Recorded results

The result CSV files under `results/recorded` are unmodified copies of the
completed computational campaigns. Their original SHA-256 checksums are in
`results/recorded/provenance.json`. The archive includes saved assignments and
the original executable fingerprints. Commercial solver installations, license
files, compiled executables, and development logs are excluded.

The source implements the final manuscript algorithms. HIG was revised after
the formal campaigns. These changes affect candidate evaluation and selective
verification, including the 0.1-second verification-call allowance. The formal
campaigns have not been rerun with this final HIG implementation. Therefore,
the archived timings and warm-start outcomes must not be described as outputs
from rebuilding this release. Reruns of the source are new computational runs
on the same archived instances. Hardware and solver versions can also change
timings and search paths.

This release does not claim an exact reconstruction of the historical
executables from source. The preserved result files permit independent
recalculation of the reported tables and figures without rerunning optimization.

## Deterministic checks

`make check-data` validates the release-file manifest, original result hashes,
benchmark metadata, sensitivity-instance fingerprints and pairing metadata,
recorded optimality counts, and comparison-model ordering.

`python3 scripts/summarize_results.py` regenerates the numerical table data.
The manuscript's deviation and gap denominators use max(1, upper_bound).
Each method's bound is used for its certified gap, while the best available
pooled bound is used for its solution deviation. Model-comparison percentage
reductions are computed per instance before averaging.

`make check` builds all non-CP programs and runs exhaustive small-instance
checks, structural checks, HIG manuscript-alignment checks, and driver tests.
`make check-cp` additionally checks the CP master and agreement among all six
initialization–method combinations. These checks establish the tested
implementation's correctness on those cases; they do not reproduce the formal
running-time comparison or establish universal numerical guarantees.

## New runs

`scripts/reproduce.py` selects the exact archived instances, records their
hashes, and runs the selected methods in a separate output directory. It does
not modify recorded results. The benchmark and sensitivity runners freeze new
executables and reject inconsistent resume requests. Model-comparison reruns
reuse the archived nominal assignment to preserve the stated rule for tied
nominal optima and credit the recorded valid two-stage lower bound explicitly.

Any newly generated or rerun results should be reported separately from the
recorded campaigns. Changes to the source or data require a new manifest and
release identifier rather than editing an archived release in place.
