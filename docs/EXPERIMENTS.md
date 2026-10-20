# Experiments and archived files

All instance data in this repository were generated for this study. No
proprietary operational records or external benchmark data are included.

## Algorithmic comparison and scalability

`data/benchmark` contains the 320 formal instances and their original
specifications. The first 160 have 20 or 30 jobs and are shared by the algorithmic
comparison. The remaining 160 have 40 or 50 jobs.

The size–fleet pairs are (20,2), (20,3), (30,2), (30,3), (40,3), (40,5), (50,3),
and (50,5). Each pair is crossed with Gamma in {1,3}, tardiness factor beta in
{0.30,0.50}, and ten replications. Nominal processing times and energy consumption
are independent discrete uniforms on {1,...,20}. The deviation level is 0.30,
capacity scale is four, swapping factor is two, and relative due-date range is
0.50. The exact values, including due dates, are stored in each instance.

The recorded results are in `results/recorded/lbbd_scalability_formal`,
`results/recorded/baseline_formal_milp`, and
`results/recorded/baseline_formal_cp`. LBBD results on the smaller instances are
reused in the comparison rather than counted as a separate experiment.

`scripts/summarize_results.py` produces `algorithmic_comparison.csv`,
`baseline_initialization.csv`, `scalability.csv`, and `unproved_instances.csv`
under `results/tables`. These correspond to the manuscript's algorithmic,
initialization, scalability, and unproved-instance tables. The cap of 1800
seconds is applied to the average time of unproved runs, as in the manuscript.

## Model comparison

`data/model_comparison/cases` contains 90 cases with n=30 and m=3. Gamma ranges
from zero to four in increments of 0.5, with ten replications per level. Every
case retains its instance, the nominal-optimal assignment, and the two-stage
assignment. A solution file stores job assignments, not a single universal
recourse schedule.

The recorded two-stage reference panel is in
`results/recorded/sensitivity_gamma_n30_m3`. Its generated inputs are in
`data/sensitivity_gamma_n30_m3`. The static and nominal-assignment evaluation
results are in `results/recorded/model_comparison_n30_m3`.

`results/tables/model_comparison.csv` retains lower and upper endpoints for
interrupted static solves. Percent reductions are averaged per instance, not
computed as a ratio of mean objectives. The static counterpart is proved
optimal in 86 of 90 cases; the remaining four retain certified intervals.

## Sensitivity figures

Each sensitivity directory includes ten original reference instances, the
prepared instances for every setting, and the factor specifications. Technology
and uncertainty changes preserve the reference due dates. Deviations are
rescaled from the original samples. Technology changes use the original job
data to calculate capacities and swapping durations.

| Study | Data and recorded-result directory | Runs |
| --- | --- | --- |
| Budget sweep | `sensitivity_gamma_n40` | 90 |
| Deviation sweep at Gamma=0.5 | `sensitivity_crossed_n40_gamma0.5` | 50 |
| Deviation sweep at Gamma=2 | `sensitivity_crossed_n40_gamma2` | 50 |
| Deviation sweep at Gamma=4 | `sensitivity_crossed_n40_gamma4` | 50 |
| Capacity and swapping duration | `sensitivity_technology_n40` | 250 |

The figure data are summarized in `results/tables/sensitivity.csv`. Native
PGFPlots snippets are supplied in `results/figures`; use them in a LaTeX document
that loads TikZ, PGFPlots, and its `groupplots` library. Their coordinates are
the means and individual values from the recorded runs.

The capacity scales are {3,3.5,4,4.5,5} and swapping factors are {1,1.5,2,2.5,3}.
The crossed uncertainty study uses Gamma in {0.5,2,4} and deviation levels in
{0.10,0.20,0.30,0.40,0.50}. There are no additional unreported capacity levels
in this archive.
