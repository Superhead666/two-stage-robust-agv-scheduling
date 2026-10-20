# Two-stage robust AGV scheduling

Source code, generated instances, and recorded computational results accompanying
**Two-Stage Robust Scheduling of Heterogeneous Automated Guided Vehicles with
Budgeted Energy Uncertainty**, by Renjie Yu and Daniel Oron.

The code implements logic-based Benders decomposition (LBBD), a hybrid iterated
greedy (HIG) warm start, and mixed-integer linear programming and constraint
programming baselines based on column-and-constraint generation (CCG). All
methods use the same exact fixed-assignment verifier over the continuous
budgeted uncertainty set. The repository also includes the static robust model
and evaluation of assignments selected using nominal energy consumption.

## Cite

Please cite the accompanying manuscript and identify the release used. No journal
or archival DOI has been assigned to this repository.

```bibtex
@misc{yu2026robustagvcode,
  author = {Yu, Renjie and Oron, Daniel},
  title = {Two-Stage Robust Scheduling of Heterogeneous Automated Guided Vehicles
           with Budgeted Energy Uncertainty: Code and Data},
  year = {2026},
  version = {1.0.0},
  url = {https://github.com/Superhead666/two-stage-robust-agv-scheduling}
}
```

Machine-readable citation information is in [CITATION.cff](CITATION.cff).

## Repository contents

| Directory or file | Contents |
| --- | --- |
| `src/common` | Instance handling, deterministic scheduling, and robust verification |
| `src/lbbd` | Assignment master, DFF separation, and scheduling conflict cuts |
| `src/hig` | HIG construction, evaluation, local search, and perturbation |
| `src/ccg` | MILP and CP scenario masters and the CCG procedure |
| `src/runner` | Timed execution, initialization, and comparison models |
| `data` | Exact generated instances and factor specifications |
| `results/recorded` | Unmodified result CSV files and saved assignments |
| `results/tables` | Table data regenerated from the recorded results |
| `results/figures` | Native PGFPlots sources for the sensitivity figures |
| `scripts` | Experiment execution, archive checks, and postprocessing |
| `tests` | Small correctness and integration checks |
| `docs` | Data format, experiment mapping, and reproducibility information |
| `SHA256SUMS` | SHA-256 checksums for the release files |

The package excludes solver distributions, license files, compiled executables,
development logs, exploratory experiment folders, and obsolete experiment
variants. Explanatory and licensing comments remain in source files.

## Requirements

- A C++17 compiler and GNU Make on macOS or Linux. The process runner uses POSIX
  process control and is not supported natively on Windows.
- Gurobi Optimizer 13.0.2, including its C++ library and a valid license.
- IBM ILOG CPLEX Optimization Studio 22.1.1, including CP Optimizer and Concert,
  for CP–CCG only.
- Python 3.10 or later. The Python scripts use only the standard library.
- A LaTeX installation with TikZ and PGFPlots to use the figure snippets. It is
  not required to build the solvers or regenerate table data.

Commercial solvers must be installed and licensed separately. They are not
covered by this repository's license. The original experiments used one solver
thread on an Apple M2 Max Mac Studio with 32 GB RAM and macOS 14.0.

## Build

The Makefile defaults match the original macOS installations.

```bash
make -j4 all
make -j4 cp
```

Set installation paths when they differ. For example, on Linux:

```bash
make -j4 all \
  GUROBI_ROOT=/opt/gurobi1302/linux64 GUROBI_LIB=gurobi130
make -j4 cp \
  GUROBI_ROOT=/opt/gurobi1302/linux64 GUROBI_LIB=gurobi130 \
  CPLEX_ROOT=/opt/ibm/ILOG/CPLEX_Studio2211 CPLEX_PLATFORM=x86-64_linux
```

Compilation on Linux is supported by the build settings but has not been
validated in this release. The clean build and tests were run on arm64 macOS.
Use the same overrides for validation targets. A C++17-compatible Gurobi C++
library is required. If necessary, rebuild that library using the instructions
supplied with Gurobi.

## Check the installation and archive

```bash
make check
make check-cp
make check-data
python3 scripts/summarize_results.py
```

`make check` tests bounds, schedule feasibility, continuous-uncertainty coverage,
conflict cuts, HIG screening and verification, and the experiment drivers.
`make check-cp` additionally tests the CP master and six end-to-end solver
configurations. These are small validation tests, not a rerun of the formal
experiments. `make check-data` and table regeneration require no solver.

The recorded outcomes include 160/160 optimal LBBD–HIG results on the baseline
comparison panel, compared with 54/160 for MILP–CCG–HIG and 71/160 for CP–CCG–HIG.
LBBD–HIG proves all 80 40-job instances and 72 of the 80 50-job instances. These
counts are checked against the saved result files.

## Reproduce the computational study

Commands below use the archived input files. They preserve the original random
data, due dates, and factor settings. All commands are run from the repository
root. See [docs/EXPERIMENTS.md](docs/EXPERIMENTS.md) for the mapping to manuscript
tables and figures.

### Algorithmic comparison

The manuscript compares five configurations on 160 instances with 20 or 30
jobs. This command runs that panel. `lbbd-basic` is also available for additional
experiments but was not part of the reported five-configuration comparison.

```bash
python3 scripts/reproduce.py benchmark \
  --jobs 20 30 \
  --methods lbbd-hig milp-ccg-basic milp-ccg-hig cp-basic cp-hig \
  --total-limit 1800 --hig-limit 60 --seed 1 \
  --output runs/algorithmic_comparison
```

### Scalability

```bash
python3 scripts/reproduce.py benchmark \
  --jobs 40 50 --methods lbbd-hig \
  --total-limit 1800 --hig-limit 60 --seed 1 \
  --output runs/scalability
```

### Value of the modelling decisions

This uses the 90 archived 30-job, three-AGV cases and their saved nominal and
two-stage assignments. Reusing the recorded nominal optimum preserves the
tie-breaking rule in the reported comparison. Imported two-stage bounds also
remain explicit. These runs compare model values, not independent method times.

```bash
python3 scripts/reproduce.py models \
  --models static-robust nominal-assignment \
  --total-limit 1800 --seed 1 \
  --output runs/model_comparison
```

### Sensitivity analyses

```bash
python3 scripts/reproduce.py sensitivity --study gamma \
  --total-limit 1800 --hig-limit 60 --seed 1 \
  --output runs/sensitivity_gamma

python3 scripts/reproduce.py sensitivity --study crossed \
  --total-limit 1800 --hig-limit 60 --seed 1 \
  --output runs/sensitivity_crossed

python3 scripts/reproduce.py sensitivity --study technology \
  --total-limit 1800 --hig-limit 60 --seed 1 \
  --output runs/sensitivity_technology
```

The three commands cover 90 budget runs, 150 budget–deviation runs, and 250
capacity–swapping-duration runs. Each setting has ten paired replications.

## License and support

The source code, documentation, and generated data in this release are
distributed under the [MIT License](LICENSE). The authors are listed in
[AUTHORS](AUTHORS). Please report reproducible issues through this repository's
GitHub issue tracker.
