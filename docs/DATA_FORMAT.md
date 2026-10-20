# Instance and solution formats

An `.agv` instance is a whitespace-delimited text file.

1. The first line contains `n m alpha Gamma`.
2. The second line contains the m battery capacities, in AGV index order.
3. Each of the following n lines contains the job's due date, followed by m
   triples `processing_time nominal_energy maximum_energy_deviation`.

The same normalized deviation xi_j is used for job j across all potential AGV
assignments. On AGV k its energy is nominal_energy_jk + deviation_jk * xi_j.
The global set is 0 <= xi_j <= 1 and sum_j xi_j <= Gamma. It is not a separate
budget for each block. Machine-specific verification uses the projection of
this global set, with the full budget, as established in the manuscript.

Processing times, due dates, and swapping duration are integer-valued in the
published instances. Energy values and Gamma may be fractional. Every AGV is
initially fully charged, and each swap restores its usable energy to its
capacity. A job includes the return to its assigned AGV's loading location.

A solver `.sol` file begins with `n upper_bound lower_bound`. The next n integer
tokens give the AGV assigned to each job, using zero-based indices. Optional
following lines describe the result. Recourse schedules can differ between
realizations, so this assignment file is not a single fixed operational schedule.
Static comparison results additionally provide a schedule CSV with job order,
swaps, completion times, and blocks.
