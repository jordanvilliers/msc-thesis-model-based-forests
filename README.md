# Replication Code for Causal Effect Estimation with Model-Based Forests under Model Misspecification

by Jordan Villiers, UZH

This repository includes the R code for the thesis "Causal Effect
Estimation with Model-Based Forests under Model Misspecification".

The core files are:

- `simulation.R`: runs the simulation study (sources `def.R` and
  `experiments.R` from `simulation_code/`), producing a new,
  timestamped `.rds` file in `results/`
- `reproducible.R`: produces every figure in the thesis from a saved
  `.rds` file in `results/` (sources `create_plots.R` from
  `simulation_code/`)

## Reproducing the thesis figures exactly

`results/res_thesis.rds` already contains the simulation output used
to produce every figure in the thesis. Simply run:

```r
source("reproducible.R")
```

## Running your own simulation

`simulation_code/def.R` sets `testonly <- TRUE` by default. In this
mode, running `simulation.R` only executes the first and last job of
the full design (via `batchtools::testJob`), just to confirm the
simulation pipeline runs correctly, rather than submitting the entire
study.

To run the **full** simulation study, set `testonly <- FALSE` in
`simulation_code/def.R`. Be aware that the full study amounts to
32,000 simulated datasets and was originally run on an HPC cluster;
running it on a personal machine is likely to take a very long time.

```r
source("simulation.R")     # writes a new, timestamped .rds to results/
source("reproducible.R")   # uses the most recently created .rds file
```

`reproducible.R` always uses the most recently created `.rds` file in
`results/`. If you run your own simulation, its output will
automatically take precedence over `res_thesis.rds`, since it will be
more recent; `res_thesis.rds` itself is never modified or deleted. To
go back to reproducing the thesis figures, remove or move your own
`.rds` file out of `results/`.