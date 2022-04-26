# iBGe

Scripts to benchmark structure learning of Bayesian networks with mixed interventional and observational data

To accompany https://arxiv.org/abs/

The file iBGe.R runs simulations for a single setting, seed and network size. This may be run locally, or called from an outside loop (like run_loop.R) for example on a cluster.

The file collateResults.R collates output files from the different runs into the sims_collated directory.
