import sgs
import matplotlib.pyplot as plt

mrast = sgs.SpatialRaster('/home/jbmeyer/extdata/mraster.tif')
srast = sgs.stratify.quantiles(mrast, num_strata={"zq90": 5})

print("running strat_random with proportional allocation")

sgs.sample.strat(
    srast,
    num_samples=100,
    num_strata=5,
    allocation="prop",
    method="random",
    plot=True,
)

print("running strat_random with equal allocation")

sgs.sample.strat(
    srast,
    num_samples=100,
    num_strata=5,
    allocation = "equal",
    method="random",
    plot=True,
)

print("running strat_random with weighted allocation")

sgs.sample.strat(
    srast,
    num_samples=100,
    num_strata=5,
    allocation="manual",
    method="random",
    weights=[0.1, 0.1, 0.1, 0.1, 0.6],
    plot=True,
)

srast = sgs.stratify.quantiles(mrast, num_strata={"zq90": 7})

print("running strat_random with proportional allocation")

sgs.sample.strat(
    srast,
    num_samples=200,
    num_strata=7,
    allocation="prop",
    method="random",
    plot=True,
)

print("running strat_random with equal allocation")

sgs.sample.strat(
    srast,
    num_samples=200,
    num_strata=7,
    allocation = "equal",
    method="random",
    plot=True
)

print("running strat_random with weighted allocation")

sgs.sample.strat(
    srast,
    num_samples=199,
    num_strata=7,
    allocation="manual",
    method="random",
    weights=[0.1, 0.15, 0.15, 0.05, 0.3, 0.05, 0.2],
    plot=True
)

