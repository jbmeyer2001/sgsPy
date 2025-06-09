from . import utils
from . import calculate
from . import sample
from . import stratify

from .utils import (
    access,
    plot,
    plot_raster,
    plot_vector,
    SpatialRaster,
    SpatialVector,
    write,
)

from .calculate import (
    allocation,
    coobs,
    lhs_optimal,
    pop,
    principal_components,
    representation,
    sample_size,
)

from .sample import (
    ahels,
    balanced,
    clhs,
    existing,
    nc,
    srs,
    strat,
    sys_strat,
    systematic,
)

from .stratify import (
    breaks,
    kmeans,
    poly,
    quantiles,
    map,
)

balanced = balanced.balanced

__all__ = list(
    set(utils.__all__) |
    set(calculate.__all__) |
    set(sample.__all__) |
    set(stratify.__all__)
)
