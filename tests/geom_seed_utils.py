import numpy as np


def rand_geom_seeds(arr, geom_key):
    """Same as np.random.rand(*arr.shape), but zeros the seeds of any CDCL polar
    ("cdcl" or "cdclsec") that is not set (all zeros).

    Perturbing an unset polar creates CL points that are out of order, which is an
    invalid polar ("* CDCL: input CL data out of order") and has no meaningful derivative.
    The same number of random values are drawn, so the random stream is unchanged.
    """
    seeds = np.random.rand(*arr.shape)
    if geom_key in ("cdcl", "cdclsec"):
        polars = seeds.reshape(-1, 6)
        unset = ~np.any(np.asarray(arr).reshape(-1, 6) != 0.0, axis=1)
        polars[unset] = 0.0
        seeds = polars.reshape(arr.shape)
    return seeds
