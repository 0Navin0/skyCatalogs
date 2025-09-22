import os
from pathlib import Path
import pandas as pd
import healpy

PACKAGE_DIR = os.path.dirname(os.path.abspath(str(Path(__file__).parent.parent)))
DATA_DIR = os.path.join(PACKAGE_DIR, "skycatalogs", "data")


# Used only by creator
def get_trilegal_hp_nrows(hp, nside=32):
    counts_path = os.path.join(DATA_DIR, "trilegal", "star_counts.parquet")
    tbl = pd.read_parquet(counts_path)
    nrows = (tbl.query("hp_ring_id == @hp")["hp_star_count"]).values[0]
    return nrows


_NSIDE = 32
_TRILEGAL_RING_NSIDE = 256
_TRILEGAL_NEST_NSIDE = 4096


# Used by creator and (indirectly) by API
def find_trilegal_subpixels(hp, n_rows, in_nside=32, n_gps=None,
                            max_query_rows=1_000_000):
    '''
    Return subhealpixels belonging to specified healpixel for all row groups

    Parameters
    ----------
    hp             int      healpixel (ring naming) id
    in_nside       int      NSIDE value for original tiling (normally 32)
    max_query_rows int      Used to determine how many row groups will be
                            written

    Returns
    -------
    out_nside  int      32,  256 or 4096
    out_ring   boolean  True if out list ids use ring naming; false for nest
    out_list  hp ids for each row group.

    So, for example, if there are 4 row groups, out_list will have 4 elememts,
    each itself a list.   If there is only to be one row group, out_list will
    be   [ [hp] ]

    '''
    # For almost all cases
    out_nside = _TRILEGAL_RING_NSIDE
    out_ring = True

    if n_rows <= max_query_rows:
        n_group = 1
    else:
        # break up into several queries
        if n_rows > 60 * max_query_rows:
            n_group = 64
        elif n_rows > 30 * max_query_rows:
            n_group = 32
        elif n_rows > 16 * max_query_rows:
            n_group = 16
        else:
            n_group = 4

    if hp == 9246:
        # This healpixel, even though it is not the healpixel with the most
        # stars, apparently includes a very dense section so that 64
        # subqueries is not fine enough.  One times out.
        n_group = 256
        out_nside = _TRILEGAL_NEST_NSIDE
        out_ring = False
    elif hp == 9119:
        n_group = 64

    if n_gps:            # Use it if it was passed in
        n_group = n_gps

    def _next_level(pixel):
        return [4 * pixel, 4 * pixel + 1, 4 * pixel + 2, 4 * pixel + 3]

    pixels = [healpy.ring2nest(_NSIDE, hp)]
    current_nside = in_nside

    while current_nside < out_nside:
        pixels = [pix for p in pixels for pix in _next_level(p)]
        current_nside = current_nside * 2

    if n_group <= 64:
        subpixels = [healpy.nest2ring(_TRILEGAL_RING_NSIDE, p) for p in pixels]
    else:
        subpixels = pixels

    per_group = len(subpixels) // n_group
    out_list = []
    for i in range(n_group):
        out_list.append(subpixels[i * per_group: (i + 1) * per_group])

    return out_nside, out_ring, out_list


# Used only by API
def get_trilegal_active(rdr, pixel, region):
    '''
    Find subpixel groups (corresponding to row groups) which intersect the
    region.  The purpose is to allow the caller to optimize by excluding
    row groups which don't intersect. No optimization is possible if pixels
    for each row group intersect the region

    Parameters
    ----------
    rdr      ParquetReader
    pixel    id of Original (typically nside=32) healpixel
    region   skycatalogs.shapes.Region

    Returns
    -------
    List of bool, length = number of row groups.

    '''
    n_gp = rdr.n_row_groups
    if n_gp == 1:
        return [True]
    if region is None:
        return [True for i in range(n_gp)]

    # Outputs specify how the healpixel is partitioned in subpixels
    # within the row groups
    sub_nside, out_ring, subs = find_trilegal_subpixels(pixel, rdr.n_rows,
                                                        n_gps=n_gp)
    # Find all subpixels of the hp which are also in the region
    region_nsub = set(region.get_intersecting_hps(sub_nside, out_ring))

    active = [bool(region_nsub.intersection(set(sub))) for sub in subs]

    return active
