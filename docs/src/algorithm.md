# Algorithm

## Footprint Oversampling

Satellite instruments observe the atmosphere or surface through irregularly shaped
footprints. Each footprint is a quadrilateral defined by 4 corner coordinates
(latitude/longitude pairs). These quadrilaterals can overlap multiple output grid
cells, so simply assigning each observation to its center coordinate would lose
spatial information.

`SatelliteGridding.jl` uses **sub-pixel oversampling** to properly distribute each
observation across the grid cells it covers. The default path treats footprints as
quadrilaterals; `CircularFootprintGridding` treats the four coordinates as a
circular-footprint bounding box or edge points.

### Step 1: Corner Sorting (CCW)

Different satellite products store corner coordinates in different orders. TROPOMI
uses counter-clockwise (CCW) order, while OCO-2 uses a different convention. Before
any geometric computation, corners are sorted into CCW order based on their angle
from the footprint centroid:

```
        3 ──── 2          Sorted CCW:
       /      /           1 = smallest angle from centroid
      /      /            2 = next angle
     4 ──── 1             Edges 1→2 and 4→3 are opposite sides
```

The sorting uses a 5-comparator sorting network (optimal for 4 elements), which
avoids allocations and works efficiently in GPU kernels.

### Step 2: Sub-Pixel Computation

The quadrilateral is subdivided into `n×n` sub-pixels by:

1. Dividing edge 1→2 into `n` equally-spaced points
2. Dividing edge 4→3 into `n` equally-spaced points
3. Connecting corresponding points across the two edges
4. Subdividing each connecting line into `n` interior points

This produces `n²` sub-pixel positions that fill the footprint. Each sub-pixel
carries weight `1/n²`.

```
    4 ────────────── 3
    |  ·  ·  ·  ·  |     n = 4
    |  ·  ·  ·  ·  |     16 sub-pixels
    |  ·  ·  ·  ·  |     Each has weight 1/16
    |  ·  ·  ·  ·  |
    1 ────────────── 2
```

### Step 3: Accumulation

Each sub-pixel is assigned to the grid cell it falls in (via `floor` of its
fractional grid index). The observation value is accumulated into that cell
with the sub-pixel's weight.

**Fast path**: If all 4 corners fall in the same grid cell, the footprint is
entirely within one cell. The full observation is assigned with weight 1 (no
sub-pixel computation needed).

**Skip path**: If the footprint's longitude extent exceeds `n` grid cells,
it's too wide for meaningful oversampling and is skipped.

### Circular Footprints

For circular products such as GOSAT, the footprint extent can come from four
bounding coordinates or from center latitude/longitude plus a radius. If a
corner-based product does not provide explicit center latitude/longitude
variables, the center is derived from the four coordinates.

The circular path samples an `n×n` grid over the inferred bounding box and keeps
only points inside the circle. In fractional grid-index space this can become an
ellipse when `dLat != dLon`, which is the right representation for a circular
footprint on a grid whose latitude and longitude cell sizes differ.

For center-plus-radius products, set `[basic] lat` and `[basic] lon`, then either
provide a per-sounding `[basic] radius` variable or a scalar `[circle] radius`.
`[circle] radius_unit` can be `degrees`, `km`, or `m`.

Use this method from Julia:

```julia
grid(config, grid_spec, time_spec,
     CircularFootprintGridding(n_oversample=80);
     outfile="gosat_sif.nc")
```

or from the CLI:

```bash
julia --project=. bin/grid.jl l2 \
    --config examples/gosat_sif_center_radius.toml \
    --footprint circle --nOversample 80 \
    --backend cpu \
    -o gosat_sif.nc
```

This method works with the sequential backend and the KA CPU/CUDA/Metal
backends. The KA path uses a separate circular index kernel and normalizes each
footprint by the number of samples inside the inferred circle/ellipse.

### Choosing `n`

The oversampling factor `n` is auto-computed if not specified:

```
n = clamp(ceil(footprint_extent / cell_size × 3), 2, 20)
```

This aims for ~3 sub-pixels per grid cell width. Typical values:
- TROPOMI at 1° grid: n ≈ 3
- TROPOMI at 0.05° grid: n ≈ 10–15
- OCO-2 at 1° grid: n ≈ 2

## Averaging Methods

### Sequential Path (Welford)

The default sequential path uses **Welford's online algorithm** for computing the
running mean and variance in a single pass:

```
for each observation x with weight w:
    total_weight += w
    delta = x - mean
    mean += (w / total_weight) * delta
    M2 += w * delta * (x - mean)    # only if compute_std=true
```

This is numerically stable even for values with large magnitudes, as it never
computes large sums. The standard deviation is: `std = sqrt(M2 / total_weight)`.

### KA Path (Sum-Based)

The KernelAbstractions path uses **weighted sums** for the mean, which is fully
parallelizable:

```
for each observation x with weight w:
    sum_data += w * x       # can use atomic add (GPU)
    sum_weight += w          # can use atomic add (GPU)

# After all files:
mean = sum_data / sum_weight
```

This enables parallel execution on both CPU and GPU backends. The trade-off is
slightly less numerical stability for Float32 with extreme values, but this is
negligible for typical satellite data ranges.

For standard deviation with the KA path, a **two-pass approach** is used:
1. Pass 1: Compute mean via sum-based accumulation
2. Pass 2: Re-read data, compute `sum(w × (x - mean)²)` per cell
3. `std = sqrt(sum_sq_dev / sum_weight)`

## Pipeline Architecture

```
Input files (NetCDF)
    │
    ▼
┌─────────────────────────┐
│  File Discovery         │  Date-based pattern matching (YYYY/MM/DD)
│  (filepatterns.jl)      │
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Quality Filtering      │  Apply filter rules from [filter] section
│  (filters.jl)           │
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Corner Sorting (CCW)   │  Normalize corner order
│  (oversampling.jl)      │
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Sub-Pixel Indices      │  Compute n×n grid cell indices per footprint
│  (oversampling.jl)      │
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Accumulation           │  Welford (sequential) or sum-based (KA)
│  (averaging.jl)         │
└──────────┬──────────────┘
           │
           ▼
┌─────────────────────────┐
│  Output (NetCDF)        │  Write gridded means, counts, optional std
│  (ncio.jl)              │
└──────────────────────────┘
```
