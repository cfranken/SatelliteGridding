# Time-slice output shared by rectangular, cubed-sphere, HEALPix, and center gridding.
_round_output(value, ::Nothing) = value
function _round_output(value, sigdigits::Int)
    rounded = round(value; sigdigits)
    # Julia's decimal rounding can return signed zero on overflow. A nonzero
    # input rounded to a positive number of significant digits cannot round to zero.
    iszero(rounded) && !iszero(value) ? oftype(value, NaN) : rounded
end

function _masked_output_slice(data, weights, fill_value::T, sigdigits, default_sigdigits;
                              min_weight=1e-10) where {T}
    result = similar(data, T)
    digits = sigdigits !== nothing && sigdigits <= (T === Float32 ? 9 : 17) ? sigdigits : nothing
    # Legacy Float32 output rounded in grid precision before narrowing.
    legacy_rounding = T === Float32 && eltype(data) !== T && digits === default_sigdigits
    @inbounds for i in eachindex(result, data, weights)
        value = T(data[i])
        if weights[i] < min_weight || !isfinite(value)
            result[i] = fill_value
            continue
        end
        rounded = legacy_rounding ? T(_round_output(data[i], digits)) : _round_output(value, digits)
        result[i] = isfinite(rounded) ? rounded : fill_value
    end
    result
end

function _write_output_slice!(nc_var, data, ct, shape)
    spatial = ntuple(_ -> Colon(), length(shape))
    # Unfold cubed-sphere accumulators or flatten HEALPix cells in NESTED
    # order only at the I/O boundary. Rectangular arrays have the same
    # shape in memory and on disk.
    nc_var[spatial..., ct] = reshape(data, shape)
    nothing
end

function _write_time_coordinates!(ds, ct, date; window_start=nothing, window_end=nothing)
    ds["time"][ct] = date
    window_start === nothing || (ds["window_start"][ct] = window_start)
    window_end === nothing || (ds["window_end"][ct] = window_end)
    nothing
end

function _write_time_slice!(ds_out, nc_vars, config, spec::AbstractGridSpec, output_options,
                            grid_data, grid_std, grid_weights, compute_std, ct, date;
                            window_start=nothing, window_end=nothing)
    shape = _output_spatial_shape(spec)
    _write_output_slice!(nc_vars["n"], grid_weights, ct, shape)
    _write_time_coordinates!(ds_out, ct, date; window_start=window_start, window_end=window_end)
    maximum(grid_weights) > 0 || return nothing
    for (column, key) in enumerate(keys(config.grid_vars))
        data = _masked_output_slice(@view(grid_data[:, :, column]), grid_weights,
                                     output_options.fill_value, output_options.significant_digits_mean,
                                     _DEFAULT_MEAN_DIGITS)
        _write_output_slice!(nc_vars[key], data, ct, shape)
        if compute_std
            std_var = nc_vars[key * "_std"]
            std = _masked_output_slice(@view(grid_std[:, :, column]), grid_weights,
                                        output_options.fill_value, output_options.significant_digits_std,
                                        _DEFAULT_STD_DIGITS)
            _write_output_slice!(std_var, std, ct, shape)
        end
    end
    nothing
end

function _emit_l2_snapshot!(ds, nc_vars, zg, config, spec::AbstractGridSpec{T}, output_options,
                            data, std, weights,
                            compute_std, ct, date; window_start, window_end, zarr_epoch) where {T}
    if output_options isa NetCDFOptions{Float32}
        # Preserve the existing rounding inputs when the grid and KA accumulator differ.
        data = _spherical_output(T, data)
        std = _spherical_output(T, std)
        weights = _spherical_output(T, weights)
    end
    ds === nothing || _write_time_slice!(ds, nc_vars, config, spec, output_options, data, std,
        weights, compute_std, ct, date; window_start=window_start, window_end=window_end)
    zg === nothing || write_zarr_day!(zg, date, data, std, weights, config.grid_vars,
                                      compute_std; epoch=zarr_epoch)
    nothing
end

function _write_center_time_slice!(ds_out, nc_vars, config, output_options, grid_data,
                                   grid_weights, veg_indices::Bool,
                                   ct::Int, d::DateTime, min_count)
    ds_out["time"][ct] = d
    nc_vars["n"][:, :, ct] = grid_weights

    co = 1
    for (key, _) in config.grid_vars
        _write_center_var!(nc_vars[key], grid_data[:, :, co],
                           grid_weights, ct, min_count, output_options)
        co += 1
    end

    if veg_indices
        for key in ("EVI", "NDVI", "NIRv", "NDWI")
            _write_center_var!(nc_vars[key], grid_data[:, :, co],
                               grid_weights, ct, min_count, output_options)
            co += 1
        end
    end
    nothing
end

function _write_center_var!(nc_var, data, weights, ct::Int, min_count, output_options)
    values = data ./ max.(weights, eps(eltype(data)))
    da = _masked_output_slice(values, weights, output_options.fill_value,
                              output_options.significant_digits_center, _DEFAULT_CENTER_DIGITS;
                              min_weight=min_count)
    nc_var[:, :, ct] = da
    nothing
end
