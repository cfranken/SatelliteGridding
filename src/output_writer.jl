# Time-slice output shared by rectangular, cubed-sphere, and center gridding.
function _masked_output_slice(data, weights, nc_var; sigdigits,
                               min_weight=1e-10, normalize=false)
    values = normalize ? data ./ max.(weights, eps(eltype(data))) : data
    result = round.(values; sigdigits=sigdigits)
    fill_value = _output_fill(nc_var)
    @inbounds for i in eachindex(result, weights)
        if weights[i] < min_weight || !isfinite(result[i])
            result[i] = fill_value
        end
    end
    result
end

function _write_output_slice!(nc_var, data, ct, shape)
    spatial = ntuple(_ -> Colon(), length(shape))
    # Unfold cubed-sphere accumulators only at the I/O boundary. Rectangular
    # arrays have the same shape in memory and on disk.
    nc_var[spatial..., ct] = reshape(data, shape)
    nothing
end

function _write_time_coordinates!(ds, ct, date; window_start=nothing, window_end=nothing)
    ds["time"][ct] = date
    window_start === nothing || (ds["window_start"][ct] = window_start)
    window_end === nothing || (ds["window_end"][ct] = window_end)
    nothing
end

function _write_time_slice!(ds_out, nc_vars, config, spec::AbstractGridSpec,
                            grid_data, grid_std, grid_weights, compute_std, ct, date;
                            window_start=nothing, window_end=nothing)
    shape = _output_spatial_shape(spec)
    _write_output_slice!(nc_vars["n"], grid_weights, ct, shape)
    _write_time_coordinates!(ds_out, ct, date; window_start=window_start, window_end=window_end)
    maximum(grid_weights) > 0 || return nothing
    for (column, key) in enumerate(keys(config.grid_vars))
        data = _masked_output_slice(@view(grid_data[:, :, column]), grid_weights,
                                     nc_vars[key]; sigdigits=8)
        _write_output_slice!(nc_vars[key], data, ct, shape)
        if compute_std
            std_var = nc_vars[key * "_std"]
            std = _masked_output_slice(@view(grid_std[:, :, column]), grid_weights,
                                        std_var; sigdigits=6)
            _write_output_slice!(std_var, std, ct, shape)
        end
    end
    nothing
end

function _emit_l2_snapshot!(ds, nc_vars, zg, config, spec, data, std, weights,
                            compute_std, ct, date; window_start, window_end, zarr_epoch)
    ds === nothing || _write_time_slice!(ds, nc_vars, config, spec, data, std,
        weights, compute_std, ct, date; window_start=window_start, window_end=window_end)
    zg === nothing || write_zarr_day!(zg, date, data, std, weights, config.grid_vars,
                                      compute_std; epoch=zarr_epoch)
    nothing
end

function _write_center_time_slice!(ds_out, nc_vars, config, grid_data,
                                   grid_weights, veg_indices::Bool,
                                   ct::Int, d::DateTime, min_count)
    ds_out["time"][ct] = d
    nc_vars["n"][:, :, ct] = grid_weights

    co = 1
    for (key, _) in config.grid_vars
        _write_center_var!(nc_vars[key], grid_data[:, :, co],
                           grid_weights, ct, min_count)
        co += 1
    end

    if veg_indices
        for key in ("EVI", "NDVI", "NIRv", "NDWI")
            _write_center_var!(nc_vars[key], grid_data[:, :, co],
                               grid_weights, ct, min_count)
            co += 1
        end
    end
    nothing
end

function _write_center_var!(nc_var, data, weights, ct::Int, min_count)
    da = _masked_output_slice(data, weights, nc_var; sigdigits=5,
                              min_weight=min_count, normalize=true)
    nc_var[:, :, ct] = da
    nothing
end
