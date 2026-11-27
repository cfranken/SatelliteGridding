_new_overlap(allocate, A, ::Tuple{}) = nothing
_new_overlap(allocate, A, names) =
    NamedTuple{names}(map(name -> allocate(_overlap_type(A, Val(name))), names))

_overlap_type(A, ::Val{:sum_w2}) = A
_overlap_type(A, ::Val{:n_soundings}) = Int32

const _SOUNDING_COUNT_OVERFLOW = typemin(Int32)

# Keep an overflow sticky so parallel updates cannot wrap back to a valid count.
@inline _add_sounding_counts(a::Int32, b::Int32) =
    0 <= a && 0 <= b <= typemax(Int32) - a ? a + b : _SOUNDING_COUNT_OVERFLOW

function _check_sounding_counts(counts)
    any(<(0), counts) &&
        throw(OverflowError("n_soundings exceeds the Int32 limit ($(typemax(Int32)))"))
    nothing
end

_map_overlap(f, ::Nothing) = nothing
_map_overlap(f, stats::NamedTuple) = map(f, stats)

_host_overlap_type(::Nothing) = Nothing
_host_overlap_type(stats::NamedTuple) =
    NamedTuple{keys(stats),Tuple{map(a -> Matrix{eltype(a)}, values(stats))...}}

@inline _record_overlap!(::Nothing, cell, w, atomic) = nothing
@inline function _record_overlap!(stats::NamedTuple, cell, w, atomic)
    _add_overlap!(get(stats, :sum_w2, nothing), cell, w * w, atomic)
    _increment_soundings!(get(stats, :n_soundings, nothing), cell, atomic)
    nothing
end

@inline _add_overlap!(::Nothing, cell, x, atomic) = nothing
@inline _add_overlap!(a, cell, x, atomic) = _add_ka_cell!(a, cell, x, atomic)

@inline _increment_soundings!(::Nothing, cell, atomic) = nothing
@inline function _increment_soundings!(a::AbstractArray, cell, ::Val{false})
    @inbounds a[cell] = _add_sounding_counts(a[cell], Int32(1))
    nothing
end
@inline function _increment_soundings!(a::AbstractArray, cell, ::Val{true})
    count = @inbounds @atomic a[cell] += Int32(1)
    # If a negative count wraps to zero, restore its overflow bit. The
    # increment that wraps is responsible for marking the cell again.
    if iszero(count)
        @inbounds @atomic a[cell] |= _SOUNDING_COUNT_OVERFLOW
    end
    nothing
end

_merge_overlap!(::Nothing, ::Nothing) = nothing
function _merge_overlap!(window::NamedTuple, day::NamedTuple)
    for (name, out) in pairs(window)
        _merge_overlap!(Val(name), out, day[name])
    end
    nothing
end

_merge_overlap!(::Val{:sum_w2}, out, input) = (out .+= input)
_merge_overlap!(::Val{:n_soundings}, out, input) = (out .= _add_sounding_counts.(out, input))
