export ∂_cen, ∂_back, int_def_trap, int_trap


############################################
"""
    ∂_cen(fx, x)

A simple centered difference estimation of the derivative of `fx` with respect to axis `x`.
"""
function ∂_cen(fx::AbstractVector{T}, x::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(fx)
    @assert length(x) == n "fx and x must have the same length"

    Fx = fill(T(NaN), n)
    n < 3 && return Fx

    @inbounds @simd for i in 2:n-1
        dx = x[i+1] - x[i-1]
        Fx[i] = (fx[i+1] - fx[i-1]) / dx
    end
    return Fx
end


"""
    ∂_back(fx, x)

A simple backward finite difference estimation of the derivative of `fx` with respect to axis `x`.
"""
function ∂_back(fx::AbstractVector{T}, x::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(fx)
    @assert length(x) == n "fx and x must have the same length"

    Fx = fill(T(NaN), n)
    n < 3 && return Fx

    @inbounds @simd for i in 2:n-1
        dx = x[ii] - x[ii-1]
        Fx[i] = (fx[ii] - fx[ii-1]) / dx
    end
    return Fx
end


"""
    int_def_trap(fx, Δx)

A simple trapezoidal method for estimating a definite integral of `fx` with separation `Δx`.
"""
function int_def_trap(fx::AbstractVector{<:AbstractFloat}, Δx::Real)
        return sum(1/2 * Δx * (fx[1:end-1] + fx[2:end]))
end


"""
    int_trap(fx, x)

A simple trapezoidal method for estimating an indefinite integral of `fx` with axis `x`.
"""
function int_trap(fx::AbstractVector{T}, x::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(x)
    @assert length(fx) == n "fx and x must have the same length"

    Fx = fill(T(NaN), n)
    n == 0 && return Fx

    Fx[1] = zero(T)
    acc = zero(T)

    @inbounds @simd for i in 2:n
        acc += (x[i] - x[i-1]) * (fx[i] + fx[i-1]) / 2
        Fx[i] = acc
    end

    return Fx
end
