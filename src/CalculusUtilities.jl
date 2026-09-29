export ∂_cen, ∂_back, int_def_trap, int_trap


############################################
"""
    ∂_cen(fx, x)

A simple centered difference estimation of the derivative of `fx` with respect to axis `x`.
"""
function ∂_cen(fx, x)
    Fx = similar(fx)
    Fx[[1,end]] .= NaN
    for ii in 2:size(Fx,1)-1
        Fx[ii] = (fx[ii+1] - fx[ii-1])/((x[ii+1] - x[ii-1]))
    end
    return Fx
end


"""
    ∂_back(fx, x)

A simple backward finite difference estimation of the derivative of `fx` with respect to axis `x`.
"""
function ∂_back(fx, x)
    Fx = similar(fx)
    Fx[[1,end]] .= NaN
    for ii in 2:size(Fx,1)
        Fx[ii] = (fx[ii] - fx[ii-1])/((x[ii] - x[ii-1]))
    end
    return Fx
end


"""
    int_def_trap(fx, Δx)

A simple trapezoidal method for estimating a definite integral of `fx` with separation `Δx`.
"""
function int_def_trap(fx, Δx)
        return sum(1/2 * Δx * (fx[1:end-1] + fx[2:end]))
end


"""
    int_trap(fx, x)

A simple trapezoidal method for estimating an indefinite integral of `fx` with axis `x`.
"""
function int_trap(fx, x)
    Fx = similar(fx)
    Fx[1] = 0
    for ii in 2:length(x)
        Fx[ii] = 1/2 * (x[ii] - x[ii-1]) * (fx[ii] + fx[ii-1])
    end
    return cumsum(Fx)
end