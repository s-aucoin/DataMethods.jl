using DimensionalData
using Rasters
import Base.Threads.@threads

export find_ctd_profiles, correct_CT_lag, remove_loops, FindLoops


############################################
"""
    find_ctd_profiles(pressure; sm = 65, Pthresh = 1, stdthresh = 4, nomotion_tol = 0.05)

Find down and upcasts from the `pressure` timeseries of a CTD (developed for the ecoCTD).

...
# Arguments
- `pressure`: time series of pressure.
- `sm`: smoothing factor for the pressure derivative.
- `Pthresh`: threshold of maximum pressure for the "surface".
- `stdthresh`: factor of the standard deviation of the derivative of pressure that is used to determine if the CTD is profiling.
- `nomotion_tol`: tolerance factor for how close to 0 the speed has to be to be considered stopped.
...
"""
function find_ctd_profiles(pressure; sm = 65, Pthresh = 1, stdthresh = 4, nomotion_tol = 0.05)

    ## First, calculate and smooth the derivative of pressure ##
    dPdt = [NaN; moving_average(∂_cen(pressure, 1:size(pressure)[1])[2:end-1], sm); NaN]

    ## Next, use a threshold of pressure (added to the minimum pressure) to determine values at the "surface" ##
    surface_inds = findall(pressure .< (minimum(pressure) + Pthresh)) # the indices of pressure that are at the "surface"

    ## Now, use the "surface" dPdt variability to determine a threshold on whether the ctd is profiling ##
    dPdt_downthresh = stdthresh*std(filter(!isnan, dPdt[surface_inds]))

    ## Use the threshold to get rough estimates of profiles ##
    dc_inds_rough = findall(dPdt .> dPdt_downthresh)

    ## Then separate the casts by non-consecutive indices ##
        downcast_start = [dc_inds_rough[1]; dc_inds_rough[findall(diff(dc_inds_rough) .> 1) .+ 1]] # need to add the first start manually
    downcast_end_guess = [dc_inds_rough[findall(diff(dc_inds_rough) .> 1)]; dc_inds_rough[end]]    # Similarly need to add the last end

    ## Now, loop through each profile and refine the estimates of the end position ##
    ncasts = size(downcast_start)[1]
    downcast_inds = Array{UnitRange{Int64}}(undef, ncasts)
    upcast_inds = Array{UnitRange{Int64}}(undef, ncasts)

    @threads for profile in 1:ncasts
        # The indices of the first guess of a downcast
        downcast_inds_guess = downcast_start[profile]:downcast_end_guess[profile]

        # Find the maximum pressure within the guess
        (val, maxP_ind) = findmax(pressure[downcast_inds_guess])
        downcast_end = downcast_inds_guess[maxP_ind] # -> Define this as the end of the downcast, start of the upcast

        # Refine the indices of the downcast
        downcast_inds[profile] = downcast_start[profile]:downcast_end


        ## Find the second zero crossing of dPdt after the profile bottom as a first guess of the end of the upcast ##
        # this works by assuming the two values adjacent to a zero crossing have opposite signs #
        upcast_start = downcast_end + 1
        upcast_end_guess = upcast_start + findall(dPdt[upcast_start:end-1].*dPdt[upcast_start+1:end] .<= 0)[2]
        upcast_inds_guess = upcast_start:upcast_end_guess

        # Find the minimum dPdt (fastest upward speed) within the upcast guess
        (val, mindPdt_ind) = findmin(dPdt[upcast_inds_guess])
        start_guess = upcast_inds_guess[mindPdt_ind]

        # Then define the end of the upcast as the first time after the maximum speed position dPdt comes within a threshold of 0
        nomotion_thresh = - nomotion_tol*(maximum(filter(!isnan, dPdt)) - minimum(filter(!isnan, dPdt)))
        upcast_end = start_guess + findfirst(dPdt[start_guess:upcast_end_guess] .>= nomotion_thresh)

        # Refine the indices of the upcast
        upcast_inds[profile] = (upcast_start + 1):upcast_end
    end

    return (; downcast_inds, upcast_inds)
end


"""
    correct_CT_lag(T, fs; τ = 0.35)
Shift the temperature `T` forward in time by `τ` seconds to account for the inherent lag between conductivity and temperature measurements for RBR CTDs.

The approach and default value of 0.35 seconds are presented in Dever et al. (2022).
"""
function correct_CT_lag(T, fs; τ = 0.35)
    τ_idx = Int(round(τ * fs)) # rounding to the nearest index (to avoid having to interpolate)
    T_corrected = fill!(similar(T), NaN) # initialize an array of NaNs
    T_corrected[1:end-(τ_idx-1)] = T[τ_idx:end] # fill in with the shifted temperatures

    return T_corrected
end



"""
    remove_loops(data; tol = 0, max_z = 0, replace=false)

Remove all values in the fields of `data` that have a higher `z` value than the previous lowest.

This is equivalent to removing points that do not monotonically decrease in `z`. Set `tol` to add a tolerance to the strict `z < z_min`.
Optionally set a starting `max_z` to something other than `0`.

If `replace=true` then the the data points are replaced with `NaN` instead of being removed.
"""
function remove_loops(data::DimStack; tol = 0, max_z = 0, replace=false)
    data = deepcopy(data)
    idx2replace = Int64[] # initalize the array of indices to replace

    # Ensure the data has a z variable
    if !hasproperty(data, :z)
        error("Input data must have property z")
    end

    # loop over each time and flag values that aren't monotonically decreasing in z
    for tidx in 1:size(data)[1]
        if data.z[tidx] < (max_z .+ tol)
            max_z = data.z[tidx]
        else
            push!(idx2replace, tidx)
        end
    end

    if replace
        # replace the flagged values with NaNs
        for val in keys(data)
            data[val][idx2replace] .= NaN
        end
        return data
    else
        return data[Not(idx2replace)]
    end

end


function remove_loops(data::RasterStack; tol = 0, max_z = 0, replace=false)
    data = deepcopy(data)
    idx2replace = Int64[] # initalize the array of indices to replace

    # Ensure the data has a z variable
    if !hasproperty(data, :z)
        error("Input data must have property z")
    end

    # loop over each time and flag values that aren't monotonically decreasing in z
    for tidx in 1:size(data)[1]
        if data.z[tidx] < (max_z .+ tol)
            max_z = data.z[tidx]
        else
            push!(idx2replace, tidx)
        end
    end

    if replace
        # replace the flagged values with NaNs
        for val in keys(data)
            data[val][idx2replace] .= NaN
        end
        return data
    else
        return data[Not(idx2replace)]
    end

end


function remove_loops(data::NamedTuple; tol = 0, max_z = 0, replace=false)
    data = deepcopy(data)
    idx2replace = Int64[] # initalize the array of indices to replace

    # Ensure the data has a z variable
    if !hasproperty(data, :z)
        error("Input data must have property z")
    end

    # loop over each time and flag values that aren't monotonically decreasing in z
    for tidx in 1:size(data.z)[1]
        if data.z[tidx] < (max_z .+ tol)
            max_z = data.z[tidx]
        else
            push!(idx2replace, tidx)
        end
    end

    if replace
        # replace the flagged values with NaNs
        for val in keys(data)
            data[val][idx2replace] .= NaN
        end
        return data
    else
        return data[Not(idx2replace)]
    end
end


"""
    FindLoops(p; tol = 0, max_z = 0)

Find all indices in `p` that are not monotonically increasing.

Optionally specify a tolerance `tol` and a starting `max_z`.
"""
function FindLoops(p; tol = 0, max_z = 0)

    # loop over each time and flag values that aren't monotonically decreasing in z
    idx2replace = Int64[] # initalize the array of indices to replace
    for tidx in eachindex(p)
        if p[tidx] > (max_z .+ tol)
            max_z = p[tidx]
        else
            push!(idx2replace, tidx)
        end
    end
    return idx2replace
end