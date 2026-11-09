module DataMethods

using Statistics
using ExtraStats
using PyCall

export IntervalAverage, interp_2d

include("CalculusUtilities.jl")
include("SignalProcessing.jl")
include("adcp_utilities/ADCPUtilities.jl")
include("ctd_utilities/CTDUtilities.jl")
include("gps_utilities/GPSUtilities.jl")



"""
    IntervalAverage(data::AbstractMatrix{T} where T<:Any, Δx, fs, D)

Average `data` in chunks of physical dimensional size `Δx` over its second dimension.
`fs` is the sample frequency of the second dimension, and has inverse units of `Δx`.
`D` is the fraction of the interval `Δx` that actually exists in `data` (the duty cycle fraction).

The purpose of this function is to take a dataset that has discontinous jumps at intervals of size `D*Δx`,
add implied "no data" of length `(1-D)*Δx`, then average over interval `Δx`.
This function replaces lower_tres_mean.
"""
function IntervalAverage(data::AbstractMatrix{T} where T<:Any, Δx, fs, D)

    Nx_chunk = Int(fs*Δx*D) # number of indices that constitues one chunk in x to be averaged
    new_Nx = Int(length(data)/(size(data, 1)*Nx_chunk)) # The new number of total indices in the x dimension

    data_reshaped = reshape(data, (size(data, 1), Nx_chunk, new_Nx)) # Reshape the data into a 3D array where the second dimension corresponds to the chunk size

    new_x = ((1:new_Nx).*Δx) .- Δx/2 # New x coordinates corresponding to the averaged data in units of Δx

    return (data = nanmean(data_reshaped; dims=2), x = new_x)

end

"""
    IntervalAverage(data::AbstractMatrix{T} where T<:Any, Δx, fs, D, δ₁, δ₂)

Average `data` in chunks of physical dimensional size `Δx` over its second dimension.
`fs` is the sample frequency of the second dimension, and has inverse units of `Δx`.
`D` is the fraction of the interval `Δx` that actually exists in `data` (the duty cycle fraction).
`δ₁` and `δ₂` are fractions of `D` for which there is data for the first and last chunks respectively.
This is to account for situations where the start or end don't have a complete duty cycle.

The purpose of this function is to take a dataset that has discontinous jumps at intervals of size `D*Δx`,
add implied "no data" of length `(1-D)*Δx`, then average over interval `Δx`.
This function replaces lower_tres_mean.
"""
function IntervalAverage(data::AbstractMatrix{T} where T<:Any, Δx, fs, D, δ₁, δ₂)

    Nx_chunk = Int(fs*Δx*D) # number of indices that constitues one chunk in x to be averaged
    Nx_first_chunk = Int(δ₁*fs*Δx*D)
    Nx_last_chunk = Int(δ₂*fs*Δx*D)

    first_chunk = data[:,1:Nx_first_chunk]
    last_chunk = data[:,end-Nx_last_chunk+1:end]
    middle_chunks = data[:,Nx_first_chunk+1:end-Nx_last_chunk]

    new_middle_Nx = Int(length(middle_chunks)/(size(middle_chunks, 1)*Nx_chunk)) # The new number of total indices in the x dimension

    data_reshaped = reshape(middle_chunks, (size(middle_chunks, 1), Nx_chunk, new_middle_Nx)) # Reshape the data into a 3D array where the second dimension corresponds to the chunk size

    first_chunk_meaned = nanmean(first_chunk; dims=2)
    last_chunk_meaned = nanmean(last_chunk; dims=2)
    data_meaned = hcat(first_chunk_meaned, nanmean(data_reshaped; dims=2), last_chunk_meaned)

    new_x = ((1:(new_middle_Nx+2)).*Δx) .- Δx/2 # New x coordinates corresponding to the averaged data in units of Δx

    return (data = data_meaned, x = new_x)

end


function IntervalAverage(data::AbstractVector{T} where T<:Any, Δx, fs, D)

    Nx_chunk = Int(fs*Δx*D) # number of indices that constitues one chunk in x to be averaged
    new_Nx = Int(length(data)/Nx_chunk) # The new number of total indices in the x dimension

    data_reshaped = reshape(data, (Nx_chunk, new_Nx))

    new_x = ((1:new_Nx).*Δx) .- Δx/2 # New x coordinates corresponding to the averaged data in units of Δx

    return (data = nanmean(data_reshaped; dims=1), x = new_x)

end


function IntervalAverage(data::AbstractVector{T} where T<:Any, Δx, fs, D, δ₁, δ₂)

    Nx_chunk = Int(fs*Δx*D) # number of indices that constitues one chunk in x to be averaged
    Nx_first_chunk = Int(δ₁*fs*Δx*D)
    Nx_last_chunk = Int(δ₂*fs*Δx*D)

    first_chunk = data[1:Nx_first_chunk]
    last_chunk = data[end-Nx_last_chunk+1:end]
    middle_chunks = data[Nx_first_chunk+1:end-Nx_last_chunk]

    new_middle_Nx = Int(length(middle_chunks)/Nx_chunk) # The new number of total indices in the x dimension

    data_reshaped = reshape(middle_chunks, (Nx_chunk, new_middle_Nx)) # Reshape the data into a 3D array where the second dimension corresponds to the chunk size

    first_chunk_meaned = nanmean(first_chunk)
    last_chunk_meaned = nanmean(last_chunk)
    data_meaned = vcat(first_chunk_meaned, nanmean(data_reshaped; dims=1), last_chunk_meaned)

    new_x = ((1:(new_middle_Nx+2)).*Δx) .- Δx/2 # New x coordinates corresponding to the averaged data in units of Δx

    return (data = data_meaned, x = new_x)

end



"""
    interp_2d(arrayjl)

Interpolate `array` in 2 dimensions using a 2D cubic spline fit.

Uses scipy to make my life easier.
"""
function interp_2d(arrayjl)
    # Use python to do this easily
    py"""
    import numpy as np
    from scipy import interpolate  

    # Import the julia array into python
    array = $arrayjl

    # Get a mask of NaNs
    array = np.ma.masked_invalid(array)

    # Make vectors of the indices of x and y of the array
    h, w = array.shape[:2]

    # Turn them into grids
    xx, yy = np.meshgrid(np.arange(w), np.arange(h))

    # Extract the known good values
    known_x = xx[~array.mask]
    known_y = yy[~array.mask]
    known_v = array[~array.mask]

    # Extract the missing values
    missing_x = xx[array.mask]
    missing_y = yy[array.mask]

    interp_values = interpolate.griddata(
        (known_x, known_y), known_v, (missing_x, missing_y),
        method='cubic')

    interp_array = array.copy()
    interp_array[missing_y, missing_x] = interp_values
    """

    return py"interp_array"

end

end
