using FFTW
using DSP
import Base.Threads.@threads

export Welch, FT_params, detrend2d, dfiltfilt, findspikes, despike

############################################
"""
    Welch(segments, N; win=nothing)

Compute the power spectrum using Welch's method and periodogram of the input time series segments.


...
# Arguments
- `segments`: must be a tuple of vectors of the segments to FFT, even if there is only one segment.
- `N`: the number of elements in each segment.
- `win`: optional window to apply to each segment.
...
"""
function Welch(segments, N; win=nothing)

    # If the window is unspecified, a rectangular window is used (no window)
    if win == nothing
        win = rect(N) # Make the window rectangular
    end

    F_segments = Matrix{Any}(undef, N, size(segments)[1]) # Initialize a matrix to append to
    @threads for idx in axes(segments, 1)                 # Iterate over each segment
        win_seg = win .* segments[idx]                    # Apply the window to the segment
        F = fftshift(fft((win_seg)/sqrt(N)))              # FFT the normalized (factor of 1/√N) windowed segment
        F_segments[:, idx] = F                            # Put the transformed segment in the matrix
    end
    I_mean = vec(mean(abs.(F_segments).^2, dims=2)) # Compute the mean over the segments for each frequency

    return I_mean, F_segments

end



"""
    FT_params(fs, s_seg; fil_frac = 0, win_s = :rect, ovlap = 0.5)

Create the parameters for computing and plotting a Welch periodogram.

...
# Arguments
- `fs`: Sample frequency of the data.
- `s_seg`: Welch window length (in inverse units of `fs`).
- `fil_frac`: Smoothing window size in % of data.
- `win_s`: Type of Welch taper window function.
- `ovlap`: Fractional amount of overlap between Welch segments.
...
"""
function FT_params(fs, s_seg; fil_frac = 0, win_s = :rect, ovlap = 0.5)

    N = Int(round(s_seg * fs))                        # Number of points in a segment (must be an integer)
    overlap = ceil(N*ovlap)                    # Amount of overlap between the segments

    win_func = getfield(DSP, win_s)
    win = win_func(N)                          # Window function

    freqs = fftshift(fftfreq(N, fs))           # (Hz) Nyquist-adjusted frequencies

    fil_w = Int(2*round((fil_frac * N)/2) + 1) # Moving average filter width (always an odd integer)
    rm_p = Int(floor(fil_w/2))                 # Number of points to remove on either end from smoothing
    freq_sm_full = freqs[1+rm_p:end-rm_p]      # Adjusted freqencies for the smoothing

    pos_idx = freq_sm_full .>= 0               # Indices of the positive frequencies
    freq_sm = freq_sm_full[pos_idx]            # Only the positive frequencies
    f_max = fs/2                               # Maximum display frequency on plots
    f_min = freq_sm[2]                         # Minimum display frequency on plots

    lims = (f_min, f_max, nothing, nothing)

    return (; N, overlap, win, freqs, fil_w, freq_sm, pos_idx, f_max, f_min, lims)
end


"""
    detrend2d(data, axis)

Removes the linear best fit line from 1d or 2d `data` along the first dimension, `axis`.

Returns both the detrended data `detr` and the fitted line `fit_ln`.
"""
function detrend2d(data, axis)

    # Transform 1d data #
    if length(size(data)) == 1
        data = data'
    end

    ## Model for least squares fitting ##
    model(t, p) = p[1] * t .+ p[2]                    # p is a vector of parameters
    p0 = [0.1, 0.1]                                   # First guess of parameters

    fit_ln = Array{eltype(data)}(undef, size(data))        # Initilize a dummy matrix to fill
    @threads for ii in 1:size(data)[1]                # Loop over the first dimension of the data
        v_temp = data[ii, :]                          # Extract the time series at the range
        nonans = findall(!isnan, v_temp)              # Find all the non-NaNs
        lfit = curve_fit(model, axis[nonans], v_temp[nonans], p0) # Fit a line to it
        fit_ln[ii,:] = model(axis, lfit.param)        # Put the trend vector into the matrix
    end

    detr = data - fit_ln                              # Remove the trend

    return (detr, fit_ln)                             # Return the detrended data and the fits

end


"""
    dfiltfilt(data, axis, filters, designmethod, fs)

Filters time series `data` of one or more dimensions along `axis` with one or more `filters` using `designmethod`.

Also re-adds the linear fit after filtering.
"""
function dfiltfilt(data, axis, filters, designmethod, fs)

    (data_detrend, fitp) = detrend2d(data, axis) # Remove the linear trend

    filt_data = copy(data_detrend)
        for ii in filters
            # applies over the 1st dimension so need to prime twice to get the proper matrix
            filt_data = filtfilt(digitalfilter(ii, designmethod; fs=fs), filt_data')'
        end

    retrend_data = filt_data + fitp              # Re-add the trend to retain mean info

    # Make 1d data a vector #
    if length(size(data)) == 1
        return vec(retrend_data)
    else
        return retrend_data
    end
end

############################################

"""
    findspikes(data; sm=35, std_thresh=2)

A simple method to find spikes in data.

The moving average of `data` is computed as a reference `ref_data`.
Values in `data` that are more than `std_thresh` standard deviations 
of `data` different from `ref_data` are replaced with the 
corresponding `ref_data` point.
"""
function findspikes(data; sm=35, std_thresh=2)

    ref_data = moving_average(data, sm) # Make a reference series

    norm = data - ref_data

    return abs.(norm) .>= std_thresh*nanstd(data)
end


"""
    despike(data; sm=35, std_thresh=2)

A simple method to despike data.

The moving average of `data` is computed as a reference `ref_data`.
Values in `data` that are more than `std_thresh` standard deviations 
of `data` different from `ref_data` are replaced with the 
corresponding `ref_data` point.
"""
function despike(data; sm=35, std_thresh=2)

    ref_data = moving_average(data, sm) # Make a reference series

    spikes = findspikes(data; sm=sm, std_thresh=std_thresh)

    data[spikes] .= ref_data[spikes]

    return data
end