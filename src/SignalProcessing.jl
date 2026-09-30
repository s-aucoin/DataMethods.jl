using FFTW
using DSP
import Base.Threads.@threads

export Welch, FT_params, findspikes, despike

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


# Example of energy conservation
```jldoctest
julia> using Random, Statistics, FFTW, DSP, DataMethods

julia> Random.seed!(1234)
TaskLocalRNG()

julia> y = randn(1000)
1000-element Vector{Float64}:
  0.9706563288552144
 -0.9792184115351997
  0.9018608835940937
 -0.03280312924463938
 -0.6007922233555612
 -1.445177115286233
  2.7074239417157804
  1.5244478634355956
  0.759804020007466
 -0.8814369061964817
  ⋮
  0.7364167064633866
  0.1919443562483462
  0.764671095596252
  0.4605476610190263
 -1.4553495505882335
 -0.731679798166314
 -0.4632854431188943
  0.5112192754137074
 -1.2911234902314417

julia> freqs = fftshift(fftfreq(length(y)))
-0.5:0.001:0.499

julia> var_seg = arraysplit(y, length(y), length(y)/2)
1-element DSP.Periodograms.ArraySplit{Vector{Float64}, Float64, Nothing}:
 [0.9706563288552144, -0.9792184115351997, 0.9018608835940937, -0.03280312924463938, -0.6007922233555612, -1.445177115286233, 2.7074239417157804, 1.5244478634355956, 0.759804020007466, -0.8814369061964817  …  -1.0486579479852471, 0.7364167064633866, 0.1919443562483462, 0.764671095596252, 0.4605476610190263, -1.4553495505882335, -0.731679798166314, -0.4632854431188943, 0.5112192754137074, -1.2911234902314417]

julia> (I_X, F_X) = Welch(var_seg, length(y))
([1.675551072924107, 0.48289321721056244, 0.43312548870576206, 0.18436518959401757, 0.2748291063463766, 2.542770369452276, 0.5146733234181756, 2.691041974711656, 0.45797673562505964, 1.4163861846888959  …  0.18444673657811564, 1.416386184688897, 0.4579767356250599, 2.691041974711655, 0.5146733234181758, 2.542770369452275, 0.2748291063463768, 0.18436518959401751, 0.4331254887057618, 0.48289321721056305], Any[-1.294430791090859 - 2.7755575615628914e-17im; 0.385197238952919 + 0.5783738447696354im; … ; -0.35683102369115904 + 0.5529892487537902im; 0.3851972389529196 - 0.5783738447696355im;;])

julia> var(y)
1.1012026004695599

julia> sum(I_X) * diff(freqs)[1]
1.1005658940486087
```
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

    win_func = getfield(FourierFuncs, win_s)
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

# Example
```
julia> using FourierFuncs, Random

julia> Random.seed!(1234)
TaskLocalRNG()

julia> x = 0.1:0.1:100
0.1:0.1:100.0

julia> y = randn(1000)
1000-element Vector{Float64}:
  0.9706563288552144
  0.871497852880908
  ⋮
  0.6530055215854716
 -1.2911234902314417

julia> ns_y = despike(y; sm=7, std_thresh=1)
1000-element Vector{Float64}:
  0.9706563288552144
  0.871497852880908
  ⋮
  0.6530055215854716
 -1.2911234902314417
```
"""
function despike(data; sm=35, std_thresh=2)

    ref_data = moving_average(data, sm) # Make a reference series

    spikes = findspikes(data; sm=sm, std_thresh=std_thresh)

    data[spikes] .= ref_data[spikes]

    return data
end