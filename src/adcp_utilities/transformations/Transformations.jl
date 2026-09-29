export rm_sidelobe_contam


include("CoordinateTransformations.jl")
include("VelocityTransformations.jl")


############################################
"""
    rm_sidelobe_contam(data, datarange, boundary)

Remove the side-lobe-contaminated measurements from `data`.

The maximum radial datarange the ADCP can measure is equal to the distance from the ADCP to the boundary.
The contaminated measurements are therefore at along-beam `datarange`s greater than or equal to `boundary`.

...
# Arguments
- `data`: the data to decontaminate.
- `datarange`: the along beam distance from the ADCP. ADCP attitude does not affect it, and it is the same for the slant beams.
- `boundary`: the distance to the visible boundary (i.e. the boundary sound is reflecting off). This is also independent of ADCP attitude.
...
"""
function rm_sidelobe_contam(data, datarange, boundary)

    # Any cell the contamination line passes through is no good
    r_edges = datarange .+ nanmean(diff(datarange, dims=1))/2 # convert the datarange vector to datarange edges

    # Find the first index of contaminated data
        contamid = vec(mapslices(findlast, (r_edges .- boundary') .< 0, dims=1)) .+ 1
    # Find the first index of non-NaN data (if the data have already been depth-shifted)
    #firstnonanid = vec(mapslices(findfirst, .!isnan.(data), dims=1))

    # Replace the contaminated data with NaNs
    nocontam_data = copy(data)
    @threads for tt in 1:size(nocontam_data)[2]
        nocontam_data[#=firstnonanid[tt] +=# contamid[tt]:end,tt] .= NaN
    end

    return nocontam_data
end