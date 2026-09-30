using CSV
using DataFrames
using Dates
using TimeZones

export add_ms, clean_gps_data, mean_gps_positions, clean_and_mean, index_by_time


############################################
"""
    add_ms(gps_time; fs_gps = 10)

Add milliseconds to repeated seconds values in array.
"""
function add_ms(gps_time; fs_gps = 10)
    # make a vector of each second that exists in the data #
    oneHz = collect(gps_time[1]:Second(1):gps_time[end])

    t_with_ms = Array{ZonedDateTime}(undef, size(gps_time))
    for ntime in 1:size(oneHz)[1] # loop through each time
        idx_of_time = findall(x -> x == oneHz[ntime], gps_time) # find all the values with that second
        with_ms = gps_time[idx_of_time] .+ Millisecond.(collect(0:(size(idx_of_time)[1] - 1)) .* 1000/fs_gps) # add ms to them
        t_with_ms[idx_of_time] = with_ms
    end
    return t_with_ms
end


"""
    clean_gps_data(gps_pth::String; time_format = "YYYYmmddHHMMSS", tz = tz"UTC", fs_gps = 10)

Add milliseconds to GPS time with only seconds and add `missing` for data that are missing and save as a file.
"""
function clean_gps_data(gps_pth::String; time_format = "YYYYmmddHHMMSS", tz = tz"UTC", fs_gps = 10, datatypes=[Int, Bool, String, String, String, String, Int, Float64, Int])

        # Load the data #
        drift_gps = CSV.read(gps_pth, DataFrame; types=datatypes)

        ## Read the time stamps ##
        gps_time = ZonedDateTime.(DateTime.("20" .* string.(drift_gps[:,"DATE"]) .* string.(drift_gps[:,"TIME"]), time_format), tz)

        ## Read the lat/lon positions ##
        gps_lat =  parse.(Float64, strip.(drift_gps[:,"LATITUDE N/S"], 'N'))
        gps_lon = -parse.(Float64, strip.(drift_gps[:,"LONGITUDE E/W"], 'W'))

        ## Add milliseconds to the times ##
        t_with_ms = add_ms(gps_time; fs_gps = fs_gps)
        data_with_ms = unique(DataFrame(time = t_with_ms, lat = gps_lat, lon = gps_lon), "time")

        ## Merge with a full time series to add the missing times ##
        fsHz = DataFrame(time = collect(gps_time[1]:Millisecond(1000/fs_gps):gps_time[end])) # full times
        gps_clean = sort(leftjoin(fsHz, data_with_ms, on = :time, makeunique=true), :time)

        ## write the cleaned data to a new file ##
        CSV.write(splitext(gps_pth)[1] * "_cleaned.csv", gps_clean)
end


"""
    mean_gps_positions(gps_pths)

Average the latitude and longitude of multiple GPSs and save as a file.
"""
function mean_gps_positions(gps_pths)
    # read the files in #
    gps_data = CSV.read.(gps_pths, DataFrame)

    # join them into one DataFrame #
    all_gps = outerjoin(gps_data..., on = :time, makeunique=true)

    # Extract the lat and lon
    lons = select(all_gps, r"lon")
    lats = select(all_gps, r"lat")

    # Take the mean while skipping missings
    transform!(lats, AsTable(:) => ByRow(t -> emptymissing(mean)(skipmissing(t))) => :lat_mean)
    transform!(lons, AsTable(:) => ByRow(t -> emptymissing(mean)(skipmissing(t))) => :lon_mean)

    # Make a new DataFrame #
    mean_gps = DataFrame(time = all_gps[:, "time"], lat = lats[:,"lat_mean"], lon = lons[:,"lon_mean"])

    ## write the mean data to a new file ##
    CSV.write(dirname(gps_pths[1]) * "/" * "GPS_mean.csv", mean_gps)
end


"""
    clean_and_mean(gps_pths::Vector{String}; time_format = "YYYYmmddHHMMSS", tz = tz"UTC", fs_gps = 10, datatypes=[Int, Bool, String, String, String, String, Int, Float64, Int])

Clean the GPS data up and average the latitude and longitude of multiple GPSs and save each as a file.
"""
function clean_and_mean(gps_pths::Vector{String}; time_format = "YYYYmmddHHMMSS", tz = tz"UTC", fs_gps = 10, datatypes=[Int, Bool, String, String, String, String, Int, Float64, Int])
    # Clean the data up and save as files #
    clean_gps_data.(gps_pths; time_format, tz, fs_gps, datatypes)

    # mean the lat and lon and save as a file #
    cleaned_fnames = first.(splitext.(gps_pths)) .* "_cleaned.csv"
    mean_gps_positions(cleaned_fnames)
end


"""
    index_by_time(times, start_times, end_times)

Find the indices to separate `times` based on `start_times` and `end_times`.
"""
function index_by_time(times, start_times, end_times)

    f(x) = findmin(abs.(times .- x)) # define a function to minimize time difference

    # use the function and map to find the closest index #
    start_idx = last.(map(x -> f(x), start_times))
    end_idx = last.(map(x -> f(x), end_times))

    return range.(start_idx, end_idx)
end