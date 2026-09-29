using Glob
using NaturalSort

export get_sig_fnames


include("transformations/Transformations.jl")
include("SignalAnalysis.jl")
include("WaveAnalysis.jl")


############################################
"""
    get_sig_fnames(filepath, file_prefix)

Get the ordered list of filepaths for Signature NetCDF file with given `file_prefix` in `filepath`.
"""
function get_sig_fnames(filepath, file_prefix)
    f_names = sort(glob(file_prefix * "*nc", filepath), lt=natural)

    # Need to reorder the data since the last file in chronological order is put first in the file list #
    raw_name = file_prefix * ".ad2cp"                      # given file name of the raw data (i.e. it will be the last file created)
    idx = findall( x -> occursin(raw_name, x), f_names) # find the index of all the nc files from the raw ad2cp file
    append!(f_names, f_names[idx])                      # add the indexed files to the end
    deleteat!(f_names, idx)                             # delete the indexed files from the list

    return f_names
end