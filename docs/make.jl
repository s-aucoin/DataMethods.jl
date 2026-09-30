using Documenter
using DataMethods

makedocs(
    sitename = "DataMethods",
    format = Documenter.HTML(),
    modules = [DataMethods],
    remotes = nothing,
    pages = ["Home" => "index.md",
            "Miscellaneous Utilities" => "misc.md",
            "ADCP Utilities" => "ADCPUtilities/adcp_utilities.md",
            "CTD Utilities" => "CTDUtilities/ctd_utilities.md",
            "GPS Utilities" => "GPSUtilities/gps_utilities.md",
            "Library" => "library.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/s-aucoin/DataMethods.jl.git",
    versions = nothing
)
