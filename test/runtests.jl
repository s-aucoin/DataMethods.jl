using Test
using DataMethods

@testset "DataMethods.jl" begin
    include("DataMethods.jl")
    include("CalculusUtilities.jl")
    include("SignalProcessing.jl")
    include("adcp_utilities/ADCPUtilities.jl")
    include("adcp_utilities/SignalAnalysis.jl")
    include("adcp_utilities/WaveAnalysis.jl")
    include("adcp_utilities/transformations/CoordinateTransformations.jl")
    include("adcp_utilities/transformations/Transformations.jl")
    include("adcp_utilities/transformations/VelocityTransformations.jl")
    include("ctd_utilities/CTDUtilities.jl")
    include("gps_utilities/GPSUtilities.jl")
end
