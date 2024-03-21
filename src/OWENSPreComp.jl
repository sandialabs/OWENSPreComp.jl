module OWENSPreComp

export properties,tw_rate

# translated version of OWENSPreComp from NREL
include("main.jl")

# methods to read OWENSPreComp input files
include("io.jl")

include("deprecated.jl")

end # module
