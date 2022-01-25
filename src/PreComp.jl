module PreComp

export properties,tw_rate

# translated version of PreComp from NREL
include("main.jl")

# methods to read PreComp input files
include("io.jl")

include("deprecated.jl")

end # module
