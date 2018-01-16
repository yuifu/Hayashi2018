include("../../../170220_RSsite/bin/GetNormalizedCoverage.jl")

pathBamFiles = ARGS[1]
pathRegions = ARGS[2]
odirPrefix = ARGS[3]

mkpath(odirPrefix)

getNormalizedCoverage(pathBamFiles, pathRegions, odirPrefix)

