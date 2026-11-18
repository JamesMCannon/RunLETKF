# make_minimum_inputs.jl
using Pkg
Pkg.activate(@__DIR__)

using JLD2, Parameters, Random, Dates, AxisKeys, StableRNGs, StatsBase
using LongwaveModePropagator, LMPTools
using ModifiedVLFInversionAlgorithms
const MVIA = ModifiedVLFInversionAlgorithms

include("common.jl")

"""
    slice_inputs(timeofday; keep=("NLK-CAL", "NML-PIN"))

Write `inputs/<timeofday>1_buildminimumpaths.jld2` from the `keep` rows of the AVID truth file.
Rows are assumed to follow `buildAVIDpaths()` order, as `observations` assumes.
"""
function slice_inputs(timeofday; keep=("NLK-CAL", "NML-PIN"))
    src = joinpath(@__DIR__, "inputs", timeofday * "1_buildAVIDpaths.jld2")
    dst = joinpath(@__DIR__, "inputs", timeofday * "1_buildminimumpaths.jld2")

    names = MVIA.pathname.(buildAVIDpaths())
    rows = [findfirst(==(k), names) for k in keep]
    MVIA.pathname.(buildminimumpaths()) == names[rows] ||
        error("buildminimumpaths() order does not match $keep")

    obsamp, obsphase, filekeys = jldopen(src, "r") do f
        f["obsamp"], f["obsphase"], keys(f)
    end
    @info "Source" src filekeys size(obsamp) size(obsphase)
    size(obsamp, 1) == length(names) ||
        error("expected $(length(names)) rows, got $(size(obsamp, 1))")

    jldsave(dst; obsamp=obsamp[rows, :], obsphase=obsphase[rows, :],
        source=basename(src), sourcerows=rows, paths=names[rows])
    @info "Wrote" dst rows
end

slice_inputs("day")